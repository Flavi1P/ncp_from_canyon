# Optimized CANYON-B implementation (vectorized + optional parallel + preloaded weights)
# - Requires mgcv and seacarb to be loaded outside the function
# - Usage: preload weights with load_canyonb_weights() and pass wgts_list for best performance

load_canyonb_weights <- function(inputsdir = "") {
  paramnames <- c('AT','CT','pH','pCO2','NO3','PO4','SiOH4')
  wlist <- setNames(lapply(paramnames, function(p) {
    read.table(paste0(inputsdir, "wgts_", p, ".txt"), stringsAsFactors = FALSE)
  }), paramnames)
  return(wlist)
}

CANYONB_fast <- function(date, lat, lon, pres, temp, psal, doxy,
                         param = NULL,
                         epres = 0.5, etemp = 0.005, epsal = 0.005, edoxy = NULL,
                         inputsdir = "", wgts_list = NULL, use_parallel = FALSE, mc.cores = NULL) {
  # NOTE: make sure library(mgcv); library(seacarb) are called outside this function
  if (is.null(edoxy)) edoxy <- 0.01 * doxy
  nol <- length(pres)
  if (length(epres) == 1) epres <- rep(epres, nol)
  if (length(etemp) == 1) etemp <- rep(etemp, nol)
  if (length(epsal) == 1) epsal <- rep(epsal, nol)
  if (length(edoxy) == 1) edoxy <- rep(edoxy, nol)
  
  paramnames <- c('AT','CT','pH','pCO2','NO3','PO4','SiOH4')
  noparams <- length(paramnames)
  inputsigma <- c(6,4,.005,NaN,2/100,2/100,2/100)
  betaipCO2 <- c(-3.114e-05,1.087e-01,-7.899e+01)
  inputsigma[3] <- sqrt(.005^2 + .01^2)
  
  if (is.null(param)) paramflag <- rep(TRUE, noparams) else paramflag <- paramnames %in% param
  if (paramflag[4] && !requireNamespace("seacarb", quietly = TRUE)) stop("seacarb required for pCO2 calculations")
  
  # preload weights if not provided
  if (is.null(wgts_list)) wgts_list <- load_canyonb_weights(inputsdir)
  
  # normalize lon range
  lon[lon > 180] <- lon[lon > 180] - 360
  date <- as.POSIXlt(date, tz = "UTC")
  year <- as.numeric(format(date, "%Y")) + as.numeric(difftime(date, as.POSIXlt(paste0(format(date, "%Y"), "-01-01 00:00"), tz = "UTC"), units = "days"))/365
  
  # polar shift (vectorized in/out from mgcv::in.out used as-is)
  plon <- c(-180,-170,-85,-80,-37,-37,143,143,180,180,-180,-180)
  plat <- c(68,66.5,66.5,80,80,90,90,68,68,90,90,68)
  arcflag <- mgcv::in.out(cbind(plat,plon), cbind(lat, lon))
  lat[arcflag] <- lat[arcflag] - sinpi((lon[arcflag] + 37)/180) * (90 - lat[arcflag]) * 0.5
  
  # construct input matrix (year, lat/90, sin/cos lon-like transforms, temp, sal, oxy, pressure transform)
  data <- cbind(year,
                lat/90,
                abs(1 - ((lon - 110) %% 360)/180),
                abs(1 - ((lon - 20) %% 360)/180),
                temp, psal, doxy,
                pres/2e4 + 1/((1 + exp(-pres/300))^3))
  
  out <- list()
  
  # helper to compute one network's outputs (vectorized across samples)
  compute_network <- function(inwgts, data, ni) {
    # inwgts: data frame of weights for one parameter
    noparsets <- ncol(inwgts) - 1
    # network results containers
    cvals <- matrix(NA_real_, nrow = nol, ncol = noparsets)
    cvalcys <- numeric(noparsets)
    inval_list <- vector("list", noparsets)
    
    for (l in seq_len(noparsets)) {
      nl1 <- as.integer(inwgts[1, l])
      nl2 <- as.integer(inwgts[2, l])
      beta <- as.numeric(inwgts[3, l])
      # unpack weights (follow original indexing)
      # w1: nl1 x ni
      off <- 4
      w1 <- matrix(as.numeric(inwgts[off + (1:(nl1 * ni)), l]), nrow = nl1, ncol = ni)
      off <- off + nl1 * ni
      b1 <- as.numeric(inwgts[off + (1:nl1), l]); off <- off + nl1
      w2 <- matrix(as.numeric(inwgts[off + (1:(nl2 * nl1)), l]), nrow = nl2, ncol = nl1); off <- off + nl2 * nl1
      b2 <- as.numeric(inwgts[off + (1:nl2), l]); off <- off + nl2
      
      nlayerflag <- 1 + as.numeric(nl2 != 0)
      
      if (nlayerflag == 1) {
        # one hidden layer
        # data_N already provided outside
        a <- data %*% t(w1) + matrix(rep(b1, each = nol), nrow = nol)
        y <- tanh(a) %*% t(w2) + matrix(rep(b2, each = nol), nrow = nol)
        cvals[, l] <- as.numeric(y)
        cvalcys[l] <- 1 / beta
        
        # input effect: inx = (1 - tanh(a)^2) %*% ( (w2) * w1 )
        tanhprime <- 1 - tanh(a)^2                      # nol x nl1
        M <- w2 * w1                                   # nl1 x ni  (w2 is 1 x nl1 for no=1)
        inx <- tanhprime %*% M                         # nol x ni
        inval_list[[l]] <- inx
        
      } else {
        # two hidden layers
        # redistribute correctly; note nl2 > 0
        w3 <- matrix(as.numeric(inwgts[off + (1:(1 * nl2)), l]), nrow = 1, ncol = nl2); off <- off + nl2
        b3 <- as.numeric(inwgts[off + (1:1), l]) # not used for derivative
        
        a <- data %*% t(w1) + matrix(rep(b1, each = nol), nrow = nol)
        b <- tanh(a) %*% t(w2) + matrix(rep(b2, each = nol), nrow = nol)
        y <- tanh(b) %*% t(w3) + matrix(rep(b3, each = nol), nrow = nol)
        cvals[, l] <- as.numeric(y)
        cvalcys[l] <- 1 / beta
        
        # input effect using chain rule vectorized:
        D1 <- 1 - tanh(a)^2  # nol x nl1
        D2 <- 1 - tanh(b)^2  # nol x nl2
        # (D2 * w3) %*% w2  -> nol x nl1
        U <- (D2 %*% diag(as.numeric(w3))) %*% w2   # nol × nl1
        V <- U * D1                                 # nol × nl1
        inx <- V %*% w1                     # nol x ni
        inval_list[[l]] <- inx
      }
    }
    list(cvals = cvals, cvalcys = cvalcys, inval_list = inval_list)
  }
  
  # loop over parameters
  for (i in seq_len(noparams)) {
    if (!paramflag[i]) next
    inwgts <- wgts_list[[paramnames[i]]]
    
    # determine ni and normalization vectors from weights (last column stores mw/sw)
    lastcol <- ncol(inwgts)
    # number of inputs in normalization stored depends on param type
    if (i > 4) {
      ni <- ncol(data) - 1
      ioffset <- -1
      mw <- as.numeric(inwgts[1:(ni + 1), lastcol])
      sw <- as.numeric(inwgts[(ni + 2):(2 * ni + 2), lastcol])
      # data_N excludes year
      data_inputs <- data[, 2:(ni + 1), drop = FALSE]
      data_N <- sweep(data_inputs, 2, mw[1:ni], "-")
      data_N <- sweep(data_N, 2, sw[1:ni], "/")
    } else {
      ni <- ncol(data)
      ioffset <- 0
      mw <- as.numeric(inwgts[1:(ni + 1), lastcol])
      sw <- as.numeric(inwgts[(ni + 2):(2 * ni + 2), lastcol])
      data_N <- sweep(data, 2, mw[1:ni], "-")
      data_N <- sweep(data_N, 2, sw[1:ni], "/")
    }
    
    # prepare wgts and betas
    noparsets <- ncol(inwgts) - 1
    wgts <- as.numeric(unlist(inwgts[4, 1:noparsets], use.names = FALSE))
    betaciw <- as.numeric(inwgts[(2 * ni + 3):nrow(inwgts), noparsets + 1])
    betaciw <- betaciw[!is.nan(betaciw)]
    
    # compute networks (optionally in parallel)
    if (use_parallel && noparsets > 1) {
      if (is.null(mc.cores)) mc.cores <- parallel::detectCores(logical = FALSE)
      res <- parallel::mclapply(seq_len(noparsets), function(l) {
        # reuse compute logic by extracting single-network table
        # construct a one-column table for compute_network
        tbl <- inwgts[, c(l, ncol(inwgts)), drop = FALSE]
        compute_network(tbl, data_N, ni)
      }, mc.cores = mc.cores)
      # combine results
      cvals_mat <- do.call(cbind, lapply(res, function(z) z$cvals))
      cvalcys_vec <- unlist(lapply(res, function(z) z$cvalcys))
      inval_all <- lapply(res, function(z) z$inval_list)
      # flatten inval_all to list of matrices and then to array
      inval_list <- do.call(c, lapply(inval_all, identity))
    } else {
      res_single <- compute_network(inwgts, data_N, ni)
      cvals_mat <- res_single$cvals
      cvalcys_vec <- res_single$cvalcys
      inval_list <- res_single$inval_list
    }
    
    # denormalize cval (output)
    cval <- sweep(cvals_mat, 2, sw[ni + 1], FUN = "*")
    cval <- sweep(cval, 2, mw[ni + 1], FUN = "+")
    
    # combine committee
    V1 <- sum(wgts); V2 <- sum(wgts^2)
    # weighted mean across committee (columns are networks)
    out_val <- rowSums(t(t(cval) * wgts)) / V1
    out[[paramnames[i]]] <- out_val
    
    # committee uncertainty components
    # CU variance
    cvalcu <- rowSums(t(t(cval) * wgts) * (cval - matrix(out_val, nrow = nol, ncol = noparsets))^2) / (V1 - V2 / V1)
    # noise variance averaged
    cvalcib <- rep(sum(wgts * cvalcys_vec) / V1, nol)
    # weight uncertainty via parameterization
    cvalciw <- (betaciw[2] + betaciw[1] * sqrt(cvalcu))^2
    
    # combine input effects: first create 3D array (nol x ni x noparsets) for inval_list
    inval_array <- array(NA_real_, dim = c(nol, ni, length(inval_list)))
    for (ll in seq_along(inval_list)) inval_array[,,ll] <- inval_list[[ll]]
    # weighted mean across networks
    inx <- rowSums(aperm(array(wgts, c(length(inval_list), nol, ni)), c(2,3,1)) * inval_array, dims = 2) / V1
    
    # rescale for normalization inside the MLP
    inx <- sweep(inx, 2, sw[ni + 1] / sw[1:ni], FUN = "*")
    # additional pressure scaling
    ddp <- 1/2e4 + 1/((1 + exp(-pres/300))^4) * exp(-pres/300) / 100
    inx[, 8 + ioffset] <- inx[, 8 + ioffset] * ddp
    cvalcin <- rowSums(inx[,(5:8) + ioffset]^2 * cbind(etemp, epsal, edoxy, epres)^2)
    
    # measurement/reference uncertainty
    if (i > 4) {
      cvalcimeas <- (inputsigma[i] * out[[paramnames[i]]])^2
    } else if (i == 4) {
      cvalcimeas <- (betaipCO2[1] * out[[paramnames[i]]]^2 + betaipCO2[2] * out[[paramnames[i]]] + betaipCO2[3])^2
    } else {
      cvalcimeas <- inputsigma[i]^2
    }
    
    # final uncertainties
    out[[paste0(paramnames[i], "_ci")]] <- sqrt(cvalcimeas + cvalcib + cvalciw + cvalcu + cvalcin)
    out[[paste0(paramnames[i], "_cim")]] <- sqrt(cvalcimeas)
    out[[paste0(paramnames[i], "_cin")]] <- sqrt(cvalcib + cvalciw + cvalcu)
    out[[paste0(paramnames[i], "_cii")]] <- sqrt(cvalcin)
    
    # pCO2: convert if needed
    if (i == 4) {
      carb_res <- seacarb::carb(flag = 15, var1 = 0.0023, var2 = out[[paramnames[i]]] * 1e-6,
                                S = 35, T = 25, P = 0, Patm = 1, Pt = 0, Sit = 0, pHscale = "T",
                                kf = "pf", k1k2 = "l", ks = "d", b = "u74", gas = "standard")
      deriv <- seacarb::derivnum('var2', flag = 15, var1 = 0.0023, var2 = out[[paramnames[i]]] * 1e-6,
                                 S = 35, T = 25, P = 0, Patm = 1, Pt = 0, Sit = 0, pHscale = "T",
                                 kf = "pf", k1k2 = "l", ks = "d", b = "u74", gas = "standard")
      out[[paramnames[i]]] <- carb_res$pCO2
      # scale uncertainties via derivative
      scale_fac <- deriv$pCO2 * 1e-6
      out[[paste0(paramnames[i], "_ci")]] <- scale_fac * out[[paste0(paramnames[i], "_ci")]]
      out[[paste0(paramnames[i], "_cim")]] <- scale_fac * out[[paste0(paramnames[i], "_cim")]]
      out[[paste0(paramnames[i], "_cin")]] <- scale_fac * out[[paste0(paramnames[i], "_cin")]]
      out[[paste0(paramnames[i], "_cii")]] <- scale_fac * out[[paste0(paramnames[i], "_cii")]]
    }
    
    # cleanup per-parameter
    rm(inwgts, mw, sw, cvals_mat, cvalcys_vec, inval_list)
  }
  
  return(out)
}

# Example usage:
# library(mgcv); library(seacarb)
# weights <- load_canyonb_weights(inputsdir = "path/to/wgts/")
# res <- CANYONB_fast(date="2014-12-09 08:45", lat=17.6, lon=-24.3, pres=180, temp=16, psal=36.1, doxy=104,
#                     wgts_list = weights, use_parallel = FALSE)
