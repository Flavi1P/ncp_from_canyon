library(tidyverse)
library(here)
library(seacarb)
library(mgcv)
library(data.table)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  in_dir       <- snakemake@input[["in_dir"]]
  out_dir      <- snakemake@output[["out_dir"]]
  utils_dir    <- file.path(dirname(snakemake@scriptdir), "scripts", "utils")
  canyon_script <- file.path(utils_dir, "fastr_canyon.R")
  if (!file.exists(canyon_script)) stop("Cannot find fastr_canyon.R at: ", canyon_script)
} else {
  args         <- commandArgs(trailingOnly = TRUE)
  in_dir       <- ifelse(length(args) > 0, args[1], "data/NorthAtlantic_2024/intermediate/doxy_profiles")
  out_dir      <- ifelse(length(args) > 1, args[2], "data/NorthAtlantic_2024/intermediate/nitrate_profiles")
  utils_dir    <- "scripts/utils"
}

# ── Load CANYON-B ─────────────────────────────────────────────────────────────
source(file.path(utils_dir, "fastr_canyon.R"))
weights <- load_canyonb_weights(inputsdir = paste0(utils_dir, "/canyon_weights/"))

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "_interp\\.csv$", full.names = TRUE)
message("Found ", length(files), " interpolated float files in ", in_dir)

# ── Process each float ────────────────────────────────────────────────────────
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
i  <- 0
processed <- c()

for (file in files) {

  i <- i + 1
  setTxtProgressBar(pb, i)

  dat <- read_csv(file, show_col_types = FALSE)

  # skip if missing required columns
  required_cols <- c("date", "lat", "lon", "depth", "temp", "sal", "oxygen")
  if (!all(required_cols %in% names(dat))) {
    message("Skipping ", basename(file), " — missing required columns.")
    next
  }

  # skip if too few valid oxygen rows to bother predicting
  if (sum(!is.na(dat$oxygen)) < 10) {
    message("Skipping ", basename(file), " — insufficient DOXY data.")
    next
  }

  # predict nitrate with CANYON-B
  tryCatch({
    nitrate <- CANYONB_fast(
      date         = dat$date,
      lat          = dat$lat,
      lon          = dat$lon,
      pres         = dat$depth,
      temp         = dat$temp,
      psal         = dat$sal,
      doxy         = dat$oxygen,
      wgts_list    = weights,
      use_parallel = FALSE,
      param        = "NO3"
    )

    dat$canyon_nitrate <- nitrate$NO3

  }, error = function(e) {
    message("CANYON-B failed for ", basename(file), ": ", e$message)
    dat$canyon_nitrate <<- NA_real_
  })

  # save
  wmo      <- gsub("argo_(.*)_interp\\.csv", "\\1", basename(file))
  out_path <- file.path(out_dir, paste0("argo_", wmo, "_nitrate.csv"))
  write_csv(dat, out_path)
  processed <- c(processed, out_path)
}

close(pb)
message("Done. Processed ", length(processed), " floats.")

# ── Write manifest ────────────────────────────────────────────────────────────
manifest <- data.frame(
  wmo  = gsub("argo_(.*)_nitrate\\.csv", "\\1", basename(processed)),
  path = processed
)
write_csv(manifest, file.path(out_dir, "nitrate_manifest.csv"))
message("Manifest → ", file.path(out_dir, "nitrate_manifest.csv"))