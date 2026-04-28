library(tidyverse)
library(zoo)
library(sf)
library(data.table)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  merged_csv      <- snakemake@input[["merged_csv"]]
  out_dir         <- snakemake@output[["out_dir"]]
  time_steps      <- snakemake@config[["ncp_time_steps"]]
  basin_name      <- snakemake@params[["basin_name"]]
  basin_polygon   <- snakemake@params[["basin_polygon"]]
  date_start      <- snakemake@config[["date_start"]]
  date_end        <- snakemake@config[["date_end"]]
  zeu_default     <- snakemake@config[["zeu_default"]]
  mld_method      <- snakemake@config[["mld_method"]]
} else {
  args          <- commandArgs(trailingOnly = TRUE)
  merged_csv    <- ifelse(length(args) > 0, args[1],
                          "data/NorthAtlantic_20s/intermediate/merged/merged_ncp.csv")
  out_dir       <- ifelse(length(args) > 1, args[2],
                          "output/NorthAtlantic_20s/ncp/IcelandBasin")
  time_steps    <- c("10 days", "15 days", "20 days")
  basin_name    <- "IcelandBasin"
  basin_polygon <- list(c(-42,54), c(-15,54), c(-15,63), c(-22,63), c(-22,57), c(-42,57))
  date_start      <- "2020-01-01"
  date_end        <- "2026-01-01"
  zeu_default     <- 40
  mld_method      <- "average"
}

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Build sf polygon from vertex list ─────────────────────────────────────────
make_basin_polygon <- function(poly_verts) {
  m <- do.call(rbind, lapply(poly_verts, as.numeric))
  if (!all(m[1, ] == m[nrow(m), ])) m <- rbind(m, m[1, ])  # force close
  st_polygon(list(m)) |> st_sfc(crs = 4326)
}

basin_poly <- make_basin_polygon(basin_polygon)
bbox       <- st_bbox(basin_poly)

# Streaming read: only pull needed columns, then bbox-filter immediately so we
# never hold the full ~15GB merged CSV in memory.
needed_cols <- c("float_wmo","prof_number","date","lon","lat",
                 "MLD","depth","canyon_nitrate")
message("Reading merged data from ", merged_csv, " (cols: ",
        paste(needed_cols, collapse = ", "), ")")
dat_dt <- fread(merged_csv, select = needed_cols, showProgress = FALSE)
message("Loaded ", nrow(dat_dt), " rows; bbox-filtering...")
dat_dt <- dat_dt[lon >= bbox[["xmin"]] & lon <= bbox[["xmax"]] &
                 lat >= bbox[["ymin"]] & lat <= bbox[["ymax"]]]
message(nrow(dat_dt), " rows in basin bbox")

# Exact polygon filter (once, at load — avoids 3× repeated st_within inside
# compute_ncp).
if (nrow(dat_dt) > 0) {
  dat_sf   <- st_as_sf(dat_dt, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  in_basin <- lengths(st_within(dat_sf, basin_poly)) > 0
  dat_dt   <- dat_dt[in_basin]
  rm(dat_sf, in_basin); gc(verbose = FALSE)
}
dat <- as_tibble(dat_dt) |> mutate(date = as.Date(date))
rm(dat_dt); gc(verbose = FALSE)
message(nrow(dat), " rows in basin polygon")

# ── Integration helper ────────────────────────────────────────────────────────
integrate_nitrate <- function(nitrate, depth, from, to) {
  idx <- depth >= from & depth <= to & !is.na(nitrate)
  if (sum(idx) < 2) return(NA_real_)
  mean(approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)), na.rm = TRUE) * (to - from)
}

# ── Main NCP function ─────────────────────────────────────────────────────────
compute_ncp <- function(
    dat,
    time_step,
    basin_poly,
    date_start, date_end,
    zeu_default    = 40,
    mld_method     = "average",
    label          = NULL,
    mld_diag_path  = NULL
) {

  if (is.null(label)) label <- time_step

  # dat is already bbox- and polygon-filtered at load time; just apply date
  # window and add zeu here.
  dat <- dat |>
    filter(date > as.Date(date_start),
           date < as.Date(date_end)) |>
    mutate(zeu = zeu_default)

  if (nrow(dat) == 0) stop("No data after filtering for: ", label)
  message(label, ": ", nrow(dat), " rows after filtering")

  # per-profile MLD/Zeu
  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |>
    filter(!is.na(MLD), !is.na(zeu))

  # IQR outlier detection on per-profile MLD (1.5 × IQR rule)
  mld_q1  <- quantile(dat_prof$MLD, 0.25, na.rm = TRUE)
  mld_q3  <- quantile(dat_prof$MLD, 0.75, na.rm = TRUE)
  mld_iqr <- mld_q3 - mld_q1
  dat_prof <- dat_prof |>
    mutate(mld_outlier = MLD < (mld_q1 - 1.5 * mld_iqr) |
                         MLD > (mld_q3 + 1.5 * mld_iqr))

  n_outliers <- sum(dat_prof$mld_outlier)
  if (n_outliers > 0)
    message(label, ": removing ", n_outliers, " MLD outlier profile(s) ",
            "(IQR fence [", round(mld_q1 - 1.5 * mld_iqr, 1), ", ",
            round(mld_q3 + 1.5 * mld_iqr, 1), "] m)")

  dat_prof_clean <- dat_prof |> filter(!mld_outlier)

  # time grid
  time_grid <- tibble(
    date_grid = seq(min(dat_prof_clean$date), max(dat_prof_clean$date), by = time_step)
  )

  # integration depth — bin MLD per time step.
  # mld_method = "average":    always use bin mean MLD
  # mld_method = "winter_max": use bin max when mean MLD > 100 m (deep winter),
  #                            bin mean otherwise
  prof_dat_smoothed <- dat_prof_clean |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid) |>
    summarise(mld_mean = mean(MLD, na.rm = TRUE),
              mld_max  = max(MLD,  na.rm = TRUE),
              zeu      = mean(zeu, na.rm = TRUE),
              .groups  = "drop") |>
    mutate(mld = case_when(
      mld_method == "winter_max" & mld_mean > 100 ~ mld_max,
      TRUE ~ mld_mean
    )) |>
    mutate(
      prev_MLD               = lag(mld),
      prev_zeu               = lag(zeu),
      NCP_integration_depth  = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
      next_integration_depth = lead(NCP_integration_depth)
    ) |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    mutate(
      NCP_integration_depth  = zoo::na.approx(NCP_integration_depth,
                                               x       = as.numeric(date_grid),
                                               na.rm   = FALSE, rule = 2),
      next_integration_depth = zoo::na.approx(next_integration_depth,
                                               x       = as.numeric(date_grid),
                                               na.rm   = FALSE, rule = 2),
      dt = as.numeric(date_grid - lag(date_grid)),
      NCP_integration_depth_smooth  = rollapply(NCP_integration_depth,
                                                 width = 3, FUN = mean,
                                                 align = "center", fill = NA),
      next_integration_depth_smooth = rollapply(next_integration_depth,
                                                 width = 3, FUN = mean,
                                                 align = "center", fill = NA)
    ) |>
    filter(date_grid %in% time_grid$date_grid) |>
    na.omit()

  # ── MLD diagnostic plot (optional) ──────────────────────────────────────────
  if (!is.null(mld_diag_path)) {
    dat_prof_plot <- dat_prof |>
      mutate(date_grid = as.Date(cut(date, breaks = time_step)))
    smoothed_plot <- prof_dat_smoothed |> filter(!is.na(mld))
    p_mld <- ggplot() +
      geom_point(data = filter(dat_prof_plot, !mld_outlier),
                 aes(x = date, y = MLD),
                 color = "steelblue", alpha = 0.4, size = 1.2) +
      geom_point(data = filter(dat_prof_plot, mld_outlier),
                 aes(x = date, y = MLD),
                 color = "#e31a1c", shape = 4, size = 2.5, stroke = 1.2) +
      geom_line(data = smoothed_plot,
                aes(x = date_grid, y = mld),
                color = "black", linewidth = 1.2) +
      geom_point(data = smoothed_plot,
                 aes(x = date_grid, y = mld,
                     fill = if_else(mld_mean > 100, "max", "mean")),
                 color = "black", size = 2.5, shape = 21) +
      scale_fill_manual(values = c(max = "#ff7f00", mean = "white"),
                        name = "Bin MLD") +
      scale_y_reverse() +
      scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
      labs(x = "Date", y = "MLD (m)",
           title = paste0(basin_name, " — MLD (× = outlier removed, orange bin = winter-max used)")) +
      theme_bw() +
      theme(legend.position = "bottom")
    ggsave(mld_diag_path, p_mld, width = 11, height = 5, dpi = 150)
    message("MLD diagnostic plot saved -> ", mld_diag_path)
  }

  # nitrate interpolation
  dat_smoothed <- dat |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    group_by(depth) |>
    mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid),
                                    na.rm = FALSE, rule = 2)) |>
    ungroup() |>
    filter(date_grid %in% time_grid$date_grid) |>
    left_join(prof_dat_smoothed, by = "date_grid") |>
    na.omit()

  # nitrate integration
  dat_final <- dat_smoothed |>
    group_by(date_grid, mld, zeu) |>
    filter(!is.na(next_integration_depth)) |>
    summarise(
      int_N_mmol_m2      = integrate_nitrate(nitrate, depth, 0,
                                              unique(NCP_integration_depth)),
      next_int_N_mmol_m2 = integrate_nitrate(nitrate, depth, 0,
                                              unique(next_integration_depth)),
      .groups = "drop"
    )

  # NCP
  ncp_results <- dat_final |>
    mutate(
      int_N_smooth      = rollapply(int_N_mmol_m2,      width = 3, FUN = mean,
                                    align = "center", fill = NA),
      next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean,
                                    align = "center", fill = NA)
    ) |>
    arrange(date_grid) |>
    mutate(
      dt                  = as.numeric(date_grid - lag(date_grid)),
      nitrate_consumption = lag(next_int_N_smooth) - int_N_smooth,
      c_consumption       = nitrate_consumption * 6.625 / dt,
      NCP                 = c_consumption,
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP)
    ) |>
    filter(!is.na(NCP)) |>
    mutate(time_step_label = label)

  return(ncp_results)
}

# ── Run all time steps ────────────────────────────────────────────────────────
message("Running NCP for basin: ", basin_name)
message("Time steps: ", paste(time_steps, collapse = ", "))

all_results <- map(seq_along(time_steps), function(i) {
  ts <- time_steps[[i]]
  tryCatch(
    compute_ncp(
      dat            = dat,
      time_step      = ts,
      basin_poly     = basin_poly,
      date_start     = date_start, date_end = date_end,
      zeu_default    = zeu_default,
      mld_method     = mld_method,
      label          = ts,
      mld_diag_path  = if (i == 1L) file.path(out_dir, "mld_timeseries.png") else NULL
    ),
    error = function(e) { message("Failed for ", ts, ": ", e$message); NULL }
  )
}) |>
  compact() |>
  bind_rows()

# ── Save results ──────────────────────────────────────────────────────────────
results_csv <- file.path(out_dir, "ncp_results.csv")

if (nrow(all_results) == 0) {
  message("No float data found in basin '", basin_name, "' — writing empty outputs.")
  write_csv(
    tibble(date_grid = as.Date(character(0)), mld = numeric(0), zeu = numeric(0),
           NCP = numeric(0), time_step_label = character(0)),
    results_csv
  )
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("No data for basin: ", basin_name),
             size = 6, hjust = 0.5, vjust = 0.5) +
    theme_void()
  ggsave(file.path(out_dir, "ncp_sensitivity.png"), p_empty,
         width = 10, height = 5, dpi = 150)
  message("Empty outputs written for basin: ", basin_name)
} else {
  write_csv(all_results, results_csv)
  message("NCP results saved -> ", results_csv)

  p_sensitivity <- ggplot(all_results) +
    geom_line(aes(x = date_grid, y = NCP, color = time_step_label)) +
    scale_color_brewer(palette = "Set1", name = "Time step") +
    labs(y = "mmol C m-2 d-1", x = "Date",
         title = paste0("NCP sensitivity to time step -- ", basin_name)) +
    theme_bw() +
    scale_x_date(date_labels = "%b %Y", breaks = "6 months")

  ggsave(file.path(out_dir, "ncp_sensitivity.png"), p_sensitivity,
         width = 10, height = 5, dpi = 150)
  message("Sensitivity plot saved -> ", file.path(out_dir, "ncp_sensitivity.png"))
}
