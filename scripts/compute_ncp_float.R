library(tidyverse)
library(zoo)
library(data.table)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  merged_csv        <- snakemake@input[["merged_csv"]]
  results_csv       <- snakemake@output[["results_csv"]]
  plot_png          <- snakemake@output[["plot_png"]]
  transect_png      <- snakemake@output[["transect_png"]]
  target_wmo        <- as.character(snakemake@params[["wmo"]])
  zeu_default       <- snakemake@params[["zeu_default"]]
  float_time_window <- snakemake@params[["float_time_window"]]  # NULL = per-profile
} else {
  args              <- commandArgs(trailingOnly = TRUE)
  target_wmo        <- ifelse(length(args) > 0, args[1], "6904240")
  merged_csv        <- ifelse(length(args) > 1, args[2],
                              "data/NorthAtlantic_20s/intermediate/merged/merged_ncp.csv")
  results_csv       <- paste0("output/NorthAtlantic_20s/ncp_float/", target_wmo, "/ncp_float.csv")
  plot_png          <- paste0("output/NorthAtlantic_20s/ncp_float/", target_wmo, "/ncp_float.png")
  transect_png      <- paste0("output/NorthAtlantic_20s/ncp_float/", target_wmo, "/nitrate_transect.png")
  zeu_default       <- 40
  float_time_window <- "15 days"   # set NULL for per-profile mode
}

use_time_window <- !is.null(float_time_window)
dir.create(dirname(results_csv), recursive = TRUE, showWarnings = FALSE)

# ── Integration helper ────────────────────────────────────────────────────────
integrate_nitrate <- function(nitrate, depth, from, to) {
  idx <- depth >= from & depth <= to & !is.na(nitrate)
  if (sum(idx) < 2) return(NA_real_)
  mean(approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)), na.rm = TRUE) * (to - from)
}

# ── Load and filter data ──────────────────────────────────────────────────────
# Streaming read: only pull needed columns, filter by WMO immediately so we
# never hold the full ~15GB merged CSV in memory.
needed_cols <- c("float_wmo","prof_number","date","lon","lat",
                 "MLD","depth","canyon_nitrate")
message("Reading merged data from ", merged_csv, " (filter wmo=", target_wmo, ")")
dat_dt <- fread(merged_csv, select = needed_cols, showProgress = FALSE)
dat_dt <- dat_dt[as.character(float_wmo) == target_wmo]
dat <- as_tibble(dat_dt) |>
  mutate(date = as.Date(date), zeu = zeu_default)
rm(dat_dt); gc(verbose = FALSE)

if (nrow(dat) == 0) stop("No data for WMO: ", target_wmo)

n_valid_nit <- sum(!is.na(dat$canyon_nitrate))
if (n_valid_nit == 0)
  stop("WMO ", target_wmo, ": all canyon_nitrate values are NA — ",
       "CANYON-B likely failed for this float (check lon/lat in the doxy profiles)")

message("WMO ", target_wmo, ": ", nrow(dat), " depth rows across ",
        n_distinct(dat$prof_number), " profiles")

# ── Per-profile MLD (shared by both modes) ───────────────────────────────────
dat_prof <- dat |>
  select(prof_number, date, MLD, zeu) |>
  distinct() |>
  filter(!is.na(MLD)) |>
  arrange(date)

mld_q1  <- quantile(dat_prof$MLD, 0.25, na.rm = TRUE)
mld_q3  <- quantile(dat_prof$MLD, 0.75, na.rm = TRUE)
mld_iqr <- mld_q3 - mld_q1
dat_prof <- dat_prof |>
  mutate(mld_outlier = MLD < (mld_q1 - 1.5 * mld_iqr) |
                       MLD > (mld_q3 + 1.5 * mld_iqr))

n_outliers <- sum(dat_prof$mld_outlier)
if (n_outliers > 0)
  message("WMO ", target_wmo, ": removing ", n_outliers, " MLD outlier profile(s)")

dat_prof_clean <- dat_prof |> filter(!mld_outlier)
dat_clean      <- dat |> filter(prof_number %in% dat_prof_clean$prof_number)

# ── Transect helper: tiles sized to their actual time span ───────────────────
transect_tile_width <- function(dates) {
  if (length(unique(dates)) < 2) return(14)
  as.numeric(median(diff(sort(unique(dates)))))
}

# ── Branch: windowed vs per-profile ──────────────────────────────────────────
if (use_time_window) {

  message("Mode: time window = ", float_time_window)

  # Time grid covering the clean profiles
  time_grid <- tibble(
    date_grid = seq(min(dat_prof_clean$date),
                    max(dat_prof_clean$date),
                    by = float_time_window)
  )

  # Bin-average MLD per window → integration depths
  prof_smoothed <- dat_prof_clean |>
    mutate(date_grid = as.Date(cut(date, breaks = float_time_window))) |>
    group_by(date_grid) |>
    summarise(mld = mean(MLD, na.rm = TRUE),
              zeu = mean(zeu, na.rm = TRUE),
              .groups = "drop") |>
    mutate(
      prev_MLD       = lag(mld),
      prev_zeu       = lag(zeu),
      NCP_depth      = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
      next_NCP_depth = lead(NCP_depth)
    ) |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    mutate(
      NCP_depth      = zoo::na.approx(NCP_depth,      x = as.numeric(date_grid),
                                       na.rm = FALSE, rule = 2),
      next_NCP_depth = zoo::na.approx(next_NCP_depth, x = as.numeric(date_grid),
                                       na.rm = FALSE, rule = 2),
      dt             = as.numeric(date_grid - lag(date_grid))
    ) |>
    filter(date_grid %in% time_grid$date_grid) |>
    na.omit()

  # Running mean of nitrate at each depth (bin average within each window)
  nitrate_binned <- dat_clean |>
    mutate(date_grid = as.Date(cut(date, breaks = float_time_window))) |>
    group_by(date_grid, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    group_by(depth) |>
    mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid),
                                     na.rm = FALSE, rule = 2)) |>
    ungroup() |>
    filter(date_grid %in% time_grid$date_grid)

  # ── Transect plot ─────────────────────────────────────────────────────────
  transect_data <- nitrate_binned |>
    filter(!is.na(nitrate), depth <= 500)

  bin_w   <- transect_tile_width(transect_data$date_grid)
  mld_pts <- prof_smoothed |> filter(!is.na(mld))

  p_transect <- ggplot(transect_data,
                       aes(x = date_grid, y = depth, fill = nitrate)) +
    geom_tile(width = bin_w) +
    geom_line(data = mld_pts, aes(x = date_grid, y = mld),
              color = "white", linewidth = 1.2, inherit.aes = FALSE) +
    scale_y_reverse() +
    scale_fill_viridis_c(option = "plasma",
                         name   = "NO₃\n(mmol m⁻³)",
                         na.value = "grey85") +
    scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
    labs(title    = paste0("Float ", target_wmo,
                           " — nitrate transect (", float_time_window, " bins)"),
         x = "Date", y = "Depth (m)") +
    theme_bw()

  # ── NCP ───────────────────────────────────────────────────────────────────
  dat_smoothed <- nitrate_binned |>
    left_join(prof_smoothed |>
                select(date_grid, mld, zeu, NCP_depth, next_NCP_depth, dt),
              by = "date_grid") |>
    na.omit()

  dat_final <- dat_smoothed |>
    group_by(date_grid, dt) |>
    filter(!is.na(next_NCP_depth)) |>
    summarise(
      int_N      = integrate_nitrate(nitrate, depth, 0, unique(NCP_depth)),
      next_int_N = integrate_nitrate(nitrate, depth, 0, unique(next_NCP_depth)),
      .groups    = "drop"
    )

  ncp_results <- dat_final |>
    arrange(date_grid) |>
    mutate(
      nitrate_consumption = lag(next_int_N) - int_N,
      NCP                 = nitrate_consumption * 6.625 / dt,
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP)
    ) |>
    filter(!is.na(NCP)) |>
    transmute(float_wmo = target_wmo, date = date_grid,
              prof_number = NA_integer_, dt, NCP)

} else {

  message("Mode: per-profile (consecutive profiles)")

  # ── Transect plot (raw profiles) ──────────────────────────────────────────
  transect_data <- dat_clean |>
    group_by(date, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(nitrate), depth <= 500)

  bin_w   <- transect_tile_width(transect_data$date)
  mld_pts <- dat_prof_clean

  p_transect <- ggplot(transect_data,
                       aes(x = date, y = depth, fill = nitrate)) +
    geom_tile(width = bin_w) +
    geom_line(data = mld_pts, aes(x = date, y = MLD),
              color = "white", linewidth = 1.2, inherit.aes = FALSE) +
    scale_y_reverse() +
    scale_fill_viridis_c(option = "plasma",
                         name   = "NO₃\n(mmol m⁻³)",
                         na.value = "grey85") +
    scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
    labs(title    = paste0("Float ", target_wmo,
                           " — nitrate transect (per-profile)"),
         x = "Date", y = "Depth (m)") +
    theme_bw()

  # ── Per-profile NCP ───────────────────────────────────────────────────────
  dat_prof_clean <- dat_prof_clean |>
    mutate(
      prev_MLD       = lag(MLD),
      prev_zeu       = lag(zeu),
      NCP_depth      = pmax(MLD, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
      next_NCP_depth = lead(NCP_depth),
      dt             = as.numeric(date - lag(date))
    )

  dat_final <- dat_clean |>
    left_join(dat_prof_clean |>
                select(prof_number, NCP_depth, next_NCP_depth, dt),
              by = "prof_number") |>
    group_by(prof_number, date, dt) |>
    filter(!is.na(next_NCP_depth)) |>
    summarise(
      int_N      = integrate_nitrate(canyon_nitrate, depth, 0, unique(NCP_depth)),
      next_int_N = integrate_nitrate(canyon_nitrate, depth, 0, unique(next_NCP_depth)),
      .groups    = "drop"
    )

  ncp_results <- dat_final |>
    arrange(date) |>
    mutate(
      nitrate_consumption = lag(next_int_N) - int_N,
      NCP                 = nitrate_consumption * 6.625 / dt,
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP)
    ) |>
    filter(!is.na(NCP)) |>
    transmute(float_wmo = target_wmo, date, prof_number, dt, NCP)
}

message("WMO ", target_wmo, ": ", nrow(ncp_results), " NCP estimates computed")

# ── Save outputs ──────────────────────────────────────────────────────────────
ggsave(transect_png, p_transect, width = 12, height = 5, dpi = 150)
message("Transect plot saved -> ", transect_png)

write_csv(ncp_results, results_csv)
message("Results saved -> ", results_csv)

mode_label <- if (use_time_window) paste0(float_time_window, " window") else "per-profile"

if (nrow(ncp_results) > 0) {
  p_ncp <- ggplot(ncp_results, aes(x = date, y = NCP)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 1.8) +
    labs(title = paste0("NCP — float ", target_wmo, "  (", mode_label, ")"),
         y = "NCP  (mmol C m⁻² d⁻¹)", x = "Date") +
    theme_bw() +
    scale_x_date(date_labels = "%b %Y", breaks = "3 months")
} else {
  p_ncp <- ggplot() +
    labs(title = paste0("No NCP estimates — float ", target_wmo)) +
    theme_bw()
}

ggsave(plot_png, p_ncp, width = 10, height = 4, dpi = 150)
message("NCP plot saved -> ", plot_png)
