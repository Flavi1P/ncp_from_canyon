library(tidyverse)
library(zoo)
library(here)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  merged_csv  <- snakemake@input[["merged_csv"]]
  out_dir     <- snakemake@output[["out_dir"]]
  time_steps  <- snakemake@config[["ncp_time_steps"]]
  lon_min     <- snakemake@config[["lon_min"]]
  lon_max     <- snakemake@config[["lon_max"]]
  lat_min     <- snakemake@config[["lat_min"]]
  lat_max     <- snakemake@config[["lat_max"]]
  date_start  <- snakemake@config[["date_start"]]
  date_end    <- snakemake@config[["date_end"]]
  mld_spar    <- snakemake@config[["mld_spar"]]
  zeu_default <- snakemake@config[["zeu_default"]]
} else {
  args        <- commandArgs(trailingOnly = TRUE)
  merged_csv  <- ifelse(length(args) > 0, args[1],
                        "data/NorthAtlantic_2024/intermediate/merged/merged_ncp.csv")
  out_dir     <- ifelse(length(args) > 1, args[2],
                        "data/NorthAtlantic_2024/ncp")
  time_steps  <- c("10 days", "15 days", "30 days")
  lon_min     <- -40;  lon_max    <- -10
  lat_min     <-  58;  lat_max    <-  65
  date_start  <- "2023-01-01"
  date_end    <- "2024-01-01"
  mld_spar    <- 0.3
  zeu_default <- 40
}

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading merged data from ", merged_csv)
dat <- read_csv(merged_csv, show_col_types = FALSE)
message("Loaded ", nrow(dat), " rows")

# ── Integration helper ────────────────────────────────────────────────────────
integrate_nitrate <- function(nitrate, depth, from, to) {
  idx <- depth >= from & depth <= to
  if (sum(idx) < 2) return(NA_real_)
  mean(approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)), na.rm = TRUE) * (to - from)
}

# ── Main NCP function ─────────────────────────────────────────────────────────
compute_ncp <- function(
    dat,
    time_step,
    lon_min, lon_max, lat_min, lat_max,
    date_start, date_end,
    zeu_default = 40,
    mld_spar    = 0.3,
    label       = NULL
) {

  if (is.null(label)) label <- time_step

  # filter
  dat <- dat |>
    filter(
      lon  >= lon_min & lon  <= lon_max,
      lat  >= lat_min & lat  <= lat_max,
      date >  as.Date(date_start),
      date <  as.Date(date_end)
    ) |>
    mutate(zeu = zeu_default)

  if (nrow(dat) == 0) stop("No data after filtering for: ", label)
  message(label, ": ", nrow(dat), " rows after filtering")

  # per-profile MLD/Zeu
  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |>
    arrange(date) |>
    filter(!is.na(MLD), !is.na(zeu))

  # smooth MLD
  dat_smooth <- dat_prof |> filter(!is.na(date), !is.na(MLD))
  fit_mld    <- loess(MLD ~ as.numeric(date), data = dat_smooth,
                      span = mld_spar, family = "symmetric")
  fit_zeu    <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)

  dat_prof <- dat_prof |>
    mutate(
      mld_smooth = pmax(predict(fit_mld,
                                newdata = data.frame(date = as.numeric(dat_prof$date))), 0),
      zeu_smooth = predict(fit_zeu, as.numeric(dat_prof$date))$y
    )

  dat <- left_join(dat,
                   select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth),
                   by = c("float_wmo", "prof_number"))

  # time grid
  time_grid <- tibble(
    date_grid = seq(min(dat_prof$date), max(dat_prof$date), by = time_step)
  )

  # integration depth
  prof_dat_smoothed <- dat |>
    group_by(float_wmo, prof_number, date) |>
    summarise(mld_smooth = unique(mld_smooth),
              zeu_smooth = unique(zeu_smooth),
              .groups    = "drop") |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid) |>
    summarise(mld = mean(mld_smooth, na.rm = TRUE),
              zeu = mean(zeu_smooth, na.rm = TRUE),
              .groups = "drop") |>
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
      int_N_mmol_m2         = integrate_nitrate(nitrate, depth, 0,
                                                 unique(NCP_integration_depth)),
      next_int_N_mmol_m2    = integrate_nitrate(nitrate, depth, 0,
                                                 unique(next_integration_depth)),
      mld_concentration     = integrate_nitrate(nitrate, depth,
                                                 unique(NCP_integration_depth) - 20,
                                                 unique(NCP_integration_depth) - 10) / 10,
      sub_mld_concentration = integrate_nitrate(nitrate, depth,
                                                 unique(NCP_integration_depth) + 20,
                                                 unique(NCP_integration_depth) + 30) / 10,
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
      diff_mld            = mld - lag(mld),
      we                  = pmax(0, diff_mld),
      delta               = sub_mld_concentration - mld_concentration,
      diff                = delta * we,
      nitrate_consumption = lag(next_int_N_smooth) - int_N_smooth,
      c_consumption       = nitrate_consumption * 6.625 / dt,
      NCP                 = c_consumption + diff,
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP),
      datenum             = as.numeric(date_grid)
    ) |>
    filter(!is.na(NCP))

  ncp_results$time_step_label <- label

  return(ncp_results)
}

# ── Run all time steps ────────────────────────────────────────────────────────
message("Running NCP for time steps: ", paste(time_steps, collapse = ", "))

all_results <- map(time_steps, function(ts) {
  tryCatch(
    compute_ncp(
      dat        = dat,
      time_step  = ts,
      lon_min    = lon_min, lon_max = lon_max,
      lat_min    = lat_min, lat_max = lat_max,
      date_start = date_start, date_end = date_end,
      mld_spar   = mld_spar,
      zeu_default= zeu_default,
      label      = ts
    ),
    error = function(e) { message("Failed for ", ts, ": ", e$message); NULL }
  )
}) |>
  compact() |>
  bind_rows()

# ── Save results ──────────────────────────────────────────────────────────────
results_csv <- file.path(out_dir, "ncp_results.csv")
write_csv(all_results, results_csv)
message("NCP results saved → ", results_csv)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_sensitivity <- ggplot(all_results) +
  geom_line(aes(x = date_grid, y = NCP,       color = time_step_label), alpha = 0.3) +
  scale_color_brewer(palette = "Set1", name = "Time step") +
  labs(y = "mmol C m⁻² d⁻¹", x = "Date",
       title = "NCP sensitivity to time step") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

ggsave(file.path(out_dir, "ncp_sensitivity.png"), p_sensitivity,
       width = 10, height = 5, dpi = 150)
message("Sensitivity plot saved → ", file.path(out_dir, "ncp_sensitivity.png"))