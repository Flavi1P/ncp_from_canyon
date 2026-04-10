suppressPackageStartupMessages({library(tidyverse); library(zoo)})

merged_csv  <- "data/NorthAtlantic_test/intermediate/merged/merged_ncp.csv"
out_dir     <- "output/NorthAtlantic_test/ncp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lon_min <- -45; lon_max <- -10; lat_min <- 54; lat_max <- 63
date_start <- "2023-01-01"; date_end <- "2024-01-01"
mld_spar <- 0.3; zeu_default <- 40

dat <- read_csv(merged_csv, show_col_types = FALSE)
message("Loaded ", nrow(dat), " rows")

# source compute_ncp helpers without running the top-level code
# (re-define the two helpers inline to avoid sourcing the full script)
integrate_nitrate <- function(nitrate, depth, from, to) {
  idx <- depth >= from & depth <= to
  if (sum(idx) < 2) return(NA_real_)
  mean(approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)), na.rm = TRUE) * (to - from)
}

compute_ncp <- function(dat, time_step, lon_min, lon_max, lat_min, lat_max,
                        date_start, date_end, zeu_default = 40, mld_spar = 0.3,
                        use_entrainment = TRUE, label = NULL) {
  if (is.null(label)) label <- time_step

  dat <- dat |>
    filter(lon >= lon_min & lon <= lon_max, lat >= lat_min & lat <= lat_max,
           date > as.Date(date_start), date < as.Date(date_end)) |>
    mutate(zeu = zeu_default)

  if (nrow(dat) == 0) stop("No data after filtering for: ", label)

  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |> arrange(date) |> filter(!is.na(MLD), !is.na(zeu))

  dat_smooth <- dat_prof
  fit_mld    <- loess(MLD ~ as.numeric(date), data = dat_smooth,
                      span = mld_spar, family = "symmetric")
  fit_zeu    <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)

  dat_prof <- dat_prof |>
    mutate(
      mld_smooth = pmax(predict(fit_mld, newdata = data.frame(date = as.numeric(dat_prof$date))), 0),
      zeu_smooth = predict(fit_zeu, as.numeric(dat_prof$date))$y
    )

  dat <- left_join(dat, select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth),
                   by = c("float_wmo", "prof_number"))

  time_grid <- tibble(date_grid = seq(min(dat_prof$date), max(dat_prof$date), by = time_step))

  prof_dat_smoothed <- dat |>
    group_by(float_wmo, prof_number, date) |>
    summarise(mld_smooth = unique(mld_smooth), zeu_smooth = unique(zeu_smooth), .groups = "drop") |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid) |>
    summarise(mld = mean(mld_smooth, na.rm = TRUE), zeu = mean(zeu_smooth, na.rm = TRUE), .groups = "drop") |>
    mutate(
      prev_MLD               = lag(mld), prev_zeu = lag(zeu),
      NCP_integration_depth  = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
      next_integration_depth = lead(NCP_integration_depth)
    ) |>
    full_join(time_grid, by = "date_grid") |> arrange(date_grid) |>
    mutate(
      NCP_integration_depth  = zoo::na.approx(NCP_integration_depth,  x = as.numeric(date_grid), na.rm = FALSE, rule = 2),
      next_integration_depth = zoo::na.approx(next_integration_depth, x = as.numeric(date_grid), na.rm = FALSE, rule = 2),
      dt = as.numeric(date_grid - lag(date_grid)),
      NCP_integration_depth_smooth  = rollapply(NCP_integration_depth,  width = 3, FUN = mean, align = "center", fill = NA),
      next_integration_depth_smooth = rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA)
    ) |>
    filter(date_grid %in% time_grid$date_grid) |> na.omit()

  dat_smoothed <- dat |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    full_join(time_grid, by = "date_grid") |> arrange(date_grid) |>
    group_by(depth) |>
    mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid), na.rm = FALSE, rule = 2)) |>
    ungroup() |>
    filter(date_grid %in% time_grid$date_grid) |>
    left_join(prof_dat_smoothed, by = "date_grid") |> na.omit()

  dat_final <- dat_smoothed |>
    group_by(date_grid, mld, zeu) |>
    filter(!is.na(next_integration_depth)) |>
    summarise(
      int_N_mmol_m2         = integrate_nitrate(nitrate, depth, 0, unique(NCP_integration_depth)),
      next_int_N_mmol_m2    = integrate_nitrate(nitrate, depth, 0, unique(next_integration_depth)),
      mld_concentration     = integrate_nitrate(nitrate, depth, unique(NCP_integration_depth) - 20, unique(NCP_integration_depth) - 10) / 10,
      sub_mld_concentration = integrate_nitrate(nitrate, depth, unique(NCP_integration_depth) + 20, unique(NCP_integration_depth) + 30) / 10,
      .groups = "drop"
    )

  ncp_results <- dat_final |>
    mutate(
      int_N_smooth      = rollapply(int_N_mmol_m2,      width = 3, FUN = mean, align = "center", fill = NA),
      next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA)
    ) |>
    arrange(date_grid) |>
    mutate(
      dt                  = as.numeric(date_grid - lag(date_grid)),
      diff_mld            = mld - lag(mld),
      we                  = pmax(0, diff_mld),
      delta               = sub_mld_concentration - mld_concentration,
      diff                = if (use_entrainment) delta * we else 0,
      nitrate_consumption = lag(next_int_N_smooth) - int_N_smooth,
      c_consumption       = nitrate_consumption * 6.625 / dt,
      NCP                 = c_consumption + diff,
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP),
      datenum             = as.numeric(date_grid)
    ) |>
    filter(!is.na(NCP))

  safe_span <- max(0.4, 6 / nrow(ncp_results))
  fit_loess <- loess(NCP ~ datenum, data = ncp_results, span = safe_span, family = "symmetric")
  ncp_results$loess_ncp       <- predict(fit_loess)
  ncp_results$entrainment     <- if (use_entrainment) "with entrainment" else "without entrainment"
  ncp_results$time_step_label <- label
  ncp_results
}

# ── Run both ──────────────────────────────────────────────────────────────────
args <- list(
  dat        = dat,
  time_step  = "10 days",
  lon_min    = lon_min, lon_max = lon_max,
  lat_min    = lat_min, lat_max = lat_max,
  date_start = date_start, date_end = date_end,
  mld_spar   = mld_spar, zeu_default = zeu_default
)

res_with    <- do.call(compute_ncp, c(args, use_entrainment = TRUE,  label = "with entrainment"))
res_without <- do.call(compute_ncp, c(args, use_entrainment = FALSE, label = "without entrainment"))

combined <- bind_rows(res_with, res_without)

# ── Save CSV ──────────────────────────────────────────────────────────────────
out_csv <- file.path(out_dir, "ncp_entrainment_comparison.csv")
write_csv(combined, out_csv)
message("Results saved -> ", out_csv)

# ── Plot ──────────────────────────────────────────────────────────────────────
# Panel 1: raw + loess NCP for both
# Panel 2: difference (with - without) = the entrainment term alone
diff_df <- inner_join(
  select(res_with,    date_grid, NCP_with    = NCP, loess_with    = loess_ncp),
  select(res_without, date_grid, NCP_without = NCP, loess_without = loess_ncp),
  by = "date_grid"
) |>
  mutate(
    entrainment_term       = NCP_with - NCP_without,
    entrainment_term_loess = loess_with - loess_without
  )

p1 <- ggplot(combined) +
  geom_line(aes(x = date_grid, y = NCP,       color = entrainment), alpha = 0.3) +
  geom_line(aes(x = date_grid, y = loess_ncp, color = entrainment), linewidth = 1.2) +
  scale_color_manual(values = c("with entrainment" = "#E41A1C", "without entrainment" = "#377EB8")) +
  labs(y = "NCP  (mmol C m-2 d-1)", x = NULL, color = NULL,
       title = "NCP with and without entrainment correction  (10-day step)") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")

p2 <- ggplot(diff_df) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(x = date_grid, y = entrainment_term),       alpha = 0.3, color = "#4DAF4A") +
  geom_line(aes(x = date_grid, y = entrainment_term_loess), linewidth = 1.2, color = "#4DAF4A") +
  labs(y = "Entrainment term  (mmol C m-2 d-1)", x = "Date",
       title = "Entrainment correction alone  (with - without)") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")

library(patchwork)
p_combined <- p1 / p2 + plot_layout(heights = c(2, 1))

out_png <- file.path(out_dir, "ncp_entrainment_comparison.png")
ggsave(out_png, p_combined, width = 10, height = 7, dpi = 150)
message("Plot saved -> ", out_png)
