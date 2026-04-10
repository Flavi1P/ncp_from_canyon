suppressPackageStartupMessages({library(tidyverse); library(zoo); library(patchwork)})

merged_csv  <- "data/NorthAtlantic_20s/intermediate/merged/merged_ncp.csv"
out_dir     <- "output/NorthAtlantic_20s/ncp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lon_min <- -45; lon_max <- -10; lat_min <- 54; lat_max <- 63
date_start  <- "2020-01-01"; date_end <- "2026-01-01"
mld_spar    <- 0.3; zeu_default <- 40
n_mc        <- 200; canyon_rmse <- 1.2
time_step   <- "10 days"

dat <- read_csv(merged_csv, show_col_types = FALSE)
message("Loaded ", nrow(dat), " rows")

# ── Helpers ───────────────────────────────────────────────────────────────────
integrate_nitrate <- function(nitrate, depth, from, to) {
  idx <- depth >= from & depth <= to
  if (sum(idx) < 2) return(NA_real_)
  mean(approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)), na.rm = TRUE) * (to - from)
}

# Shared setup: filter, smooth MLD, build time grid, compute integration depths.
# Returns a list so both the deterministic and MC paths reuse the same objects.
prepare_dat <- function(dat, time_step, lon_min, lon_max, lat_min, lat_max,
                        date_start, date_end, zeu_default, mld_spar) {
  dat <- dat |>
    filter(lon >= lon_min & lon <= lon_max, lat >= lat_min & lat <= lat_max,
           date > as.Date(date_start), date < as.Date(date_end)) |>
    mutate(zeu = zeu_default)

  if (nrow(dat) == 0) stop("No data after filtering")

  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |> arrange(date) |> filter(!is.na(MLD), !is.na(zeu))

  fit_mld <- loess(MLD ~ as.numeric(date), data = dat_prof,
                   span = mld_spar, family = "symmetric")
  fit_zeu <- smooth.spline(as.numeric(dat_prof$date), dat_prof$zeu, spar = 0.6)

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

  list(dat = dat, time_grid = time_grid, prof_dat_smoothed = prof_dat_smoothed)
}

ncp_from_smoothed <- function(dat_smoothed, use_entrainment) {
  dat_final <- dat_smoothed |>
    group_by(date_grid, mld, zeu) |>
    filter(!is.na(next_integration_depth)) |>
    summarise(
      int_N_mmol_m2         = integrate_nitrate(nitrate, depth, 0, unique(NCP_integration_depth)),
      next_int_N_mmol_m2    = integrate_nitrate(nitrate, depth, 0, unique(next_integration_depth)),
      mld_concentration     = integrate_nitrate(nitrate, depth, unique(NCP_integration_depth) - 20, unique(NCP_integration_depth) - 10) / 10,
      sub_mld_concentration = integrate_nitrate(nitrate, depth, unique(NCP_integration_depth) + 20, unique(NCP_integration_depth) + 30) / 10,
      .groups = "drop"
    ) |>
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
      NCP                 = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP)
    ) |>
    filter(!is.na(NCP)) |>
    select(date_grid, NCP)
}

build_mean_nitrate <- function(dat, time_step, time_grid, prof_dat_smoothed) {
  dat |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    full_join(time_grid, by = "date_grid") |> arrange(date_grid) |>
    group_by(depth) |>
    mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid), na.rm = FALSE, rule = 2)) |>
    ungroup() |>
    filter(date_grid %in% time_grid$date_grid) |>
    left_join(prof_dat_smoothed, by = "date_grid") |>
    na.omit()
}

# ── Prepare shared objects ────────────────────────────────────────────────────
message("Preparing data ...")
prep <- prepare_dat(dat, time_step, lon_min, lon_max, lat_min, lat_max,
                    date_start, date_end, zeu_default, mld_spar)
dat_f             <- prep$dat
time_grid         <- prep$time_grid
prof_dat_smoothed <- prep$prof_dat_smoothed
message(nrow(dat_f), " rows after filtering; ", nrow(time_grid), " time bins")

# ── Deterministic runs ────────────────────────────────────────────────────────
message("Running deterministic NCP ...")
mean_nitrate <- build_mean_nitrate(dat_f, time_step, time_grid, prof_dat_smoothed)

det_with    <- ncp_from_smoothed(mean_nitrate, use_entrainment = TRUE)  |> mutate(run = "with entrainment")
det_without <- ncp_from_smoothed(mean_nitrate, use_entrainment = FALSE) |> mutate(run = "without entrainment")

loess_smooth <- function(df) {
  df <- df |> mutate(datenum = as.numeric(date_grid))
  span <- max(0.4, 6 / nrow(df))
  fit  <- loess(NCP ~ datenum, data = df, span = span, family = "symmetric")
  df$loess_ncp <- predict(fit)
  df
}
det_with    <- loess_smooth(det_with)
det_without <- loess_smooth(det_without)

# ── Monte Carlo runs ──────────────────────────────────────────────────────────
run_mc <- function(use_entrainment, label) {
  message(label, ": running ", n_mc, " MC iterations ...")
  dat_binned   <- dat_f |> mutate(date_grid = as.Date(cut(date, breaks = time_step)))
  profile_pool <- dat_binned |> select(date_grid, float_wmo, prof_number) |> distinct()

  mc_ncp <- map(seq_len(n_mc), function(i) {
    tryCatch({
      resampled <- profile_pool |>
        group_by(date_grid) |>
        slice_sample(prop = 1, replace = TRUE) |>
        ungroup() |>
        mutate(draw_id = row_number())

      dat_resampled <- resampled |>
        left_join(dat_binned |> select(date_grid, float_wmo, prof_number, depth, canyon_nitrate),
                  by = c("date_grid", "float_wmo", "prof_number"),
                  relationship = "many-to-many") |>
        group_by(draw_id) |>
        mutate(nitrate = canyon_nitrate + rnorm(1, 0, canyon_rmse)) |>
        ungroup()

      nitrate_binned <- dat_resampled |>
        group_by(date_grid, depth) |>
        summarise(nitrate = mean(nitrate, na.rm = TRUE), .groups = "drop")

      dat_smoothed <- nitrate_binned |>
        full_join(time_grid, by = "date_grid") |> arrange(date_grid) |>
        group_by(depth) |>
        mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid), na.rm = FALSE, rule = 2)) |>
        ungroup() |>
        filter(date_grid %in% time_grid$date_grid) |>
        left_join(prof_dat_smoothed, by = "date_grid") |>
        na.omit()

      ncp_from_smoothed(dat_smoothed, use_entrainment = use_entrainment) |> mutate(iter = i)
    }, error = function(e) NULL)
  }) |> compact() |> bind_rows()

  n_ok <- n_distinct(mc_ncp$iter)
  message(label, ": ", n_ok, "/", n_mc, " iterations succeeded")

  mc_ncp |>
    group_by(date_grid) |>
    summarise(NCP_q05 = quantile(NCP, 0.05, na.rm = TRUE),
              NCP_q95 = quantile(NCP, 0.95, na.rm = TRUE),
              .groups = "drop") |>
    mutate(run = label)
}

mc_with    <- run_mc(TRUE,  "with entrainment")
mc_without <- run_mc(FALSE, "without entrainment")

# ── Combine ───────────────────────────────────────────────────────────────────
det_all <- bind_rows(det_with, det_without)
mc_all  <- bind_rows(mc_with,  mc_without)

plot_df <- det_all |>
  left_join(mc_all, by = c("date_grid", "run"))

# ── Save CSV ──────────────────────────────────────────────────────────────────
out_csv <- file.path(out_dir, "ncp_entrainment_mc_comparison.csv")
write_csv(plot_df, out_csv)
message("Results saved -> ", out_csv)

# ── Plot ──────────────────────────────────────────────────────────────────────
colours <- c("with entrainment" = "#E41A1C", "without entrainment" = "#377EB8")

p1 <- ggplot(plot_df, aes(x = date_grid, color = run, fill = run)) +
  geom_ribbon(aes(ymin = NCP_q05, ymax = NCP_q95), alpha = 0.15, color = NA) +
  geom_line(aes(y = NCP),       alpha = 0.25) +
  geom_line(aes(y = loess_ncp), linewidth = 1.1) +
  scale_color_manual(values = colours) +
  scale_fill_manual( values = colours) +
  labs(y = "NCP  (mmol C m-2 d-1)", x = NULL, color = NULL, fill = NULL,
       title    = "NCP with and without entrainment correction  (10-day step, 2020-2026)",
       subtitle = paste0("Shading = 5-95th percentile across ", n_mc, " MC iterations")) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

# Difference panel: loess_with - loess_without = entrainment term
diff_df <- inner_join(
  select(det_with,    date_grid, loess_with    = loess_ncp, NCP_with    = NCP),
  select(det_without, date_grid, loess_without = loess_ncp, NCP_without = NCP),
  by = "date_grid"
) |>
  mutate(entrainment_term       = NCP_with - NCP_without,
         entrainment_term_loess = loess_with - loess_without)

p2 <- ggplot(diff_df, aes(x = date_grid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = entrainment_term),       alpha = 0.25, color = "#4DAF4A") +
  geom_line(aes(y = entrainment_term_loess), linewidth = 1.1, color = "#4DAF4A") +
  labs(y = "Entrainment term  (mmol C m-2 d-1)", x = "Date",
       title = "Entrainment correction alone  (with - without)") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

p_out <- p1 / p2 + plot_layout(heights = c(2, 1))

out_png <- file.path(out_dir, "ncp_entrainment_mc_comparison.png")
ggsave(out_png, p_out, width = 14, height = 8, dpi = 150)
message("Plot saved -> ", out_png)
