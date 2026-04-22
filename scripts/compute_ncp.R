library(tidyverse)
library(zoo)
library(sf)

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
    basin_poly,
    date_start, date_end,
    zeu_default    = 40,
    label          = NULL,
    mld_diag_path  = NULL
) {

  if (is.null(label)) label <- time_step

  # fast bbox pre-filter before expensive st_within
  bbox     <- st_bbox(basin_poly)
  dat <- dat |> filter(lon >= bbox["xmin"], lon <= bbox["xmax"],
                       lat >= bbox["ymin"], lat <= bbox["ymax"])

  # exact polygon filter
  dat_sf   <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  in_basin <- lengths(st_within(dat_sf, basin_poly)) > 0
  dat <- dat[in_basin, ] |>
    filter(
      date > as.Date(date_start),
      date < as.Date(date_end)
    ) |>
    mutate(zeu = zeu_default)

  if (nrow(dat) == 0) stop("No data after filtering for: ", label)
  message(label, ": ", nrow(dat), " rows after filtering")

  # per-profile MLD/Zeu
  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |>
    filter(!is.na(MLD), !is.na(zeu))

  # time grid
  time_grid <- tibble(
    date_grid = seq(min(dat_prof$date), max(dat_prof$date), by = time_step)
  )

  # integration depth — bin-average MLD and Zeu per time step
  prof_dat_smoothed <- dat_prof |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid) |>
    summarise(mld = mean(MLD, na.rm = TRUE),
              zeu = mean(zeu, na.rm = TRUE),
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

  # ── MLD diagnostic plot (optional) ──────────────────────────────────────────
  if (!is.null(mld_diag_path)) {
    dat_prof_plot <- dat_prof |>
      mutate(date_grid = as.Date(cut(date, breaks = time_step)))
    p_mld <- ggplot() +
      geom_point(data = dat_prof_plot,
                 aes(x = date, y = MLD, color = factor(float_wmo)),
                 alpha = 0.5, size = 1.2, show.legend = TRUE) +
      geom_line(data = prof_dat_smoothed |> filter(!is.na(mld)),
                aes(x = date_grid, y = mld),
                color = "black", linewidth = 1.2, linetype = "solid") +
      geom_point(data = prof_dat_smoothed |> filter(!is.na(mld)),
                 aes(x = date_grid, y = mld),
                 color = "black", size = 2.5, shape = 21, fill = "white") +
      scale_y_reverse() +
      scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
      labs(x = "Date", y = "MLD (m)",
           title = paste0(basin_name, " — MLD: individual profiles vs bin average (", time_step, ")"),
           color = "Float WMO") +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7))
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
write_csv(all_results, results_csv)
message("NCP results saved -> ", results_csv)

# ── Plots ─────────────────────────────────────────────────────────────────────
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
