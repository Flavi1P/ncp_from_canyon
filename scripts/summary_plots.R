library(tidyverse)
library(lubridate)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  basin_unc_csvs <- as.character(unlist(snakemake@input[["basin_unc_csvs"]]))
  basin_res_csvs <- as.character(unlist(snakemake@input[["basin_res_csvs"]]))
  float_csvs     <- as.character(unlist(snakemake@input[["float_csvs"]]))
  basin_names    <- as.character(unlist(snakemake@params[["basin_names"]]))
  float_wmos     <- as.character(unlist(snakemake@params[["float_wmos"]]))
  out_dir        <- snakemake@output[["out_dir"]]
  marker_path    <- snakemake@output[["marker"]]
} else {
  run            <- "Single_float_test"
  basin_names    <- c("AtlanticSAZ_PFZ", "AtlanticSubtropical", "IndianKerguelen")
  basin_unc_csvs <- file.path("output", run, "ncp", basin_names, "ncp_uncertainty.csv")
  basin_res_csvs <- file.path("output", run, "ncp", basin_names, "ncp_results.csv")
  float_wmos     <- c("5904845", "5905099", "5906206")
  float_csvs     <- file.path("output", run, "ncp_float", float_wmos, "ncp_float.csv")
  out_dir        <- file.path("output", run, "summary")
  marker_path    <- file.path(out_dir, "_done.txt")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ───────────────────────────────────────────────────────────────────
read_tagged <- function(path, tag, tag_col = "region") {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d[[tag_col]] <- tag
  d
}

ts_days <- function(label) suppressWarnings(as.numeric(str_extract(label, "[0-9]+")))

# Pick the longest available time step as the default reference (e.g. "30 days").
pick_ref_ts <- function(labels) {
  if (length(labels) == 0) return(NA_character_)
  u <- unique(labels)
  u[which.max(ts_days(u))]
}

# ── Read basin uncertainty + results CSVs ─────────────────────────────────────
basin_unc <- map2_dfr(basin_unc_csvs, basin_names, read_tagged)
basin_res <- map2_dfr(basin_res_csvs, basin_names, read_tagged)

if (!is.null(basin_unc) && nrow(basin_unc) > 0)
  basin_unc <- basin_unc |> filter(!is.na(NCP_mean))
if (!is.null(basin_res) && nrow(basin_res) > 0)
  basin_res <- basin_res |> filter(!is.na(NCP))

ref_ts <- if (!is.null(basin_unc) && nrow(basin_unc) > 0)
  pick_ref_ts(basin_unc$time_step_label) else NA_character_

message("Loaded basin_unc: ", ifelse(is.null(basin_unc), 0, nrow(basin_unc)),
        " rows, basin_res: ", ifelse(is.null(basin_res), 0, nrow(basin_res)),
        " rows. Reference timestep: ", ref_ts)

# ── Derived datasets ──────────────────────────────────────────────────────────
ancp_basins <- if (!is.null(basin_unc) && nrow(basin_unc) > 0) {
  basin_unc |>
    mutate(year     = year(date_grid),
           timestep = ts_days(time_step_label),
           tot_ncp  = NCP_mean * timestep,
           tot_q05  = NCP_q05  * timestep,
           tot_q95  = NCP_q95  * timestep,
           tot_sd   = NCP_sd   * timestep) |>
    group_by(region, year, time_step_label) |>
    summarise(ANCP     = sum(tot_ncp, na.rm = TRUE) / 1000,
              ANCP_q05 = sum(tot_q05, na.rm = TRUE) / 1000,
              ANCP_q95 = sum(tot_q95, na.rm = TRUE) / 1000,
              ANCP_sd  = sqrt(sum(tot_sd**2, na.rm = TRUE)) / 1000,
              conf_int = sqrt(sum(tot_sd**2, na.rm = TRUE)) / 1000 * 1.96,
              .groups  = "drop")
} else NULL

climatology <- if (!is.null(basin_unc) && nrow(basin_unc) > 0) {
  basin_unc |>
    mutate(month = month(date_grid)) |>
    group_by(region, month, time_step_label) |>
    summarise(NCP = mean(NCP_mean, na.rm = TRUE),
              sd  = mean(NCP_sd,   na.rm = TRUE),
              q05 = mean(NCP_q05,  na.rm = TRUE),
              q95 = mean(NCP_q95,  na.rm = TRUE),
              .groups = "drop")
} else NULL

ANCP_clima <- if (!is.null(climatology) && nrow(climatology) > 0 && !is.na(ref_ts)) {
  climatology |>
    filter(time_step_label == ref_ts) |>
    group_by(region) |>
    summarise(ancp = sum(NCP * 30.5, na.rm = TRUE) / 1000,
              aq05 = sum(q05 * 30.5, na.rm = TRUE) / 1000,
              aq95 = sum(q95 * 30.5, na.rm = TRUE) / 1000,
              sd   = sqrt(sum((sd * 30.5)**2, na.rm = TRUE)) / 1000,
              .groups = "drop")
} else NULL

ANCP_clima_timestep <- if (!is.null(climatology) && nrow(climatology) > 0) {
  climatology |>
    group_by(region, time_step_label) |>
    summarise(ancp = sum(NCP * 30.5, na.rm = TRUE) / 1000,
              aq05 = sum(q05 * 30.5, na.rm = TRUE) / 1000,
              aq95 = sum(q95 * 30.5, na.rm = TRUE) / 1000,
              sd   = sqrt(sum((sd * 30.5)**2, na.rm = TRUE)) / 1000,
              .groups = "drop")
} else NULL

# ── sANCP variants (NCP clipped at 0 — only positive contributions) ───────────
sancp_basins <- if (!is.null(basin_unc) && nrow(basin_unc) > 0) {
  basin_unc |>
    mutate(year     = year(date_grid),
           timestep = ts_days(time_step_label),
           tot_ncp  = pmax(NCP_mean, 0) * timestep,
           tot_q05  = pmax(NCP_q05,  0) * timestep,
           tot_q95  = pmax(NCP_q95,  0) * timestep,
           tot_sd   = NCP_sd * timestep) |>
    group_by(region, year, time_step_label) |>
    summarise(sANCP     = sum(tot_ncp, na.rm = TRUE) / 1000,
              sANCP_q05 = sum(tot_q05, na.rm = TRUE) / 1000,
              sANCP_q95 = sum(tot_q95, na.rm = TRUE) / 1000,
              sANCP_sd  = sqrt(sum(tot_sd**2, na.rm = TRUE)) / 1000,
              conf_int  = sqrt(sum(tot_sd**2, na.rm = TRUE)) / 1000 * 1.96,
              .groups   = "drop")
} else NULL

sANCP_clima <- if (!is.null(climatology) && nrow(climatology) > 0 && !is.na(ref_ts)) {
  climatology |>
    filter(time_step_label == ref_ts) |>
    mutate(NCP_pos = pmax(NCP, 0),
           q05_pos = pmax(q05, 0),
           q95_pos = pmax(q95, 0)) |>
    group_by(region) |>
    summarise(sancp = sum(NCP_pos * 30.5, na.rm = TRUE) / 1000,
              aq05  = sum(q05_pos * 30.5, na.rm = TRUE) / 1000,
              aq95  = sum(q95_pos * 30.5, na.rm = TRUE) / 1000,
              sd    = sqrt(sum((sd * 30.5)**2, na.rm = TRUE)) / 1000,
              .groups = "drop")
} else NULL

sANCP_clima_timestep <- if (!is.null(climatology) && nrow(climatology) > 0) {
  climatology |>
    mutate(NCP_pos = pmax(NCP, 0),
           q05_pos = pmax(q05, 0),
           q95_pos = pmax(q95, 0)) |>
    group_by(region, time_step_label) |>
    summarise(sancp = sum(NCP_pos * 30.5, na.rm = TRUE) / 1000,
              aq05  = sum(q05_pos * 30.5, na.rm = TRUE) / 1000,
              aq95  = sum(q95_pos * 30.5, na.rm = TRUE) / 1000,
              sd    = sqrt(sum((sd * 30.5)**2, na.rm = TRUE)) / 1000,
              .groups = "drop")
} else NULL

# ── ANCP partition: positive vs negative components per year ──────────────────
# `positive_sum`  = sum(pmax(NCP_mean, 0) * dt)  → equals sANCP
# `negative_sum`  = sum(pmin(NCP_mean, 0) * dt)  → ≤ 0
# net = positive_sum + negative_sum             → equals net ANCP
ancp_partition <- if (!is.null(basin_unc) && nrow(basin_unc) > 0) {
  basin_unc |>
    mutate(year     = year(date_grid),
           timestep = ts_days(time_step_label),
           pos_part = pmax(NCP_mean, 0) * timestep / 1000,
           neg_part = pmin(NCP_mean, 0) * timestep / 1000) |>
    group_by(region, year, time_step_label) |>
    summarise(positive_sum = sum(pos_part, na.rm = TRUE),
              negative_sum = sum(neg_part, na.rm = TRUE),
              .groups = "drop") |>
    mutate(net = positive_sum + negative_sum)
} else NULL

rate_clima <- if (!is.null(basin_res) && nrow(basin_res) > 0 &&
                  all(c("NCP","mld") %in% colnames(basin_res))) {
  basin_res |>
    filter(!is.na(NCP), !is.na(mld), mld > 0) |>
    mutate(rate  = NCP / mld,
           month = month(date_grid)) |>
    group_by(region, month) |>
    summarise(rate_mean = mean(rate, na.rm = TRUE),
              rate_sd   = sd(rate,   na.rm = TRUE),
              .groups   = "drop")
} else NULL

# ── Float ANCP ────────────────────────────────────────────────────────────────
float_data <- if (length(float_csvs) > 0) {
  map2_dfr(float_csvs, float_wmos, read_tagged, tag_col = "wmo")
} else NULL

float_ancp <- if (!is.null(float_data) && nrow(float_data) > 0) {
  float_data |>
    mutate(NCP_pos = pmax(NCP, 0),
           month   = month(date),
           year    = year(date)) |>
    group_by(wmo, month, year) |>
    summarise(NCP    = mean(NCP_pos, na.rm = TRUE),
              NCP_sd = sd(NCP_pos,   na.rm = TRUE),
              .groups = "drop") |>
    mutate(tot_ncp = NCP    * 30.5,
           tot_sd  = NCP_sd * 30.5) |>
    group_by(wmo, year) |>
    summarise(ancp = sum(tot_ncp, na.rm = TRUE) / 1000,
              asd  = sqrt(sum(tot_sd**2, na.rm = TRUE)) / 1000,
              .groups = "drop")
} else NULL

# ── Save plot helper ──────────────────────────────────────────────────────────
save_plot <- function(p, name, w = 12, h = 6) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 200)
  message("Saved -> ", path)
}

n_basins <- length(basin_names)
basin_h  <- max(4, 2.8 * n_basins)

# ── Plots: basin timeseries (all basins overlaid, longest timestep) ──────────
if (!is.null(basin_unc) && nrow(basin_unc) > 0 && !is.na(ref_ts)) {
  p <- ggplot(filter(basin_unc, time_step_label == ref_ts)) +
    geom_line(aes(x = date_grid, y = NCP_mean, color = region)) +
    geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = region),
                alpha = 0.2) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Basin") +
    scale_x_date(date_labels = "%b %Y", breaks = "6 months") +
    labs(y = "NCP (mmol C m-2 d-1)", x = "Date",
         title = paste0(ref_ts, " time window — basin comparison"),
         subtitle = "Shading = MC q05-q95") +
    theme_bw(base_size = 14)
  save_plot(p, "timeseries_all_basins_ref_ts.png", w = 14, h = 6)
}

# ── Plots: basin timeseries faceted, all timesteps ───────────────────────────
if (!is.null(basin_unc) && nrow(basin_unc) > 0) {
  p <- ggplot(basin_unc) +
    geom_line(aes(x = date_grid, y = NCP_mean, color = time_step_label)) +
    geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95,
                    fill = time_step_label), alpha = 0.2) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_color_brewer(palette = "Set2", name = "Time bin") +
    scale_x_date(date_labels = "%b %Y", breaks = "6 months") +
    labs(y = "NCP (mmol C m-2 d-1)", x = "Date") +
    facet_wrap(~ region, ncol = 1, scales = "free_y") +
    theme_bw(base_size = 14)
  save_plot(p, "timeseries_per_basin_per_timestep.png", w = 14, h = basin_h)
}

# ── Plots: yearly ANCP barplot ───────────────────────────────────────────────
if (!is.null(ancp_basins) && nrow(ancp_basins) > 0) {
  p <- ggplot(ancp_basins) +
    geom_col(aes(x = year, y = ANCP, fill = time_step_label),
             position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(x = year, ymin = ANCP - conf_int, ymax = ANCP + conf_int,
                      group = time_step_label),
                  position = position_dodge(width = 0.9), width = 0.3) +
    scale_fill_brewer(palette = "Set2", name = "Time bin") +
    labs(y = "ANCP (mol C m-2 yr-1)", x = "Year") +
    facet_wrap(~ region, ncol = 1) +
    theme_bw(base_size = 14)
  save_plot(p, "ancp_barplot_per_basin_per_timestep.png", w = 14, h = basin_h)
}

# ── Plots: yearly sANCP barplot (positive-only NCP) ──────────────────────────
if (!is.null(sancp_basins) && nrow(sancp_basins) > 0) {
  p <- ggplot(sancp_basins) +
    geom_col(aes(x = year, y = sANCP, fill = time_step_label),
             position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(x = year, ymin = sANCP - conf_int, ymax = sANCP + conf_int,
                      group = time_step_label),
                  position = position_dodge(width = 0.9), width = 0.3) +
    scale_fill_brewer(palette = "Set2", name = "Time bin") +
    labs(y = "sANCP (mol C m-2 yr-1)", x = "Year",
         title = "Yearly sANCP — positive NCP only (NCP < 0 clipped to 0)") +
    facet_wrap(~ region, ncol = 1) +
    theme_bw(base_size = 14)
  save_plot(p, "sancp_barplot_per_basin_per_timestep.png", w = 14, h = basin_h)
}

# ── Plots: ANCP partition (positive vs negative, diverging stack) ────────────
if (!is.null(ancp_partition) && nrow(ancp_partition) > 0 && !is.na(ref_ts)) {
  part_long <- ancp_partition |>
    filter(time_step_label == ref_ts) |>
    pivot_longer(c(positive_sum, negative_sum),
                 names_to = "component", values_to = "ANCP_part")
  net_pts <- ancp_partition |> filter(time_step_label == ref_ts)

  p <- ggplot(part_long) +
    geom_col(aes(x = year, y = ANCP_part, fill = component)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_point(data = net_pts,
               aes(x = year, y = net),
               shape = 23, fill = "white", color = "black", size = 2.8, stroke = 0.8) +
    scale_fill_manual(values = c(positive_sum = "#1f78b4",
                                 negative_sum = "#e31a1c"),
                      labels = c(positive_sum = "Positive (sANCP)",
                                 negative_sum = "Negative (loss)"),
                      name = NULL) +
    labs(y = "ANCP component (mol C m-2 yr-1)", x = "Year",
         title = paste0("Yearly ANCP partition: positive vs negative (",
                        ref_ts, " timestep)"),
         subtitle = "White diamond = net ANCP  |  positive bar height = sANCP") +
    facet_wrap(~ region, ncol = 1, scales = "free_y") +
    theme_bw(base_size = 14)
  save_plot(p, "ancp_partition_per_basin_per_year.png", w = 14, h = basin_h)
}

# ── Plots: monthly climatology faceted ───────────────────────────────────────
if (!is.null(climatology) && nrow(climatology) > 0) {
  p <- ggplot(climatology) +
    geom_line(aes(x = month, y = NCP, color = time_step_label)) +
    geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = time_step_label),
                alpha = 0.2) +
    scale_color_brewer(palette = "Set2", name = "Time step") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_x_continuous(breaks = 1:12) +
    labs(y = "NCP (mmol C m-2 d-1)", x = "Month") +
    facet_wrap(~ region, ncol = 1) +
    theme_bw(base_size = 14)
  save_plot(p, "climatology_per_basin_per_timestep.png", w = 12, h = basin_h)
}

# ── Plots: monthly climatology overlaid (longest timestep) ───────────────────
if (!is.null(climatology) && nrow(climatology) > 0 && !is.na(ref_ts)) {
  p <- ggplot(filter(climatology, time_step_label == ref_ts)) +
    geom_line(aes(x = month, y = NCP, color = region)) +
    geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = region),
                alpha = 0.2) +
    scale_color_brewer(palette = "Set1", name = "Basin") +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_x_continuous(breaks = 1:12) +
    labs(y = "NCP (mmol C m-2 d-1)", x = "Month",
         title = paste0("Monthly climatology — ", ref_ts, " timestep")) +
    theme_bw(base_size = 14)
  save_plot(p, "climatology_all_basins_ref_ts.png", w = 12, h = 6)
}

# ── Plots: ANCP from climatology — single bar per basin ──────────────────────
if (!is.null(ANCP_clima) && nrow(ANCP_clima) > 0) {
  subt <- paste(sprintf("%s: %.2f mol C m-2 yr-1", ANCP_clima$region,
                        ANCP_clima$ancp), collapse = "  |  ")
  p <- ggplot(ANCP_clima) +
    geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8) +
    geom_errorbar(aes(x = region, ymin = ancp - sd, ymax = ancp + sd),
                  width = 0.3) +
    scale_fill_brewer(palette = "Set1", name = "Basin") +
    labs(y = "ANCP (mol C m-2 yr-1)",
         title = paste0("Annual NCP from monthly climatology (", ref_ts, ")"),
         subtitle = subt) +
    theme_bw(base_size = 14)
  save_plot(p, "ancp_clima_per_basin.png", w = 10, h = 6)
}

# ── Plots: ANCP from climatology — dodged by timestep ────────────────────────
if (!is.null(ANCP_clima_timestep) && nrow(ANCP_clima_timestep) > 0) {
  p <- ggplot(ANCP_clima_timestep) +
    geom_col(aes(x = region, y = ancp, fill = time_step_label),
             position = position_dodge(width = 0.9), alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", name = "Time bin") +
    labs(y = "ANCP (mol C m-2 yr-1)",
         title = "Annual NCP from monthly climatology — by timestep") +
    theme_bw(base_size = 14)
  save_plot(p, "ancp_clima_per_basin_per_timestep.png", w = 10, h = 6)
}

# ── Plots: sANCP from climatology — single bar per basin ─────────────────────
if (!is.null(sANCP_clima) && nrow(sANCP_clima) > 0) {
  subt <- paste(sprintf("%s: %.2f mol C m-2 yr-1", sANCP_clima$region,
                        sANCP_clima$sancp), collapse = "  |  ")
  p <- ggplot(sANCP_clima) +
    geom_col(aes(x = region, y = sancp, fill = region), alpha = 0.8) +
    geom_errorbar(aes(x = region, ymin = sancp - sd, ymax = sancp + sd),
                  width = 0.3) +
    scale_fill_brewer(palette = "Set1", name = "Basin") +
    labs(y = "sANCP (mol C m-2 yr-1)",
         title = paste0("Annual sANCP from monthly climatology (", ref_ts,
                        ", positive only)"),
         subtitle = subt) +
    theme_bw(base_size = 14)
  save_plot(p, "sancp_clima_per_basin.png", w = 10, h = 6)
}

# ── Plots: sANCP from climatology — dodged by timestep ───────────────────────
if (!is.null(sANCP_clima_timestep) && nrow(sANCP_clima_timestep) > 0) {
  p <- ggplot(sANCP_clima_timestep) +
    geom_col(aes(x = region, y = sancp, fill = time_step_label),
             position = position_dodge(width = 0.9), alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", name = "Time bin") +
    labs(y = "sANCP (mol C m-2 yr-1)",
         title = "Annual sANCP from monthly climatology — by timestep",
         subtitle = "NCP < 0 clipped to 0 before integration") +
    theme_bw(base_size = 14)
  save_plot(p, "sancp_clima_per_basin_per_timestep.png", w = 10, h = 6)
}

# ── Plots: NCP/MLD rate climatology ──────────────────────────────────────────
if (!is.null(rate_clima) && nrow(rate_clima) > 0) {
  ncols <- if (n_basins <= 2) n_basins else 2
  nrows <- ceiling(n_basins / ncols)
  p <- ggplot(rate_clima) +
    geom_line(aes(x = month, y = rate_mean)) +
    geom_ribbon(aes(x = month,
                    ymin = rate_mean - rate_sd, ymax = rate_mean + rate_sd),
                alpha = 0.4) +
    scale_x_continuous(breaks = 1:12) +
    labs(y = "NCP / MLD  (mmol C m-3 d-1)", x = "Month",
         title = "Volumetric NCP rate climatology") +
    facet_wrap(~ region, ncol = ncols, scales = "free_y") +
    theme_bw(base_size = 14)
  save_plot(p, "rate_climatology_per_basin.png",
            w = 6 * ncols, h = max(4, 3.5 * nrows))
}

# ── Plots: float yearly ANCP ─────────────────────────────────────────────────
if (!is.null(float_ancp) && nrow(float_ancp) > 0) {
  p <- ggplot(float_ancp) +
    geom_col(aes(x = year, y = ancp, fill = wmo),
             position = position_dodge(width = 0.9), alpha = 0.85) +
    scale_fill_brewer(palette = "Set1", name = "WMO") +
    labs(y = "ANCP (mol C m-2 yr-1)", x = "Year",
         title = "Float-mode yearly ANCP",
         subtitle = "From monthly clipped-positive NCP * 30.5 d (per WMO)") +
    facet_wrap(~ wmo, scales = "free_x") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")
  save_plot(p, "float_ancp_per_year.png", w = 12, h = 6)
}

# ── Marker file ──────────────────────────────────────────────────────────────
writeLines(c(paste("Summary plots generated", Sys.time()),
             paste("Basins:", paste(basin_names, collapse = ", ")),
             paste("Floats:", paste(float_wmos, collapse = ", ")),
             paste("Reference timestep:", ref_ts)),
           marker_path)
message("All summary plots written under ", out_dir)
