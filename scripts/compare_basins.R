library(tidyverse)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  csv_files   <- unlist(snakemake@input)
  basin_names <- snakemake@params[["basin_names"]]
  out_fig     <- snakemake@output[["fig"]]
} else {
  run         <- "NorthAtlantic_20s"
  basin_names <- c("IcelandBasin", "IrmingerSea")
  csv_files   <- file.path("output", run, "ncp", basin_names, "ncp_uncertainty.csv")
  out_fig     <- file.path("output", run, "ncp", "basins_comparison.png")
}

dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)

# ── Load and tag each basin ───────────────────────────────────────────────────
dat <- map2_dfr(csv_files, basin_names, function(f, bn) {
  read_csv(f, show_col_types = FALSE) |> mutate(basin = bn)
})

message("Loaded ", nrow(dat), " rows across ", n_distinct(dat$basin), " basins")

# ── One panel per time step, basins overlaid ──────────────────────────────────
p <- ggplot(dat, aes(x = date_grid, color = basin, fill = basin)) +
  geom_ribbon(aes(ymin = NCP_q05, ymax = NCP_q95), alpha = 0.15, color = NA) +
  geom_line(aes(y = NCP_mean), linewidth = 1) +
  facet_wrap(~time_step_label, ncol = 1) +
  scale_color_brewer(palette = "Dark2", name = "Basin") +
  scale_fill_brewer( palette = "Dark2", name = "Basin") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    y        = "NCP  (mmol C m-2 d-1)",
    x        = "Date",
    title    = "NCP comparison across basins",
    subtitle = "Shading = 5-95th percentile of Monte Carlo uncertainty"
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

ggsave(out_fig, p,
       width  = 12,
       height = 4 * n_distinct(dat$time_step_label),
       dpi    = 150)
message("Basin comparison figure saved -> ", out_fig)
