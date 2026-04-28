library(tidyverse)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  float_csvs <- snakemake@input[["float_csvs"]]
  wmo_list   <- as.character(snakemake@params[["wmo_list"]])
  fig_path   <- snakemake@output[["fig"]]
} else {
  args       <- commandArgs(trailingOnly = TRUE)
  fig_path   <- ifelse(length(args) > 0, args[1],
                       "output/NorthAtlantic_20s/ncp_float/floats_comparison.png")
  float_csvs <- list.files(dirname(dirname(fig_path)), "ncp_float.csv",
                            recursive = TRUE, full.names = TRUE)
  wmo_list   <- basename(dirname(float_csvs))
}

dir.create(dirname(fig_path), recursive = TRUE, showWarnings = FALSE)

all_ncp <- map_dfr(float_csvs, read_csv, show_col_types = FALSE) |>
  mutate(float_wmo = as.character(float_wmo))

if (nrow(all_ncp) == 0)
  stop("No NCP data available across all selected floats — nothing to plot")

# ── Weekly-binned ensemble mean ───────────────────────────────────────────────
mean_ncp <- all_ncp |>
  mutate(date_week = as.Date(cut(date, breaks = "1 week"))) |>
  group_by(date_week) |>
  summarise(NCP_mean = mean(NCP, na.rm = TRUE), .groups = "drop")

# ── Comparison plot ───────────────────────────────────────────────────────────
p <- ggplot(all_ncp, aes(x = date, y = NCP,
                          color = float_wmo, group = float_wmo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(
    data        = mean_ncp,
    aes(x = date_week, y = NCP_mean),
    inherit.aes = FALSE,
    color       = "black",
    linewidth   = 1.5
  ) +
  labs(
    title    = "NCP — per-float comparison  (black = ensemble mean)",
    y        = "NCP  (mmol C m⁻² d⁻¹)",
    x        = "Date",
    color    = "WMO"
  ) +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months")

ggsave(fig_path, p, width = 12, height = 5, dpi = 150)
message("Comparison plot saved -> ", fig_path)
