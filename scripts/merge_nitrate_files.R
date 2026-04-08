library(tidyverse)
library(gsw)
library(here)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  in_dir   <- snakemake@input[["in_dir"]]
  out_csv  <- snakemake@output[["merged_csv"]]
  out_map  <- snakemake@output[["map_png"]]
  lon_min  <- snakemake@config[["lon_min"]]
  lon_max  <- snakemake@config[["lon_max"]]
  lat_min  <- snakemake@config[["lat_min"]]
  lat_max  <- snakemake@config[["lat_max"]]
} else {
  args     <- commandArgs(trailingOnly = TRUE)
  in_dir   <- ifelse(length(args) > 0, args[1], "data/NorthAtlantic_2024/intermediate/nitrate_profiles")
  out_csv  <- ifelse(length(args) > 1, args[2], "data/NorthAtlantic_2024/intermediate/merged/merged_ncp.csv")
  out_map  <- ifelse(length(args) > 2, args[3], "data/NorthAtlantic_2024/intermediate/merged/float_map.png")
  lon_min  <- -45; lon_max <- -10
  lat_min  <-  54; lat_max <-  63
}

# ── Setup ─────────────────────────────────────────────────────────────────────
out_dir <- dirname(out_csv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Merge all per-float CSV files ─────────────────────────────────────────────
files <- list.files(in_dir, pattern = "_nitrate\\.csv$", full.names = TRUE)
# exclude the manifest
files <- files[!grepl("manifest", files)]
message("Merging ", length(files), " float files from ", in_dir)

dat <- map_dfr(files, function(file) {
  dat_temp <- read_csv(file, show_col_types = FALSE)
  dat_temp$float_wmo <- strsplit(basename(file), "_")[[1]][2]
  dat_temp
})

message("Total rows before filtering: ", nrow(dat))

# ── Derived variables ─────────────────────────────────────────────────────────
dat <- dat |>
  filter(!is.na(lon)) |>
  mutate(
    sa     = gsw::gsw_SA_from_SP(sal,  depth, lon, lat),
    ct     = gsw_CT_from_t(sa,    temp,  depth),
    sigma0 = gsw_sigma0(sa, ct)
  )

# ── MLD per profile ───────────────────────────────────────────────────────────
prof_dat <- dat |>
  group_by(prof_number, float_wmo, date) |>
  arrange(depth) |>
  summarise(
    MLD = castr::mld(sigma0, depth,
                  ref.depths  = 0:5,
                  criteria    = 0.03,
                  default.depth = 300),
    .groups = "drop"
  )

dat_all <- left_join(dat, prof_dat, by = c("prof_number", "float_wmo", "date"))

# ── Float track map ───────────────────────────────────────────────────────────
map_data <- dat_all |>
  select(float_wmo, lon, lat, date) |>
  distinct()

p_map <- ggplot(map_data) +
  geom_point(aes(x = lon, y = lat, color = date)) +
  geom_rect(
    aes(xmin = lon_min, xmax = lon_max,
        ymin = lat_min, ymax = lat_max),
    fill = NA, color = "red", inherit.aes = FALSE
  ) +
  coord_quickmap() +
  labs(title = "Float tracks", x = "Longitude", y = "Latitude") +
  theme_minimal()

ggsave(out_map, p_map, width = 8, height = 6, dpi = 150)
message("Map saved → ", out_map)

# ── MLD diagnostic plot ───────────────────────────────────────────────────────
p_mld <- ggplot(prof_dat) +
  geom_point(aes(x = date, y = -MLD)) +
  labs(title = "Mixed Layer Depth", x = "Date", y = "MLD (m)") +
  theme_minimal()

mld_plot_path <- file.path(out_dir, "mld_timeseries.png")
ggsave(mld_plot_path, p_mld, width = 8, height = 4, dpi = 150)
message("MLD plot saved → ", mld_plot_path)

# ── Select final columns and save ─────────────────────────────────────────────
merged_out <- dat_all |>
  select(float_wmo, prof_number, date, lon, lat,
         canyon_nitrate, sa, ct, sigma0, MLD,
         sal, oxygen, depth)

write_csv(merged_out, out_csv)
message("Merged CSV saved → ", out_csv, " (", nrow(merged_out), " rows)")