library(ncdf4)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(readr)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  sprof_dir     <- snakemake@params[["sprof_dir"]]
  shared_dir    <- snakemake@params[["shared_dir"]]
  manifest_path <- snakemake@output[["manifest"]]
} else {
  args          <- commandArgs(trailingOnly = TRUE)
  sprof_dir     <- ifelse(length(args) > 0, args[1], "data/raw")
  shared_dir    <- ifelse(length(args) > 1, args[2], "data/shared")
  manifest_path <- ifelse(length(args) > 2, args[3],
                          "data/NorthAtlantic_seas_comparison/intermediate/doxy_profiles/processing_manifest.csv")
}

# ── Setup ─────────────────────────────────────────────────────────────────────
shared_doxy_dir <- file.path(shared_dir, "doxy_profiles")
dir.create(shared_doxy_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)

files <- list.files(sprof_dir, pattern = "_Sprof.nc", full.names = TRUE)
message("Found ", length(files), " Sprof files in ", sprof_dir)

pb <- txtProgressBar(min = 0, max = length(files), style = 3)
i  <- 0

for (file in files) {

  i <- i + 1
  setTxtProgressBar(pb, i)

  nc  <- nc_open(file)
  wmo <- gsub(" ", "", ncvar_get(nc, "PLATFORM_NUMBER")[1])

  filename <- file.path(shared_doxy_dir, paste0("argo_", wmo, "_interp.csv"))

  if (file.exists(filename)) {
    message("Skipping ", wmo, " (already in shared)")
    nc_close(nc)
    next
  }

  if (!"DOXY_ADJUSTED" %in% names(nc$var)) {
    message("DOXY_ADJUSTED not found: ", file)
    nc_close(nc)
    next
  }

  # extract metadata
  lon  <- ncvar_get(nc, "LONGITUDE")
  lat  <- ncvar_get(nc, "LATITUDE")
  juld <- ncvar_get(nc, "JULD")
  date <- as.Date(juld, origin = "1950-01-01")

  # extract profiles
  depth  <- ncvar_get(nc, "PRES")
  temp   <- ncvar_get(nc, "TEMP")
  sal    <- ncvar_get(nc, "PSAL")
  oxygen <- ncvar_get(nc, "DOXY_ADJUSTED")
  qc     <- ncvar_get(nc, "DOXY_ADJUSTED_QC")
  nc_close(nc)

  # dimensions
  n_levels   <- dim(depth)[1]
  n_profiles <- dim(depth)[2]

  # repeat metadata across depth levels
  lon_rep  <- rep(lon,  each = n_levels)
  lat_rep  <- rep(lat,  each = n_levels)
  date_rep <- rep(date, each = n_levels)

  # QC filter — keep flags 1 and 2 only (ASCII 49, 50)
  doxy_qc <- as.numeric(charToRaw(paste(qc, collapse = "")))
  oxygen  <- as.vector(oxygen)
  oxygen[!doxy_qc %in% c(49, 50)] <- NA

  df <- tibble(
    lon    = lon_rep,
    lat    = lat_rep,
    date   = date_rep,
    depth  = as.vector(depth),
    temp   = as.vector(temp),
    sal    = as.vector(sal),
    oxygen = oxygen
  )

  df_nona <- filter(df, !is.na(temp)) |>
    mutate(depth = round(depth))

  df_interp <- df_nona %>%
    group_by(lon, lat, date) %>%
    arrange(depth, .by_group = TRUE) %>%
    reframe(depth = seq(0, 2000, by = 1)) %>%
    left_join(df_nona, by = c("lon", "lat", "date", "depth")) %>%
    arrange(lon, lat, date, depth) %>%
    group_by(lon, lat, date) %>%
    mutate(across(
      c(oxygen, temp, sal),
      ~ {
          vals <- .x
          if (sum(!is.na(vals)) > 10) {
            tryCatch(
              na.approx(vals, x = depth, na.rm = FALSE, rule = 2),
              error = function(e) rep(NA_real_, length(vals))
            )
          } else {
            rep(NA_real_, length(vals))
          }
        }
    )) %>%
    ungroup() |>
    group_by(lon, lat, date, depth) |>
    summarise(
    temp = mean(temp, na.rm = TRUE),
    sal = mean(sal, na.rm = TRUE),
    oxygen = mean(oxygen, na.rm = TRUE),
    .groups = "drop") |>
    arrange(date) |>
    mutate(prof_number = dense_rank(date))

  write_csv(df_interp, filename)
  message(filename, " written")
}

close(pb)
message("Done. Processed ", i, " floats.")

# Write manifest pointing to all shared files
shared_files <- list.files(shared_doxy_dir, pattern = "_interp\\.csv$", full.names = TRUE)
manifest <- data.frame(
  wmo  = gsub(".*argo_(.*)_interp\\.csv", "\\1", basename(shared_files)),
  path = shared_files
)
write_csv(manifest, manifest_path)
message("Manifest -> ", manifest_path)
