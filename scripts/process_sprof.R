library(ncdf4)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(readr)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  sprof_dir <- snakemake@input[["sprof_dir"]]
  out_dir   <- snakemake@output[["out_dir"]]
} else {
  # CLI fallback for testing
  args      <- commandArgs(trailingOnly = TRUE)
  sprof_dir <- ifelse(length(args) > 0, args[1], "data/NorthAtlantic_test/raw")
  out_dir   <- ifelse(length(args) > 1, args[2], "data/NorthAtlantic_test/intermediate/doxy_profiles")
}

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(sprof_dir, pattern = "_Sprof.nc", full.names = TRUE)
message("Found ", length(files), " Sprof files in ", sprof_dir)

pb <- txtProgressBar(min = 0, max = length(files), style = 3)
i  <- 0

for (file in files) {

  i <- i + 1
  setTxtProgressBar(pb, i)
  nc <- nc_open(file)

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
  wmo    <- gsub(" ", "", ncvar_get(nc, "PLATFORM_NUMBER")[1])

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
      ~ if (sum(!is.na(.x)) > 10) {
          na.approx(.x, x = depth, na.rm = FALSE, rule = 2)
        } else {
          rep(NA_real_, length(.x))
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

  # save
  filename <- file.path(out_dir, paste0("argo_", wmo, "_interp.csv"))
  write_csv(df_interp, filename)
  message(filename, " written")
}

close(pb)
message("Done. Processed ", i, " floats.")

# Signal completion to Snakemake by writing a manifest
manifest <- data.frame(
  wmo  = gsub(".*argo_(.*)_interp.csv", "\\1",
               list.files(out_dir, pattern = "_interp.csv")),
  path = list.files(out_dir, pattern = "_interp.csv", full.names = TRUE)
)
write_csv(manifest, file.path(out_dir, "processing_manifest.csv"))