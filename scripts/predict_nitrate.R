library(tidyverse)
library(here)
library(seacarb)
library(mgcv)
library(data.table)
library(furrr)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  in_manifest   <- snakemake@input[["manifest"]]
  shared_dir    <- snakemake@params[["shared_dir"]]
  out_manifest  <- snakemake@output[["manifest"]]
  n_cores       <- snakemake@threads
  utils_dir     <- file.path(dirname(snakemake@scriptdir), "scripts", "utils")
  canyon_script <- file.path(utils_dir, "fastr_canyon.R")
  if (!file.exists(canyon_script)) stop("Cannot find fastr_canyon.R at: ", canyon_script)
} else {
  args         <- commandArgs(trailingOnly = TRUE)
  in_manifest  <- ifelse(length(args) > 0, args[1],
                         "data/NorthAtlantic_seas_comparison/intermediate/doxy_profiles/processing_manifest.csv")
  shared_dir   <- ifelse(length(args) > 1, args[2], "data/shared")
  out_manifest <- ifelse(length(args) > 2, args[3],
                         "data/NorthAtlantic_seas_comparison/intermediate/nitrate_profiles/nitrate_manifest.csv")
  n_cores      <- as.integer(Sys.getenv("NCORES", unset = max(1L, parallel::detectCores() - 1L)))
  utils_dir    <- "scripts/utils"
}

# ── Load CANYON-B ─────────────────────────────────────────────────────────────
source(file.path(utils_dir, "fastr_canyon.R"))
weights <- load_canyonb_weights(inputsdir = paste0(utils_dir, "/canyon_weights/"))

# ── Setup ─────────────────────────────────────────────────────────────────────
shared_nitrate_dir <- file.path(shared_dir, "nitrate_profiles")
dir.create(shared_nitrate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_manifest), recursive = TRUE, showWarnings = FALSE)

# Read file list from the doxy manifest
in_df <- readr::read_csv(in_manifest, show_col_types = FALSE)
files <- in_df$path
message("Found ", length(files), " interpolated float files in manifest")
message("Using ", n_cores, " parallel worker(s)")

# ── Parallel worker ───────────────────────────────────────────────────────────
process_float <- function(file, shared_nitrate_dir, weights, utils_dir) {

  # Re-source CANYON-B inside each worker (furrr multisession = separate R processes)
  source(file.path(utils_dir, "fastr_canyon.R"))
  dir.create(shared_nitrate_dir, recursive = TRUE, showWarnings = FALSE)

  wmo      <- gsub("argo_(.*)_interp\\.csv", "\\1", basename(file))
  out_path <- file.path(shared_nitrate_dir, paste0("argo_", wmo, "_nitrate.csv"))

  if (file.exists(out_path)) {
    message("Skipping ", wmo, " (nitrate already in shared)")
    return(out_path)
  }

  dat <- readr::read_csv(file, show_col_types = FALSE)

  required_cols <- c("date", "lat", "lon", "depth", "temp", "sal", "oxygen")
  if (!all(required_cols %in% names(dat))) {
    message("Skipping ", basename(file), " -- missing required columns.")
    return(NULL)
  }
  if (sum(!is.na(dat$oxygen)) < 10) {
    message("Skipping ", basename(file), " -- insufficient DOXY data.")
    return(NULL)
  }

  tryCatch({
    nitrate <- CANYONB_fast(
      date         = dat$date,
      lat          = dat$lat,
      lon          = dat$lon,
      pres         = dat$depth,
      temp         = dat$temp,
      psal         = dat$sal,
      doxy         = dat$oxygen,
      wgts_list    = weights,
      use_parallel = FALSE,   # within-float parallelism via mclapply does not work on Windows
      param        = "NO3"
    )
    dat$canyon_nitrate    <- nitrate$NO3
    dat$canyon_nitrate_ci <- nitrate$NO3_ci
  }, error = function(e) {
    message("CANYON-B failed for ", basename(file), ": ", e$message)
    dat$canyon_nitrate    <<- NA_real_
    dat$canyon_nitrate_ci <<- NA_real_
  })

  tryCatch({
    readr::write_csv(dat, out_path)
    out_path
  }, error = function(e) {
    message("Write failed for ", basename(file), ": ", e$message)
    NULL
  })
}

# ── Run in parallel ───────────────────────────────────────────────────────────
plan(multisession, workers = n_cores)

processed <- future_map(
  files,
  process_float,
  shared_nitrate_dir = shared_nitrate_dir,
  weights            = weights,
  utils_dir          = utils_dir,
  .options           = furrr_options(seed = TRUE),
  .progress          = TRUE
) |>
  compact() |>
  unlist()

plan(sequential)

message("Done. Processed ", length(processed), " floats.")

# ── Write manifest ────────────────────────────────────────────────────────────
manifest <- data.frame(
  wmo  = gsub("argo_(.*)_nitrate\\.csv", "\\1", basename(processed)),
  path = processed
)
readr::write_csv(manifest, out_manifest)
message("Manifest -> ", out_manifest)
