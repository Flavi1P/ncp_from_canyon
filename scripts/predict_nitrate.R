library(tidyverse)
library(here)
library(seacarb)
library(mgcv)
library(data.table)
library(furrr)

# ── Snakemake / CLI input handling ────────────────────────────────────────────
if (exists("snakemake")) {
  in_dir        <- snakemake@input[["in_dir"]]
  out_dir       <- snakemake@output[["out_dir"]]
  n_cores       <- snakemake@threads
  utils_dir     <- file.path(dirname(snakemake@scriptdir), "scripts", "utils")
  canyon_script <- file.path(utils_dir, "fastr_canyon.R")
  if (!file.exists(canyon_script)) stop("Cannot find fastr_canyon.R at: ", canyon_script)
} else {
  args      <- commandArgs(trailingOnly = TRUE)
  in_dir    <- ifelse(length(args) > 0, args[1],
                      "data/NorthAtlantic_2024/intermediate/doxy_profiles")
  out_dir   <- ifelse(length(args) > 1, args[2],
                      "data/NorthAtlantic_2024/intermediate/nitrate_profiles")
  n_cores   <- as.integer(Sys.getenv("NCORES", unset = max(1L, parallel::detectCores() - 1L)))
  utils_dir <- "scripts/utils"
}

# ── Load CANYON-B ─────────────────────────────────────────────────────────────
source(file.path(utils_dir, "fastr_canyon.R"))
weights <- load_canyonb_weights(inputsdir = paste0(utils_dir, "/canyon_weights/"))

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "_interp\\.csv$", full.names = TRUE)
message("Found ", length(files), " interpolated float files in ", in_dir)
message("Using ", n_cores, " parallel worker(s)")

# ── Parallel worker ───────────────────────────────────────────────────────────
process_float <- function(file, out_dir, weights, utils_dir) {

  # Re-source CANYON-B inside each worker (furrr multisession = separate R processes)
  source(file.path(utils_dir, "fastr_canyon.R"))

  # Ensure output directory exists in this worker's context
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    dat$canyon_nitrate <- nitrate$NO3
  }, error = function(e) {
    message("CANYON-B failed for ", basename(file), ": ", e$message)
    dat$canyon_nitrate <<- NA_real_
  })

  wmo      <- gsub("argo_(.*)_interp\\.csv", "\\1", basename(file))
  out_path <- file.path(out_dir, paste0("argo_", wmo, "_nitrate.csv"))
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
  out_dir   = out_dir,
  weights   = weights,
  utils_dir = utils_dir,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
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
readr::write_csv(manifest, file.path(out_dir, "nitrate_manifest.csv"))
message("Manifest -> ", file.path(out_dir, "nitrate_manifest.csv"))
