#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Hmsc))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript make_cv_hpc_inits.R <cv_dir> <samples> <thin> <nchains> [transient] [adapt_fraction] [verbose]",
  "",
  "Example:",
  "  Rscript make_cv_hpc_inits.R cv/red 250 1 4",
  sep = "\n"
)

if (length(args) < 4) stop(usage, call. = FALSE)

cv_dir <- args[[1]]
samples <- as.integer(args[[2]])
thin <- as.integer(args[[3]])
nchains <- as.integer(args[[4]])
transient <- if (length(args) >= 5) as.integer(args[[5]]) else ceiling(0.5 * samples * thin)
adapt_fraction <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.4
verbose <- if (length(args) >= 7) as.integer(args[[7]]) else 100L

if (!dir.exists(cv_dir)) stop(sprintf("cv_dir not found: %s", cv_dir), call. = FALSE)
if (any(is.na(c(samples, thin, nchains, transient))) || samples < 1 || thin < 1 || nchains < 1) {
  stop("samples, thin, nchains and transient must be positive integers", call. = FALSE)
}

settings_path <- file.path(cv_dir, "settings_cv.rds")
if (!file.exists(settings_path)) stop("settings_cv.rds not found. Run make_cv_models.R first.", call. = FALSE)
settings <- readRDS(settings_path)
nfolds <- settings$nfolds

settings$samples <- samples
settings$thin <- thin
settings$nchains <- nchains
settings$transient <- transient
settings$adapt_fraction <- adapt_fraction
settings$verbose <- verbose
saveRDS(settings, settings_path)

for (k in seq_len(nfolds)) {
  message(sprintf("Initializing fold %d / %d for Hmsc-HPC", k, nfolds))

  model_path <- file.path(cv_dir, sprintf("model_fold_%02d.rds", k))
  if (!file.exists(model_path)) stop(sprintf("Fold model missing: %s", model_path), call. = FALSE)

  hM <- readRDS(model_path)
  adaptNf <- rep(ceiling(adapt_fraction * samples * thin), hM$nr)

  init_obj <- sampleMcmc(
    hM = hM,
    samples = samples,
    thin = thin,
    adaptNf = adaptNf,
    transient = transient,
    nChains = nchains,
    verbose = verbose,
    engine = "HPC"
  )

  saveRDS(init_obj, file.path(cv_dir, sprintf("init_fold_%02d.rds", k)))
}

message("Done. Wrote Hmsc-HPC init files to: ", cv_dir)
