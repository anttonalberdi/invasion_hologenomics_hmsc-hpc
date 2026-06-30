#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Hmsc)
  library(abind)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript collect_cv_predictions.R <model_file|cv_dir/full_model.rds> <cv_dir> <post_dir> <output_file> <samples> <thin> <nchains> [transient] [expected]",
  "",
  "Examples:",
  "  Rscript collect_cv_predictions.R hmsc/model_red cv/red output_cv/red hmsc/cv_red.rds 250 1 4",
  "  Rscript collect_cv_predictions.R hmsc/model_grey cv/grey output_cv/grey hmsc/cv_grey.rds 250 1 4",
  sep = "\n"
)

if (length(args) < 7) stop(usage, call. = FALSE)

model_file <- args[[1]]
cv_dir <- args[[2]]
post_dir <- args[[3]]
output_file <- args[[4]]
samples <- as.integer(args[[5]])
thin <- as.integer(args[[6]])
nchains <- as.integer(args[[7]])
transient <- if (length(args) >= 8) as.integer(args[[8]]) else ceiling(0.5 * samples * thin)
expected <- if (length(args) >= 9) as.logical(args[[9]]) else TRUE

if (!file.exists(model_file)) stop(sprintf("Model file not found: %s", model_file), call. = FALSE)
if (!dir.exists(cv_dir)) stop(sprintf("cv_dir not found: %s", cv_dir), call. = FALSE)
if (!dir.exists(post_dir)) stop(sprintf("post_dir not found: %s", post_dir), call. = FALSE)

load_hmsc_model <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Hmsc")) stop("RDS file does not contain an Hmsc object", call. = FALSE)
    return(obj)
  }

  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)

  preferred <- c("m", "model", "hM")
  for (p in preferred) {
    if (p %in% nm && inherits(get(p, envir = e), "Hmsc")) return(get(p, envir = e))
  }

  candidates <- nm[vapply(nm, function(x) inherits(get(x, envir = e), "Hmsc"), logical(1))]
  if (length(candidates) == 1) return(get(candidates[[1]], envir = e))

  stop(
    sprintf(
      "Could not identify a single Hmsc object in %s. Objects found: %s",
      path,
      paste(nm, collapse = ", ")
    ),
    call. = FALSE
  )
}

make_validation_design <- function(hM, val) {
  dfPi <- as.data.frame(
    matrix(NA, sum(val), hM$nr),
    stringsAsFactors = TRUE
  )
  colnames(dfPi) <- hM$rLNames

  for (r in seq_len(hM$nr)) {
    dfPi[, r] <- factor(hM$dfPi[val, r])
  }

  dfPi
}

read_chain_list <- function(paths) {
  chainList <- vector("list", length(paths))

  for (i in seq_along(paths)) {
    if (!file.exists(paths[[i]])) stop(sprintf("Missing posterior file: %s", paths[[i]]), call. = FALSE)
    chainList[[i]] <- readRDS(paths[[i]])[["list"]][[1]]
  }

  chainList
}

model <- load_hmsc_model(model_file)
partition <- readRDS(file.path(cv_dir, "partition.rds"))
settings <- readRDS(file.path(cv_dir, "settings_cv.rds"))
nfolds <- settings$nfolds

if (length(partition) != model$ny) {
  stop("partition length does not match model$ny", call. = FALSE)
}

postN <- samples * nchains
cv <- array(NA_real_, dim = c(model$ny, model$ns, postN))

for (k in seq_len(nfolds)) {
  message(sprintf("Collecting CV fold %d / %d", k, nfolds))

  val <- partition == k
  fold_model_path <- file.path(cv_dir, sprintf("model_fold_%02d.rds", k))
  fold_model <- readRDS(fold_model_path)

  fold_paths <- file.path(
    post_dir,
    sprintf("fold_%02d_chain_%02d.rds", k, 0:(nchains - 1))
  )

  chainList <- read_chain_list(fold_paths)

  fold_fit <- importPosteriorFromHPC(
    fold_model,
    chainList,
    samples,
    thin,
    transient
  )

  postList <- poolMcmcChains(fold_fit$postList)

  XVal <- switch(
    class(model$X)[1L],
    matrix = model$X[val, , drop = FALSE],
    list = lapply(model$X, function(a) a[val, , drop = FALSE]),
    stop(sprintf("Unsupported model$X class: %s", class(model$X)[1L]), call. = FALSE)
  )

  XRRRVal <- if (model$ncRRR > 0) {
    model$XRRR[val, , drop = FALSE]
  } else {
    NULL
  }

  dfPiVal <- make_validation_design(model, val)

  pred <- predict(
    fold_fit,
    post = postList,
    X = XVal,
    XRRR = XRRRVal,
    studyDesign = dfPiVal,
    Yc = NULL,
    mcmcStep = 1,
    expected = expected
  )

  cv[val, , ] <- abind(pred, along = 3)
}

partition_cv <- partition

if (grepl("\\.rds$", output_file, ignore.case = TRUE)) {
  saveRDS(cv, output_file)
} else {
  save(cv, partition_cv, file = output_file)
}

message("Done. Wrote CV predictions to: ", output_file)
