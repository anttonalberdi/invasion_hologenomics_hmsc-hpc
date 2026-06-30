#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Hmsc))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript make_cv_models.R <model_file> <cv_dir> <nfolds> [cv_column] [seed]",
  "",
  "Examples:",
  "  Rscript make_cv_models.R hmsc/model_red  cv/red  5",
  "  Rscript make_cv_models.R hmsc/model_grey cv/grey 5",
  "  Rscript make_cv_models.R hmsc/model_red  cv/red  5 sample 1",
  sep = "\n"
)

if (length(args) < 3) stop(usage, call. = FALSE)

model_file <- args[[1]]
cv_dir     <- args[[2]]
nfolds     <- as.integer(args[[3]])
cv_column  <- if (length(args) >= 4) args[[4]] else NA_character_
seed       <- if (length(args) >= 5) as.integer(args[[5]]) else 1L

if (is.na(nfolds) || nfolds < 2) stop("<nfolds> must be an integer >= 2", call. = FALSE)
if (!file.exists(model_file)) stop(sprintf("Model file not found: %s", model_file), call. = FALSE)

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

save_as_loadable_m <- function(hM, file) {
  m <- hM
  save(m, file = file)
}

make_train_model <- function(hM, train) {
  dfPi <- as.data.frame(
    matrix(NA, sum(train), hM$nr),
    stringsAsFactors = TRUE
  )
  colnames(dfPi) <- hM$rLNames

  for (r in seq_len(hM$nr)) {
    dfPi[, r] <- factor(hM$dfPi[train, r])
  }

  XTrain <- switch(
    class(hM$X)[1L],
    matrix = hM$X[train, , drop = FALSE],
    list = lapply(hM$X, function(a) a[train, , drop = FALSE]),
    stop(sprintf("Unsupported hM$X class: %s", class(hM$X)[1L]), call. = FALSE)
  )

  XRRRTrain <- if (hM$ncRRR > 0) {
    hM$XRRR[train, , drop = FALSE]
  } else {
    NULL
  }

  hM1 <- Hmsc(
    Y = hM$Y[train, , drop = FALSE],
    X = XTrain,
    XRRR = XRRRTrain,
    ncRRR = hM$ncRRR,
    XSelect = hM$XSelect,
    distr = hM$distr,
    studyDesign = dfPi,
    Tr = hM$Tr,
    C = hM$C,
    ranLevels = hM$rL
  )

  setPriors(
    hM1,
    V0 = hM$V0,
    f0 = hM$f0,
    mGamma = hM$mGamma,
    UGamma = hM$UGamma,
    aSigma = hM$aSigma,
    bSigma = hM$bSigma,
    nu = hM$nu,
    a1 = hM$a1,
    b1 = hM$b1,
    a2 = hM$a2,
    b2 = hM$b2,
    rhopw = hM$rhowp
  )

  hM1$YScalePar <- hM$YScalePar
  hM1$YScaled <- (
    hM1$Y - matrix(hM1$YScalePar[1, ], hM1$ny, hM1$ns, byrow = TRUE)
  ) / matrix(hM1$YScalePar[2, ], hM1$ny, hM1$ns, byrow = TRUE)

  hM1$XInterceptInd <- hM$XInterceptInd
  hM1$XScalePar <- hM$XScalePar

  switch(
    class(hM$X)[1L],
    matrix = {
      hM1$XScaled <- (
        hM1$X - matrix(hM1$XScalePar[1, ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
      ) / matrix(hM1$XScalePar[2, ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
    },
    list = {
      hM1$XScaled <- list()
      for (zz in seq_len(length(hM1$X))) {
        hM1$XScaled[[zz]] <- (
          hM1$X[[zz]] - matrix(hM1$XScalePar[1, ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
        ) / matrix(hM1$XScalePar[2, ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
      }
    }
  )

  if (hM1$ncRRR > 0) {
    hM1$XRRRScalePar <- hM$XRRRScalePar
    hM1$XRRRScaled <- (
      hM1$XRRR - matrix(hM1$XRRRScalePar[1, ], hM1$ny, hM1$ncORRR, byrow = TRUE)
    ) / matrix(hM1$XRRRScalePar[2, ], hM1$ny, hM1$ncORRR, byrow = TRUE)
  }

  hM1$TrInterceptInd <- hM$TrInterceptInd
  hM1$TrScalePar <- hM$TrScalePar
  hM1$TrScaled <- (
    hM1$Tr - matrix(hM1$TrScalePar[1, ], hM1$ns, hM1$nt, byrow = TRUE)
  ) / matrix(hM1$TrScalePar[2, ], hM1$ns, hM1$nt, byrow = TRUE)

  hM1
}

set.seed(seed)
model <- load_hmsc_model(model_file)

dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)

partition <- if (is.na(cv_column) || cv_column %in% c("", "NA", "NULL")) {
  createPartition(model, nfolds = nfolds)
} else {
  createPartition(model, nfolds = nfolds, column = cv_column)
}

saveRDS(model, file.path(cv_dir, "full_model.rds"))
saveRDS(partition, file.path(cv_dir, "partition.rds"))

fold_table <- data.frame(
  unit = seq_along(partition),
  fold = as.integer(partition)
)
write.table(
  fold_table,
  file = file.path(cv_dir, "fold_table.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

settings <- list(
  model_file = model_file,
  cv_dir = cv_dir,
  nfolds = nfolds,
  cv_column = cv_column,
  seed = seed,
  ny = model$ny,
  ns = model$ns,
  created_at = as.character(Sys.time())
)
saveRDS(settings, file.path(cv_dir, "settings_cv.rds"))

for (k in seq_len(nfolds)) {
  message(sprintf("Creating fold %d / %d", k, nfolds))
  train <- partition != k
  fold_model <- make_train_model(model, train)

  rds_path <- file.path(cv_dir, sprintf("model_fold_%02d.rds", k))
  rdata_path <- file.path(cv_dir, sprintf("model_fold_%02d", k))

  saveRDS(fold_model, rds_path)
  save_as_loadable_m(fold_model, rdata_path)
}

message("Done. Wrote CV fold models to: ", cv_dir)
