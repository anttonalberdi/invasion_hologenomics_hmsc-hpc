# hmsc.R
library(Hmsc)
library(jsonify)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript hmsc.R <unfitted_model_path> <output_init_path>")
}

model_path <- args[1]
init_path <- args[2]

# Pass shell variables into R, use Sys.getenv()
samples <- as.integer(Sys.getenv("samples"))
thin <- as.integer(Sys.getenv("thin"))
transient <- as.integer(Sys.getenv("transient"))
nChains <- as.integer(Sys.getenv("nChains"))
verbose <- as.logical(Sys.getenv("verbose"))
cInd <- as.integer(Sys.getenv("cInd"))
unfitted_model <- Sys.getenv("unfitted_model")

# Load unfitted model
load(file = model_path)

# Run init_obj
init_obj <- sampleMcmc(m,
              samples=samples,
              thin=thin,
              transient=transient,
              nChains=nChains,
              verbose=verbose,
              engine="HPC")

# Save init_obj
saveRDS(init_obj, file=init_path)
