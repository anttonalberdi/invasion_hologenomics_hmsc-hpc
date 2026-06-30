# Invasion hologenomics using hmsc-hpc

With crossvalidation using GPU

## Create fold-specific Hmsc models

- hmsc/model_red     Path to the original Hmsc model file
- cv/red             Output directory where CV fold models will be saved
- 5                  Number of cross-validation folds

```{sh}
conda activate hmsc-hpc
Rscript scripts/make_cv_models.R hmsc/model_red  cv/red  5
Rscript scripts/make_cv_models.R hmsc/model_grey cv/grey 5
```

## Initialize fold models

- cv/red             Directory containing the CV fold models
- 250                Number of posterior samples per chain
- 100                  Thinning interval
- 4                  Number of MCMC chains

```{sh}
Rscript scripts/make_cv_hpc_inits.R cv/red  250 100 4
Rscript scripts/make_cv_hpc_inits.R cv/grey 250 100 4
```

## Run GPU jobs

- cv/red             Directory containing the Hmsc-HPC init files
- output_cv/red      Directory where GPU posterior outputs will be saved
- 250                Number of posterior samples per chain
- 100                Thinning interval
- 4                  Number of MCMC chains

```{sh}
sbatch --array=0-19 scripts/run_cv_hmsc_gpu.sh cv/red  output_cv/red  250 100 4
sbatch --array=0-19 scripts/run_cv_hmsc_gpu.sh cv/grey output_cv/grey 250 100 4
```

## Collect the fold predictions using the GPU-generated fold posteriors

- cv/grey            Directory containing the Hmsc-HPC init files
- output_cv/grey     Directory where GPU posterior outputs will be saved
- 250                Number of posterior samples per chain
- 100                  Thinning interval
- 4                  Number of MCMC chains

```{sh}
Rscript scripts/collect_cv_predictions.R hmsc/model_red  cv/red  output_cv/red  hmsc/cv_model_red_250_1.rds  250 100 4
Rscript scripts/collect_cv_predictions.R hmsc/model_grey cv/grey output_cv/grey hmsc/cv_model_grey_250_1.rds 250 100 4
```