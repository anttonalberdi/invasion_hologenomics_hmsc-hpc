#!/bin/bash
#SBATCH --job-name=hmsc_cv_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=32gb
#SBATCH --time=120:00:00
# Set the array at submission, e.g. --array=0-19 for 5 folds x 4 chains.

set -euo pipefail

if [ "$#" -lt 4 ]; then
  echo "Usage: sbatch --array=0-<nfolds*nchains-1> run_cv_hmsc_gpu.sh <cv_dir> <post_dir> <samples> <thin> [nchains] [transient] [verbose]" >&2
  echo "Example: sbatch --array=0-19 run_cv_hmsc_gpu.sh cv/red output_cv/red 250 1 4" >&2
  exit 1
fi

CV_DIR="$1"
POST_DIR="$2"
SAMPLES="$3"
THIN="$4"
NCHAINS="${5:-4}"
TRANSIENT="${6:-$(( (SAMPLES * THIN + 1) / 2 ))}"
VERBOSE="${7:-100}"

mkdir -p "$POST_DIR"

NFOLDS=$(find "$CV_DIR" -maxdepth 1 -name 'init_fold_*.rds' | wc -l | tr -d ' ')
if [ "$NFOLDS" -lt 1 ]; then
  echo "No init_fold_*.rds files found in $CV_DIR. Run make_cv_hpc_inits.R first." >&2
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
FOLD=$(( TASK_ID / NCHAINS + 1 ))
CHAIN=$(( TASK_ID % NCHAINS ))

if [ "$FOLD" -gt "$NFOLDS" ]; then
  echo "Task $TASK_ID maps to fold $FOLD, but only $NFOLDS folds exist. Exiting." >&2
  exit 0
fi

FOLD2=$(printf "%02d" "$FOLD")
CHAIN2=$(printf "%02d" "$CHAIN")

INPUT_PATH="${CV_DIR}/init_fold_${FOLD2}.rds"
OUTPUT_PATH="${POST_DIR}/fold_${FOLD2}_chain_${CHAIN2}.rds"

if [ ! -f "$INPUT_PATH" ]; then
  echo "Missing input: $INPUT_PATH" >&2
  exit 1
fi

echo "Running fold ${FOLD2}, chain ${CHAIN2}"

srun python3 -m hmsc.run_gibbs_sampler \
  --input "$INPUT_PATH" \
  --output "$OUTPUT_PATH" \
  --samples "$SAMPLES" \
  --transient "$TRANSIENT" \
  --thin "$THIN" \
  --verbose "$VERBOSE" \
  --chain "$CHAIN"
