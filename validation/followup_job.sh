#!/bin/bash
#$ -S /bin/bash
#$ -shell yes
#$ -cwd
#$ -j y
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G

source /etc/profile
set -eo pipefail
module load R/4.4
export R_LIBS_USER="$HOME/HighlanderLab/jortiz/pca-benchmark/R-library/4.4"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# qsub this script from the benchmarks directory, passing exact or genotyped.
experiment="${1:?Specify exact or genotyped}"
shift
python3 validation/run_followups.py "$experiment" \
  --package-root ../randPedGenPCA "$@"
