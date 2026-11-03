#!/bin/bash
#$ -S /bin/bash
#$ -shell yes
#$ -N pca_lineprof
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

Rscript -e 'stopifnot(requireNamespace("lineprof", quietly = TRUE))'
dataset_dir="output/followup_genotyped_20261031_182624_492600"
test -s "$dataset_dir/data_g500.rds"
test -s "$dataset_dir/data_g4000.rds"
profile_dir="output/lineprof_${JOB_ID:?Run through qsub}"
mkdir "$profile_dir"
for count in 500 4000; do
  Rscript validation/profile_lineprof.R "$dataset_dir/data_g${count}.rds" \
    ../randPedGenPCA "$profile_dir/g${count}" FALSE 4 20 50 123
done
echo "Finished. lineprof profiles: $profile_dir"
