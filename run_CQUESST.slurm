#!/bin/bash
#SBATCH --job-name="CQUESST_MTT"
#SBATCH --time=168:00:00
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:0
#SBATCH -a 1
module load gcc
export STAN_NUM_THREADS=32
time ./CQUEST/CQUESST sample num_samples=50000 num_warmup=20000 save_warmup=1 thin=10 adapt engaged=1 algorithm=hmc stepsize=1.0 data file=data_cmdStan init=init_params random seed=1 output file=output.csv
