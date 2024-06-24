#!/bin/bash
#SBATCH --job-name="CQUESST_MTT_MPI"
#SBATCH --time=168:00:00
#SBATCH --mem=4GB
#SBATCH -N 4
#SBATCH --tasks-per-node=10
#SBATCH --gres=gpu:0
#SBATCH -a 1
module purge
module load slurm SC openmpi  gcc  libxml2 python
export STAN_NUM_THREADS=32
time mpirun -np 40 --bind-to none --mca orte_tmpdir_base <insert_your_path_here>/cmdstan ./CQUESST1/CQUESST sample num_samples=50000 num_warmup=20000 save_warmup=1 thin=10 adapt engaged=1 algorithm=hmc stepsize=1.0 data file=/scratch3/pag06e/cmdstan/CQUESST1/data_cmdStan init=<insert_your_path_here>/cmdstan/CQUESST1/init_params random seed=1 output file=<insert_your_path_here>/cmdstan/CQUESST1/output.csv
