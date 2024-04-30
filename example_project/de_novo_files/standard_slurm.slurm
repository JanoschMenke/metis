#!/bin/bash
#SBATCH -N 1 --gpus-per-node=T4:1 
#SBATCH -t 0-23:00:00
#SBATCH --output={slurm_path}/{run_name}_slurm.out

source activate reinvent
rm -rf {slurm_path}/results/*
python /path/to/input.py {slurm_path}/new_run.json