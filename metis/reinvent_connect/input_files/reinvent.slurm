#!/bin/bash
#SBATCH --job-name=myTest           # Job name
#SBATCH --account=project_2008745   # Billing project, has to be defined!
#SBATCH --time=02:00:00             # Max. duration of the job
#SBATCH --mem-per-cpu=16G            # Memory to reserve per core
#SBATCH --partition=small           # Job queue (partition)
#SBATCH --output=/users/nahalyas/Test_my_code/test_slurm.out

module load tykky
export PATH=/projappl/project_2008745/nahalyas/tykky-reinvent-hitl/bin:$PATH
rm -rf /users/nahalyas/Test_my_code/results/*
python /users/nahalyas/Test_my_code/Reinvent/input.py /users/nahalyas/Test_my_code/new_run.json