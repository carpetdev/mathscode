#!/bin/bash
#SBATCH --job-name=simple_job   # Job name
#SBATCH --time=01:00:00         # Request runtime (hh:mm:ss)
#SBATCH --mem=1G                # Request memory
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=1       # Number of cores per task

# Load any necessary modules
module load <module_name>

# Execute your application
./example.bin