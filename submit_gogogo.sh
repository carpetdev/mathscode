#!/bin/bash
#SBATCH --job-name=potts12   # Job name
#SBATCH --time=10:00:00         # Request runtime (hh:mm:ss)
#SBATCH --mem=16G                # Request memory
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=100       # Number of cores per task

# Load any necessary modules
module load julia

# Execute your application
julia -t 100 gogogo.jl