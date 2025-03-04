#!/bin/bash
#SBATCH --job-name=first_annotation
#SBATCH --output=/data/CARD_MPU/users/trujilloae/CARD_unified_workflow/annotation.txt
#SBATCH --error=/data/CARD_MPU/users/trujilloae/CARD_unified_workflow/annotation.log
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=6        # Number of threads
#SBATCH --mem=2000GB                 # Memory per node
#SBATCH --partition=largemem
#SBATCH --time=48:00:00          # Time limit hrs:min:sec



module purge
module load apptainer
module load snakemake/7.7.0

# Pull profile, this will only run once, and is required for running on Biowulf
git clone https://github.com/NIH-HPC/snakemake_profile.git

# Pull the containers
apptainer pull envs/single_cell_cpu.sif oras://quay.io/adamcatchingdti/single_cell_cpu

# Load singularity
module load singularity/4.1.5

# Bind external directories on Biowulf
. /usr/local/current/singularity/app_conf/sing_binds

# RUN SCRIPT
snakemake --cores all --profile snakemake_profile --use-singularity