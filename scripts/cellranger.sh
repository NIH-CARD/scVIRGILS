#!/bin/bash

#SBATCH --array=0-20  # Adjust this based on the number of samples
#SBATCH --gres=lscratch:2000  # Request 2000 MB of local scratch space
#SBATCH --mem=64G  # Request the actual memory needed per task
#SBATCH --cpus-per-task=8  # Allocate CPU resources
#SBATCH --time=48:00:00  # Request 48 hours of runtime
#SBATCH --output=/data/CARD_MPU/users/veeraraghavank/marsan_2023/cellranger/logs/cellranger-%a.out  # Log file output
#SBATCH --error=/data/CARD_MPU/users/veeraraghavank/marsan_2023/cellranger/logs/cellranger-%a.err  # Error log output

# Load Cell Ranger
module load cellranger/8.0.1

# Define an array of sample names
samples=("CTRL_1" "CTRL_2" "CTRL_3" "CTRL_4" "CTRL_5" "CTRL_6" "CTRL_7" "CTRL_8" "CTRL_9" "CTRL_10"
"FTD_1" "FTD_2" "FTD_3" "FTD_4" "FTD_5" "FTD_6" "FTD_7" "FTD_8" "FTD_9")

# Set variables for the current sample and directories
sample=${samples[$SLURM_ARRAY_TASK_ID]}
fastq_dir="set/working/directory/${sample}"
output_dir="set/working/directory//${sample}_output"
scratch_dir="set/working/directory/lscratch/$SLURM_JOB_ID/${sample}_count"

# Ensure scratch directory exists
mkdir -p $scratch_dir

# Run Cell Ranger count, directing outputs to scratch
cellranger count --id="${sample}_count" \
  --fastqs="$fastq_dir" \
  --transcriptome=/fdb/cellranger/refdata-gex-GRCh38-2024-A \
  --localcores=8 \
  --localmem=64 \
  --create-bam true \
  --output-dir="$scratch_dir"

# Verify and move output files to the final directory
if [[ -f "$scratch_dir/outs/filtered_feature_bc_matrix.h5" && -f "$scratch_dir/outs/possorted_genome_bam.bam" && -f "$scratch_dir/outs/possorted_genome_bam.bam.bai" ]]; then
    mkdir -p "$output_dir"
    mv "$scratch_dir/outs/filtered_feature_bc_matrix.h5" "$output_dir/"
    mv "$scratch_dir/outs/possorted_genome_bam.bam" "$output_dir/"
    mv "$scratch_dir/outs/possorted_genome_bam.bam.bai" "$output_dir/"
else
    echo "Error: Expected output files not found for sample $sample"
fi
