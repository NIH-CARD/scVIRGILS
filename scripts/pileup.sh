#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --array=0-303
#SBATCH --gres=lscratch:2000
#SBATCH --output=/data/CARD_singlecell/SN_atlas/cellbender/logs/pileup-%a.out

module load samtools
module load bcftools

# Make a temporary location

mkdir /lscratch/$SLURM_JOB_ID/%a
export TMPDIR=/lscratch/$SLURM_JOB_ID

# Define the directories to run cell bender within
output_file_base=$(echo /data/CARD_singlecell/Brain_atlas/SN_Multiome/*/Multiome/*/outs/)
# Convert the locations of the directories to an array
out_dirs=(`echo ${output_file_base}`);

# Read in the sample sheet for each batch
bcftools mpileup -f /fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa ${out_dirs[$SLURM_ARRAY_TASK_ID]}"atac_possorted_bam.bam" | bcftools call -mv -Ov -o ${out_dirs[$SLURM_ARRAY_TASK_ID]}"variants.vcf"
#samtools mpileup -f /fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa ${out_dirs[$SLURM_ARRAY_TASK_ID]}"atac_possorted_bam.bam" > $TMPDIR"/output.pileup"
#bcftools call -mv -Ov $TMPDIR/"output.pileup" > ${out_dirs[$SLURM_ARRAY_TASK_ID]}"variants.vcf"