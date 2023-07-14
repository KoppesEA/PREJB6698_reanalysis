#!/bin/bash
#
##SBATCH --job-name=EKfastqc2
#SBATCH -n 1
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --output=fastqc.out

module load FastQC/0.11.5
mkdir pretrim_fastqc

fastqc -o ./pretrim_fastqc/ ./raw_fastq/120-1_subline_1_S1_R1_001.fastq.gz
fastqc -o ./pretrim_fastqc/ ./raw_fastq/120-1_subline_1_S1_R2_001.fastq.gz
fastqc -o ./pretrim_fastqc/ ./raw_fastq/120-1_subline_2_S2_R1_001.fastq.gz
fastqc -o ./pretrim_fastqc/ ./raw_fastq/120-1_subline_2_S2_R2_001.fastq.gz