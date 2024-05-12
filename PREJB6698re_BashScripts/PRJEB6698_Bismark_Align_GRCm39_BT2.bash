#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH -J BSMRK_ARRAY
#SBATCH --output=Bismark_PRJEB6698-%A_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --array=0-9 # 10 paired-end Samples 

#BASH SCRIPT TO Align Trimmed RRBS DATA FROM PRJEB6698
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
TRIMDIR=$WKDIR/fastq_trim
GENREF=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus/GRCm39_ref/
echo "BT2 generated genome directory is var GENREF: " $GENREF 
OUTDIR=$WKDIR/Bismark

names=($(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt | grep "_1.fastq.gz"))
ERR=${names[${SLURM_ARRAY_TASK_ID}]}
echo "ERR _R1 Record is: " $ERR
ERR_base=`basename $ERR _1.fastq.gz`
echo "ERR Basename is: " $ERR_base
TRIM_FASTQ1=${ERR_base}_1_val_1.fq.gz
echo "ERR RRBS Trimmed Fastq R1 Record is: " $TRIM_FASTQ1
TRIM_FASTQ2=${ERR_base}_2_val_2.fq.gz
echo "ERR RRBS Trimmed Fastq R2 Record is: " $TRIM_FASTQ2
echo "Starting Bismark for ERR Record: " $ERR_base
	
module load bismark/0.20.0

bismark $GENREF -1 ${TRIMDIR}/$TRIM_FASTQ1 -2 ${TRIMDIR}/$TRIM_FASTQ2
