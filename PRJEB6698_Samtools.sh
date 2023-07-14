#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH -J SAM_ARRAY
#SBATCH --output=SAM-RIP-%A_%a.txt
#SBATCH --cpus-per-task=8
#SBATCH --mem=12g
#SBATCH --array=0-9 # 10 paired-end Samples 

#BASH SCRIPT TO HAVE PRETRIM QC OF RRBS DATA FROM PRJEB6698
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
BAMDIR=$WKDIR/Bismark

module load gcc/8.2.0
module load samtools/1.14

names=($(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt | grep "_1.fastq.gz"))
ERR=${names[${SLURM_ARRAY_TASK_ID}]}
echo "ERR _R1 Record is: " $ERR
ERR_base=`basename $ERR _1.fastq.gz`
echo "ERR Basename is: " $ERR_base
BAM_IN=${ERR_base}_1_val_1_bismark_bt2_pe.bam
echo "Bismark BAM file is : " $BAM_IN
BAM_OUT=${ERR_base}_1_val_1_bismark_bt2_pe_sorted.bam

samtools sort \
	-@ 8 \
	-o $BAMDIR/$BAM_OUT \
	-O BAM \
	--write-index \
	$BAMDIR/$BAM_IN
done

