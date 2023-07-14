#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH -J BSMRK_ARRAY
#SBATCH --output=Bismark_MethExtract-%A_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --array=0-9 # 10 paired-end Samples 

#BASH SCRIPT TO EXTRACT METHYLATION FOR EACH POSITION
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
GENREF=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus/GRCm39_ref/
BAMDIR=$WKDIR/Bismark
OUTDIR=$BAMDIR

names=($(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt | grep "_1.fastq.gz"))
ERR=${names[${SLURM_ARRAY_TASK_ID}]}
echo "ERR _R1 Record is: " $ERR
ERR_base=`basename $ERR _1.fastq.gz`
echo "ERR Basename is: " $ERR_base
BAM_IN=${ERR_base}_1_val_1_bismark_bt2_pe.bam
echo "Bismark BAM file is : " $BAM_IN
echo "Starting Bismark Meth Extractor for BAM Record: " $ERR_base "using genome in: " $GENREF

module load bismark/0.20.0

rm -r $OUTDIR/${ERR_base}_methextract
mkdir $OUTDIR/${ERR_base}_methextract

bismark_methylation_extractor \
	--parallel 16 \
	--paired-end \
	--gzip \
	--output $OUTDIR/${ERR_base}_methextract \
	--bedGraph \
	--comprehensive \
	--cytosine_report \
	--genome_folder $GENREF \
	$BAMDIR/$BAM_IN
	