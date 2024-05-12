#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J pretrimQC
#SBATCH --output=pretrimQC_PRJEB6698-%A_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=12g 

#BASH SCRIPT TO HAVE PRETRIM QC OF RRBS DATA FROM PRJEB6698
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
FQDIR=$WKDIR/ERR_fastq_download
OUTDIR=$WKDIR/pretrim_fastqc

rm -r $OUTDIR
mkdir $OUTDIR

module load fastqc/0.11.9

for ERR in $(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt)
do
	echo $ERR
	ERR_base=`basename $ERR`
	echo $ERR_base
	fastqc -o $OUTDIR ${FQDIR}/$ERR_base
done