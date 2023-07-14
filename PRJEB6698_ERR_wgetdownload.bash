#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J PBJ-TIME
#SBATCH --output=PRJEB6698.txt
#SBATCH --cpus-per-task=1

#BASH SCRIPT TO WGET ENA SEQ DATA FROM PRJEB6698
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
OUTDIR=$WKDIR/ERR_fastq_download

for ERR in $(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt)
do
	echo $ERR
	ERR_base=`basename $ERR .fastq.gz`
	echo $ERR_base
	wget -nv -O $OUTDIR/${ERR_base}.fastq.gz $ERR
done