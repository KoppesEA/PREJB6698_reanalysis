#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J RRBStrimQC
#SBATCH --output=RRBStrimQC_PRJEB6698-%A_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=12g 

#BASH SCRIPT TO HAVE PRETRIM QC OF RRBS DATA FROM PRJEB6698
WKDIR=/ix1/mmann/KoppesEA/PRJEB6698
FQDIR=$WKDIR/ERR_fastq_download
OUTDIR=$WKDIR/fastq_trim

rm -r $OUTDIR
mkdir $OUTDIR

module load trimgalore/0.6.5

for ERR in $(cat $WKDIR/ERR_Acc_List_PRJEB6698.txt | grep "_1.fastq.gz") #filter for _1 record paired end then get basename
do
	echo $ERR
	ERR_base=`basename $ERR _1.fastq.gz`
	echo $ERR_base
	INPUT_FASTQ1=${ERR_base}_1.fastq.gz
	echo $INPUT_FASTQ1
	INPUT_FASTQ2=${ERR_base}_2.fastq.gz
	echo $INPUT_FASTQ2
	trim_galore \
	--rrbs \
	--paired \
	--output_dir $OUTDIR \
	--fastqc \
	--gzip \
	${FQDIR}/$INPUT_FASTQ1 ${FQDIR}/$INPUT_FASTQ2
done