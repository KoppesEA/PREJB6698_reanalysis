#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J BSgen-%A_%a
#SBATCH --output=PRJEB6698.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g 

GenDir=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus/GRCm39_ref

module load bismark/0.20.0

bismark_genome_preparation

USAGE: bismark_genome_preparation \
--bowtie2



[options] <arguments>


OPTIONS:

--help/--man             Displays this help filea and exits.

--version                Displays version information and exits.

--verbose                Print verbose output for more details or debugging.

--path_to_bowtie </../>  The full path to the Bowtie 1 or Bowtie 2 installation on your system
                         (depending on which aligner/indexer you intend to use). Unless this path
                         is specified it is assumed that Bowtie is in the PATH.

--bowtie2                This will create bisulfite indexes for Bowtie 2. (Default: ON).

--bowtie1                This will create bisulfite indexes for Bowtie 1. (Default: OFF).

--single_fasta           Instruct the Bismark Indexer to write the converted genomes into
                         single-entry FastA files instead of making one multi-FastA file (MFA)
                         per chromosome. This might be useful if individual bisulfite converted
                         chromosomes are needed (e.g. for debugging), however it can cause a
                         problem with indexing if the number of chromosomes is vast (this is likely
                         to be in the range of several thousand files; the operating system can
                         only handle lists up to a certain length, and some newly assembled
                         genomes may contain 20000-500000 contigs of scaffold files which do exceed
                         this list length limit).

--genomic_composition    Calculate and extract the genomic sequence composition for mono and di-nucleotides
                         and write the genomic composition table 'genomic_nucleotide_frequencies.txt' to the
                         genome folder. This may be useful later on when using bam2nuc or the Bismark option
                         --nucleotide_coverage.


ARGUMENTS:

<path_to_genome_folder>  The path to the folder containing the genome to be bisulfite converted.
                         The Bismark Genome Preparation expects one or more fastA files in the folder
                         (with the file extension: .fa or .fasta (also ending in .gz)). Specifying this path is mandatory.


This script was last modified on 26 April 2018.
