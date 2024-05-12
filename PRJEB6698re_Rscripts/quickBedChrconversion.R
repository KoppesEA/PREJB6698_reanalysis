## Script to import UCSC bed files and change from UCSC to ENSEMBL genome names to 1--base

## readr dependencies
library(readr)
library(dplyr)
library(stringr)

## Download UCSC CGI or other BED track files here: https://genome.ucsc.edu/cgi-bin/hgTables
## UCSC and Ensembl have different chromosome names
## UCSC and Ensembl use 0-based numbering for .bed files (don't need start +1 to convert)
CPG_island_bed <-"./GRCm39_Annot/GRCm39_UCSC_CGI.bed" #not formatted as bed, has "#bin" as column 1
CPGi_df <- read_tsv(CPG_island_bed) %>%
  select(chrom, chromStart, chromEnd, name, length, cpgNum, gcNum, perCpg, perGc, obsExp)

## Download GRCm39 assembly chromosome tsv file here: https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000001635.27/
## UCSC chrom names are chr1-19, chrX, chrY, chrMT plus UCSC style contigs
## Ensembl names are 1-19, X, Y, MT plus Genbank style contigs
GRCm39_chrNames <- read_tsv("GCF_000001635.27.tsv") %>%
  select(4, 7, 13)
names(GRCm39_chrNames) <- c("Chrom", "GenBank", "UCSC")

GRCm39_chrNames_conv <-bind_rows(
  GRCm39_chrNames[1:21, c(1,3)] %>% rename(Ensembl = Chrom),
  GRCm39_chrNames[22:61, c(2,3)] %>% rename(Ensembl = GenBank)
  )

## Inner Join CGI dataframe to GRCm39 name conversion df then select Ensembl names and cols 2-4, then rename 1st col back to chrom
CPGi_df_rename <- left_join(CPGi_df, GRCm39_chrNames_conv, # use inner or left join to filter to chr1-19 +X +Y (no MT, contigs, alts)
          join_by(chrom == UCSC)) %>%
  select(Ensembl, chromStart, chromEnd, name, length, cpgNum, gcNum, perCpg, perGc, obsExp) %>%
  rename(chrom = Ensembl)

## Export as tsv with .bed file name
write_tsv(CPGi_df_rename, "GRCm39_CGI_Ensembl.bed")

