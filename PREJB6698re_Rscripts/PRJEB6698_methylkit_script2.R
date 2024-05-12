## MethylKit (re)-Analysis of DNMT1 Tet-OFF (Mcgraw et al. 2015, NAR)
## MethylKit Script2 for 100bp-windowed analysis
## Erik Koppes, PhD
## Magee Womens Research Institute
## March 2024


## Dependencies
library(readr)
library(methylKit) 
packageVersion("methylKit") #v1.26.0 on CRC

## make a file list

MAIN_DIR <- "/ix1/mmann/KoppesEA/PRJEB6698"
BIS_DIR <- paste0(MAIN_DIR, "/Bismark")
METH_DIR <- paste0(MAIN_DIR, "/methylKit") ##wd

## Make a meta data table with PRJEB6698 Sample_names and Metadata
PRJEB6698_list <- c("ERR560526", "ERR560527", "ERR560528", "ERR560529", "ERR560530",
                    "ERR560531", "ERR560532", "ERR560533", "ERR560534", "ERR560538")
DIR_list <- as.list(paste0(BIS_DIR, "/", PRJEB6698_list, "_methextract/"))
PRJEB6698_Covfile_list <- as.list(paste0(PRJEB6698_list, "_1_val_1_bismark_bt2_pe.bismark.cov.gz"))
PRJEB6698_CpGreport_list <- as.list(paste0(PRJEB6698_list, "_1_val_1_bismark_bt2_pe.CpG_report.txt.gz"))
PRJEB6698_Covfile_list_loc <- as.list(paste0(BIS_DIR, "/", PRJEB6698_list, "_methextract/",PRJEB6698_list, "_1_val_1_bismark_bt2_pe.bismark.cov.gz"))

##Import metadata TSV (see script to obtain from .xml)
metaData <- read_tsv(paste0(MAIN_DIR, "/PRJEB6698_metadata.tsv"),
                     col_names = c("ERR_ID", "SAMEA_ID", "Treatment"),
                     col_select = 1:3)

## Import of low coverage, min=3 reads per CpG, for 100bp window calc
PRJEB6698_MethObj_lowCov = methRead(location = PRJEB6698_Covfile_list_loc,
                                    sample.id = as.list(metaData$Treatment),
                                    assembly="GRCm39",
                                    pipeline = "bismarkCoverage",
                                    treatment = c(0,0,1,1,1,2,2,2,2,0), #0 is control no-dox, 1 is d0 , 2 is d21 
                                    sep = "\t",
                                    context="CpG",
                                    mincov = 3
)

## 100bp window tiling --takes a few minutes---
PRJEB6698_MethObj_100win = tileMethylCounts(PRJEB6698_MethObj_lowCov,win.size=100,step.size=100,cov.bases = 10)
head(PRJEB6698_MethObj_100win[[1]],3)

## Filtering
PRJEB6698_MethObj_100win <- filterByCoverage(PRJEB6698_MethObj_100win,
                                          lo.count=3, ##already done with input
                                          lo.perc=NULL,
                                          hi.count=NULL,
                                          hi.perc=99.9) ##upperlimit cutoff implemented
## Normalization
PRJEB6698_MethObj_100win_norm <- normalizeCoverage(PRJEB6698_MethObj_100win, method = "median")

## Merge Data
PRJEB6698_MethObj_100win_norm <- unite(PRJEB6698_MethObj_100win_norm, destrand=FALSE)

## get percent methylation matrix
pm <- percMethylation(PRJEB6698_Meth)
pm1 <- percMethylation(PRJEB6698_MethObj_100win_norm)

## calculate standard deviation of CpGs
sds <- matrixStats::rowSds(pm)
sds1 <- matrixStats::rowSds(pm1)

## Visualize the distribution of the per-CpG standard deviation
## to determine a suitable cutoff
hist(sds, breaks = 100)
hist(sds1, breaks = 100)

## keep only CpG with standard deviations larger than 2%, ## This introduces NAs into matrix to be removed by na.omit
PRJEB6698_MethObj_100win_norm <- PRJEB6698_MethObj_100win_norm[sds > 2]
PRJEB6698_MethObj_100win_norm <- na.omit(PRJEB6698_MethObj_100win_norm)

## This leaves us with this number of CpG 100-bp windows (23,241)
nrow(PRJEB6698_MethObj_100win_norm)

## cluster all samples using correlation distance and return a tree object for plclust
hc1 <- clusterSamples(PRJEB6698_MethObj_100win_norm, dist="correlation", method="ward.D", plot=FALSE)

## cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(PRJEB6698_MethObj_100win_norm, dist="correlation", method="ward.D", plot=TRUE)

## screeplot of principal component analysis
PCASamples(PRJEB6698_MethObj_100win_norm, screeplot=TRUE)

## principal component anlaysis of all samples
PCASamples(PRJEB6698_MethObj_100win_norm)

## Split each Meth object by sample 2-way comparisons to identify intergroup differences; also move control 10 to position 3
PRJEB6698_Meth_100win_ConvD0 = reorganize(PRJEB6698_MethObj_100win_norm,sample.ids=metaData$Treatment[c(1:2, 10, 3:5)],treatment=c(0,0,0, 1,1,1))
PRJEB6698_Meth_100win_ConvD21 = reorganize(PRJEB6698_MethObj_100win_norm,sample.ids=metaData$Treatment[c(1:2, 10, 6:9)],treatment=c(0,0,0, 1,1,1,1))
PRJEB6698_Meth_100win_D0vD21 = reorganize(PRJEB6698_MethObj_100win_norm,sample.ids=metaData$Treatment[c(3:5, 6:9)],treatment=c(0,0,0,1,1,1,1))

####Differential Methylation
# calculate differential methylation p-values and q-values
# Takes some time to calculate
# if more than two groups calculated as abs(max difference between two groups; rather than to control)
PRJEB6698_Diff_100bp=calculateDiffMeth(PRJEB6698_MethObj_100win_norm)
PRJEB6698_Diff_100bp_ConvD0=calculateDiffMeth(PRJEB6698_Meth_100win_ConvD0)
PRJEB6698_Diff_100bp_ConvD21=calculateDiffMeth(PRJEB6698_Meth_100win_ConvD21)
PRJEB6698_Diff_100bp_D0vD21=calculateDiffMeth(PRJEB6698_Meth_100win_D0vD21)

# check how data part of methylDiff object looks like
head(PRJEB6698_Diff_100bp)
head(PRJEB6698_Diff_100bp_ConvD0)
head(PRJEB6698_Diff_100bp_ConvD21)
head(PRJEB6698_Diff_100bp_D0vD21)

# get differentially methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_100bp <- getMethylDiff(PRJEB6698_Diff_100bp, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_100bp_ConvD0 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD0, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_100bp_ConvD21 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD21, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_100bp_D0vD21 <- getMethylDiff(PRJEB6698_Diff_100bp_D0vD21, difference=25, qvalue=0.01)

# get differentially hypo methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hypo_100bp <- getMethylDiff(PRJEB6698_Diff_100bp, difference=25, qvalue=0.01, type="hypo")
PRJEB6698_Diff_25p_hypo_100bp_ConvD0 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD0, difference=25, qvalue=0.01, type="hypo")
PRJEB6698_Diff_25p_hypo_100bp_ConvD21 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD21, difference=25, qvalue=0.01, type="hypo")
PRJEB6698_Diff_25p_hypo_100bp_D0vD21 <- getMethylDiff(PRJEB6698_Diff_100bp_D0vD21, difference=25, qvalue=0.01, type="hypo")

# get differentially hyper methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hyper_100bp <- getMethylDiff(PRJEB6698_Diff_100bp, difference=25, qvalue=0.01, type="hyper")
PRJEB6698_Diff_25p_hyper_100bp_ConvD0 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD0, difference=25, qvalue=0.01, type="hyper")
PRJEB6698_Diff_25p_hyper_100bp_ConvD21 <- getMethylDiff(PRJEB6698_Diff_100bp_ConvD21, difference=25, qvalue=0.01, type="hyper")
PRJEB6698_Diff_25p_hyper_100bp_D0vD21 <- getMethylDiff(PRJEB6698_Diff_100bp_D0vD21, difference=25, qvalue=0.01, type="hyper")

# read-in transcript locations to be used in annotation
# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation
# files. Download appropriate annotation files from UCSC (or other sources) in BED format
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))

# annotate differentially methylated Cs with promoter/exon/intron using annotation data
annotate.WithGenicParts(myDiff25p,gene.obj)