## MethylKit (re)-Analysis of DNMT1 Tet-OFF (Mcgraw et al. 2015, NAR)
## Erik Koppes July 12 2023

## Dependencies
library(readr)
library(methylKit)

## make a file list

MAIN_DIR <- "/ix1/mmann/KoppesEA/PRJEB6698"
BIS_DIR <- paste0(MAIN_DIR, "/Bismark")
METH_DIR <- paste0(MAIN_DIR, "/methylKit") ##wd

##actually make a meta data table with PRJEB6698 Sample_names and Metadata
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


##unsure if better to import cov or 
# _1_val_1_bismark_bt2_pe.bismark.cov.gz #coverage files will be unstranded
# _1_val_1_bismark_bt2_pe.CpG_report.txt.gz

## Generate PRJEB6698 methRead object
PRJEB6698_MethObj <-
        methRead(
                   location = PRJEB6698_Covfile_list_loc, ##list with directory and name of each file
                   sample.id = as.list(metaData$Treatment), # use treatment ID from metadata
                   assembly = "GRCm39",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = TRUE,
                   skip = 0,
                   sep = "\t",
                   context = "CpG",
                   resolution = "base",
                   treatment = c(1,1,0,0,0,2,2,2,2,1), #1 is control no-dox, 0 is d0 , 2 is d21 
                   dbdir = getwd(),
                   mincov = 10
                 )

# Get a histogram of the methylation percentage per sample
# Here for sample 1
getMethylationStats(PRJEB6698_MethObj[[1]], plot=TRUE, both.strands=FALSE)

# Get a histogram of the read coverage per sample
getCoverageStats(PRJEB6698_MethObj[[5]], plot=TRUE, both.strands=FALSE)
# Get percentile data by setting plot=FALSE
getCoverageStats(PRJEB6698_MethObj[[2]], plot=FALSE, both.strands=FALSE)

## Filtering
PRJEB6698_MethObj_fil <- filterByCoverage(PRJEB6698_MethObj,
                               lo.count=10, ##already done with input
                               lo.perc=NULL,
                               hi.count=NULL,
                               hi.perc=99.9) ##upperlimit cutoff implemented
## Normalization
PRJEB6698_MethObj_norm <- normalizeCoverage(PRJEB6698_MethObj_fil, method = "median")

## Merge Data
PRJEB6698_Meth <- unite(PRJEB6698_MethObj_norm, destrand=FALSE)

## get percent methylation matrix
pm <- percMethylation(PRJEB6698_Meth)

## calculate standard deviation of CpGs
sds <- matrixStats::rowSds(pm)

## Visualize the distribution of the per-CpG standard deviation
## to determine a suitable cutoff
hist(sds, breaks = 100)

## keep only CpG with standard deviations larger than 2%, ##WHY??
PRJEB6698_Meth <- PRJEB6698_Meth[sds > 2]

## This leaves us with this number of CpG sites (509,514)
nrow(PRJEB6698_Meth)

## cluster all samples using correlation distance and return a tree object for plclust
hc <- clusterSamples(PRJEB6698_Meth, dist="correlation", method="ward", plot=FALSE)

## cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(PRJEB6698_Meth, dist="correlation", method="ward", plot=TRUE)

## screeplot of principal component analysis.
PCASamples(PRJEB6698_Meth, screeplot=TRUE)

## principal component anlaysis of all samples.
PCASamples(PRJEB6698_Meth)

####Differential Methylation
# calculate differential methylation p-values and q-values
PRJEB6698_Diff=calculateDiffMeth(PRJEB6698_Meth)
## more than 2-groups detected but there are actually 3 can i use 0,1,2 for differential??

# check how data part of methylDiff object looks like
head(PRJEB6698_Diff)

# get differentially methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01)

# get differentially hypo methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hypo <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01, type="hypo") 

# get differentially hyper methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hyper <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01, type="hyper")
