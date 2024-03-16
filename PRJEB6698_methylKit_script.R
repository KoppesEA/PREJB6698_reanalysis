## MethylKit (re)-Analysis of DNMT1 Tet-OFF (Mcgraw et al. 2015, NAR)
## Erik Koppes, PhD
## Magee Womens Research Institute
## March 2024

# Install methylKit as needed
BiocManager::install("methylKit")
BiocManager::install("genomation")

## Dependencies
library(readr)
library(methylKit) 
library(genomation)
packageVersion("methylKit") #v1.28.0 on M1 Macbook Pro
packageVersion("genomation") #v1.34.0 on M1 Macbook Pro

## make a file list

MAIN_DIR <- "./"
BIS_DIR <- paste0(MAIN_DIR, "Bismark")
METH_DIR <- paste0(MAIN_DIR, "methylKit") ##wd

##actually make a meta data table with PRJEB6698 Sample_names and Metadata
PRJEB6698_list <- c("ERR560526", "ERR560527", "ERR560528", "ERR560529", "ERR560530",
                    "ERR560531", "ERR560532", "ERR560533", "ERR560534", "ERR560538")
DIR_list <- as.list(paste0(BIS_DIR, "/", PRJEB6698_list, "_methextract/"))
PRJEB6698_Covfile_list <- as.list(paste0(PRJEB6698_list, "_1_val_1_bismark_bt2_pe.bismark.cov.gz"))
PRJEB6698_CpGreport_list <- as.list(paste0(PRJEB6698_list, "_1_val_1_bismark_bt2_pe.CpG_report.txt.gz"))
PRJEB6698_Covfile_list_loc <- as.list(paste0(BIS_DIR, "/", PRJEB6698_list, "_methextract/",PRJEB6698_list, "_1_val_1_bismark_bt2_pe.bismark.cov.gz"))

##Import metadata TSV (see script to obtain from .xml)
metaData <- read_tsv(paste0(MAIN_DIR, "PRJEB6698_metadata.tsv"),
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
                   treatment = c(0,0,1,1,1,2,2,2,2,0), #0 is control no-dox, 1 is d0 , 2 is d21 
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
hc <- clusterSamples(PRJEB6698_Meth, dist="correlation", method="ward.D", plot=FALSE)

## cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(PRJEB6698_Meth, dist="correlation", method="ward.D", plot=TRUE)

## screeplot of principal component analysis.
PCASamples(PRJEB6698_Meth, screeplot=TRUE)

## principal component anlaysis of all samples.
PCASamples(PRJEB6698_Meth)

## Split each Meth object by sample 2-way comparisons to identify intergroup differences; also move control 10 to position 3
PRJEB6698_Meth_ConvD0 = reorganize(PRJEB6698_Meth,sample.ids=metaData$Treatment[c(1:2, 10, 3:5)],treatment=c(0,0,0, 1,1,1))
PRJEB6698_Meth_ConvD21 = reorganize(PRJEB6698_Meth,sample.ids=metaData$Treatment[c(1:2, 10, 6:9)],treatment=c(0,0,0, 1,1,1,1))
PRJEB6698_Meth_D0vD21 = reorganize(PRJEB6698_Meth,sample.ids=metaData$Treatment[c(3:5, 6:9)],treatment=c(0,0,0,1,1,1,1))

####Differential Methylation
# calculate differential methylation p-values and q-values
# if more than two groups calculated as abs(max difference between two groups; rather than to control)
# may take a few minutes to run
PRJEB6698_Diff=calculateDiffMeth(PRJEB6698_Meth)
PRJEB6698_Diff_ConvD0=calculateDiffMeth(PRJEB6698_Meth_ConvD0)
PRJEB6698_Diff_ConvD21=calculateDiffMeth(PRJEB6698_Meth_ConvD21)
PRJEB6698_Diff_D0vD21=calculateDiffMeth(PRJEB6698_Meth_D0vD21)

## more than 2-groups detected and there are actually 3; differential tested to average mean me

# check how data part of methylDiff object looks like
head(PRJEB6698_Diff)
head(PRJEB6698_Diff_ConvD0)
head(PRJEB6698_Diff_ConvD21)
head(PRJEB6698_Diff_D0vD21)


# get differentially methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_ConvD0 <- getMethylDiff(PRJEB6698_Diff_ConvD0, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_ConvD21 <- getMethylDiff(PRJEB6698_Diff_ConvD21, difference=25, qvalue=0.01)
PRJEB6698_Diff_25p_D0vD21 <- getMethylDiff(PRJEB6698_Diff_D0vD21, difference=25, qvalue=0.01)

# get differentially hypo methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hypo <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01, type="hypo") #no non-missing
PRJEB6698_Diff_25p_hypo_ConvD0 <- getMethylDiff(PRJEB6698_Diff_ConvD0, difference=25, qvalue=0.01, type="hypo")
PRJEB6698_Diff_25p_hypo_ConvD21 <- getMethylDiff(PRJEB6698_Diff_ConvD21, difference=25, qvalue=0.01, type="hypo")
PRJEB6698_Diff_25p_hypo_D0vD21 <- getMethylDiff(PRJEB6698_Diff_D0vD21, difference=25, qvalue=0.01, type="hypo") #no-nonmissing
# get differentially hyper methylated regions with 25% difference and qvalue<0.01
PRJEB6698_Diff_25p_hyper <- getMethylDiff(PRJEB6698_Diff, difference=25, qvalue=0.01, type="hyper")
PRJEB6698_Diff_25p_hyper_ConvD0 <- getMethylDiff(PRJEB6698_Diff_ConvD0, difference=25, qvalue=0.01, type="hyper") #no non-missing
PRJEB6698_Diff_25p_hyper_ConvD21 <- getMethylDiff(PRJEB6698_Diff_ConvD21, difference=25, qvalue=0.01, type="hyper")
PRJEB6698_Diff_25p_hyper_D0vD21 <- getMethylDiff(PRJEB6698_Diff_D0vD21, difference=25, qvalue=0.01, type="hyper")

# read-in transcript locations to be used in annotation
# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation
# files. Download appropriate annotation files from UCSC (or other sources) in BED format
CGI_BED <- paste0(MAIN_DIR, "GRCm39_CGI_Ensembl.bed")
GRCm39_CGIgr <- readGeneric(CGI_BED, header = TRUE, sep = "\t", meta.cols = list(name=4))

# annotate differentially methylated Cs with promoter/exon/intron using annotation data
annotateWithFeature(as(PRJEB6698_Diff_25p,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hypo,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hyper,"GRanges"), GRCm39_CGIgr)

annotateWithFeature(as(PRJEB6698_Diff_25p_ConvD0,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hypo_ConvD0,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hyper_ConvD0,"GRanges"), GRCm39_CGIgr)

annotateWithFeature(as(PRJEB6698_Diff_25p_ConvD21,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hypo_ConvD21,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hyper_ConvD21,"GRanges"), GRCm39_CGIgr)

annotateWithFeature(as(PRJEB6698_Diff_25p_hyper_D0vD21,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hypo_D0vD21,"GRanges"), GRCm39_CGIgr)
annotateWithFeature(as(PRJEB6698_Diff_25p_hyper_D0vD21,"GRanges"), GRCm39_CGIgr)

# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

annotateWithGeneParts()
annotateWithFeatures()