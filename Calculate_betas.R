rm(list=ls())
setwd(dir="/Users/rudneva/myproject/Analysis/")
source("/Users/rudneva/scripts/MethyationDataPreprocessingFiltering.R")
today=paste(strsplit(date(), " ")[[1]][c(2,4,6)], collapse = "-")
today
library("minfi")
library("limma")
library("Rtsne")
library("weights")

load("rgSet.RData")
info=read.table("Subgroup_Annotations.txt", header = T, sep = "\t", stringsAsFactors = F, na.strings = "N/A")

material=info[match(sampleNames(rgSet), info$idat),c("idat", "prep")]; colnames(material)=c("id", "mat");
rownames(material)=seq(1:dim(material)[1])

bVals=PreprocessFilterMethylationData(rgSet,material)
betas=bVals

save(betas, file="betas.RData")
