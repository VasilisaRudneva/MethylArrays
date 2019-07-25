rm(list=ls())
setwd(dir="/Users/rudneva/Documents/myproject/Analysis")

suppressMessages(library("minfi")) 
suppressMessages(library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")) 
suppressMessages(library("IlluminaHumanMethylationEPICmanifest")) 
suppressMessages(library("minfiData"))

####################################################################################################################################
info=read.table("../Subgroup_Annotations.txt", header = T, sep = "\t", stringsAsFactors = F, na.strings = "N/A")
dim(info); head(info)

targets=read.table("../Data/path_to_idats.all.txt"); dim(targets); colnames(targets)=c("Sample_Plate", "folder","Basename")
targets=cbind(targets,targets$Basename,rep("ACNS0331", dim(targets)[1])); colnames(targets)=c("Sample_Plate", "folder","Basename","Sample_Name", "Trial"); targets$Basename=as.character(targets$Basename)
head(targets)
dim(targets)

targets.tmp=targets[targets$Basename %in% info[info$Manual_prediction %in% c("MB", "NA"),"idat"],]
head(targets.tmp)
dim(targets.tmp)
targets=targets.tmp
targets$folder=as.character(targets$folder)
head(targets)
dim(targets)

rm(targets.tmp)

folder1=unique(targets$folder)[1]
targets.tmp=targets[targets$folder==folder1,]
rgSet1= read.metharray.exp(base=paste0("/Volumes/Illumina/", targets.tmp$Sample_Plate, "/", folder1), 
                           targets =targets.tmp, force=T)
rgSet1$ArrayTypes=rep("IlluminaHumanMethylationEPIC", dim(rgSet1)[2])

folder2=unique(targets$folder)[2]
targets.tmp=targets[targets$folder==folder2,]
rgSet2= read.metharray.exp(base=paste0("/Volumes/Illumina/", targets.tmp$Sample_Plate, "/", folder2), 
                           targets =targets.tmp, force=T)
rgSet2$ArrayTypes=rep("IlluminaHumanMethylationEPIC", dim(rgSet2)[2])

rgSetold=combineArrays(rgSet1, rgSet2)
rgSetold

rgSet1
rgSet2

rm(rgSet1)
rm(rgSet2)

dim(rgSetold)

len=length(unique(targets$folder))
for (ind in 3:len){
  #print(ind)
  folder=unique(targets$folder)[ind]
  #print(folder)
  targets.tmp=targets[targets$folder==folder,]
  this_rgset = read.metharray.exp(base=paste0("/Volumes/Illumina/", targets.tmp$Sample_Plate, "/", folder), 
                                  targets =targets.tmp, force=T)
  this_rgset$ArrayTypes=rep("IlluminaHumanMethylationEPIC", dim(this_rgset)[2])
  
  rgSetnew=combineArrays(rgSetold, this_rgset)
  print(paste(ind, folder, dim(rgSetnew), sep=" "))
  rm(this_rgset)
  rm(rgSetold)
  rgSetold=rgSetnew
  rm(rgSetnew)
}

rgSet=rgSetold
rm(rgSetold)

save(rgSet, file="rgSet.RData")
