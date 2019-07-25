rm(list=ls())
today=paste(strsplit(date(), " ")[[1]][c(2,3,5)], collapse = "-")
today
print("Loading dependencies")
library(minfi)

print("Loading dataset and info")
load("rgSet.RData")
info=read.table("Subgroup_Annotations.txt", header = T, sep = "\t", stringsAsFactors = F, na.strings = "N/A")
info=info[,c("idat", "Trial", "Batch", "PatientUSI", "batch", "sample", "prep", "array", "prediction")]; colnames(info)[9]="MNP_Subgroup"
head(info); dim(info)

for (i in 1:6){
  print(i)
  sel.samples=sampleNames(rgSet)[((i-1)*100+1):(i*100)]
  rgSet1=rgSet[,sel.samples]
  mset=preprocessRaw(rgSet1); rm(rgSet1)
  qc <- getQC(mset)
  head(qc)
  addQC(mset, qc=qc)
  
  RSet <- ratioConvert(mset, what = "both", keepCN = TRUE); rm(mset)
  GRset <- mapToGenome(RSet)
  
  sex_prediction=data.frame(cbind(colnames(GRset), getSex(GRset, cutoff = -2)$predictedSex), stringsAsFactors = F) ; rm(mSetRaw)
  colnames(sex_prediction)=c("idat", "predicted_sex")
  out.tmp=merge(sex_prediction, info, by="idat")
  head(out.tmp); dim(out.tmp)
  known=read.table("SexandAgePtLists.forR.txt", header = T, sep = "\t", stringsAsFactors = F)
  colnames(out.tmp)[5]="USI_PatientID"
  results=merge(known, out.tmp, by="USI_PatientID")
  write.table(results, file=paste0("Sex_Predictions.",i, ".txt"), quote = F, sep = "\t", col.names = T, row.names = F)
}
i=7
sel.samples=sampleNames(rgSet)[601:dim(rgSet)[2]]
rgSet1=rgSet[,sel.samples]
mset=preprocessRaw(rgSet1); rm(rgSet1)
RSet <- ratioConvert(mset, what = "both", keepCN = TRUE); rm(mset)
GRset <- mapToGenome(RSet)
sex_prediction=data.frame(cbind(colnames(GRset), getSex(GRset, cutoff = -2)$predictedSex), stringsAsFactors = F) ; rm(mSetRaw)
colnames(sex_prediction)=c("idat", "predicted_sex")
out.tmp=merge(sex_prediction, info, by="idat")
head(out.tmp); dim(out.tmp)
known=read.table("SexandAgePtLists.forR.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(out.tmp)[5]="USI_PatientID"
results=merge(known, out.tmp, by="USI_PatientID")
write.table(results, file=paste0("ex_Predictions.",i, ".txt"), quote = F, sep = "\t", col.names = T, row.names = F)


### making plots
submitted_dna=read.xls("Summary_DNA_extraction.xlsx", sheet = 2)
submitted_dna=submitted_dna[c("PatientUSI", "Total.DNA..ng.")]
colnames(submitted_dna)[1]="USI_PatientID"
result1=read.xls("Sex_Predictions.xlsx", sheet = 1)
results=merge(result1, submitted_dna, by="USI_PatientID"); head(results); dim(results)

concordant=results[results$Gender==results$predicted_sex,]; dim(concordant)
different=results[results$Gender!=results$predicted_sex,]; dim(different) 

df1=cbind(concordant[c("Gender", "predicted_sex", "Total.DNA..ng.")], rep("Conc", dim(concordant)[1]))
colnames(df1)[4]="Status"
df2=cbind(different[c("Gender", "predicted_sex", "Total.DNA..ng.")], rep("Diff", dim(different)[1]))
colnames(df2)[4]="Status"

df=data.frame(rbind(df1, df2))
library("reshape")
library("ggplot2")

ddf=melt(df)


pdf(paste0("Sample_swap/", today,".Boxplots_TotalDNA.pdf"), paper="a4", onefile=TRUE, useDingbats=FALSE)
ggplot(ddf, aes(x=Status, y=value, color=Gender)) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Total DNA (ng)")+ 
  ggtitle(paste0("Colored by Gender; t-test p-value = ", round(t.test(ddf[ddf$Status=="Conc",]$value, ddf[ddf$Status=="Diff",]$value, alternative = "greater")$p.value, digits = 3)))
ggplot(ddf, aes(x=Status, y=value, color=predicted_sex)) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Total DNA (ng)")+ 
  ggtitle(paste0("Colored by Predicted Sex; t-test p-value = ", round(t.test(ddf[ddf$Status=="Conc",]$value, ddf[ddf$Status=="Diff",]$value, alternative = "greater")$p.value, digits = 3)))
dev.off()


  
### Add QC
rm(list=ls())
today=paste(strsplit(date(), " ")[[1]][c(2,3,5)], collapse = "-")
today
library(minfi)
library("reshape")
library("ggplot2")
require("gdata")

load("rgSet.RData")
info=read.table("Subgroup_Annotations.txt", header = T, sep = "\t", stringsAsFactors = F, na.strings = "N/A")
info=info[,c("idat", "Trial", "Batch", "PatientUSI", "batch", "sample", "prep", "array", "prediction")]; colnames(info)[9]="MNP_Subgroup"
head(info); dim(info)
badSampleCutoff=10.5
mset=preprocessRaw(rgSet); 
qc <- getQC(mset)
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < badSampleCutoff))
badSamples=rownames(qc)[whichBad]
plotQC(qc)

submitted_dna=read.xls("Summary_DNA_extraction.xlsx", sheet = 2)
submitted_dna=submitted_dna[c("PatientUSI", "Total.DNA..ng.")]
colnames(submitted_dna)[1]="USI_PatientID"
result1=read.xls("Sex_Predictions.xlsx", sheet = 1)
results=merge(result1, submitted_dna, by="USI_PatientID"); 
results$QC=rep("Pass", dim(results)[1])
results[results$idat %in% badSamples,"QC"]="Fail"
head(results); dim(results)
concordant=results[results$Gender==results$predicted_sex,]; dim(concordant)
different=results[results$Gender!=results$predicted_sex,]; dim(different) 
df1=cbind(concordant[c("Gender", "predicted_sex", "Total.DNA..ng.", "QC")], rep("Conc", dim(concordant)[1]))
colnames(df1)[5]="Status"
df2=cbind(different[c("Gender", "predicted_sex", "Total.DNA..ng.", "QC")], rep("Diff", dim(different)[1]))
colnames(df2)[5]="Status"
df=data.frame(rbind(df1, df2))

ddf=melt(df)

pdf(paste0("Sample_swap/", today,".Boxplots_TotalDNA_QC.pdf"), paper="a4", onefile=TRUE, useDingbats=FALSE)
ggplot(ddf, aes(x=Status, y=value, color=QC)) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Total DNA (ng)")+ 
  ggtitle(paste0("Colored by QC results; Pass (n=577) vs Fail (n=73) total DNA amount t-test p-value = ", round(t.test(ddf[ddf$QC=="Pass",]$value, ddf[ddf$QC=="Fail",]$value, alternative = "greater")$p.value, digits = 3)))
dev.off()


### Patients with suspected sample swap
pt.ids=c("PARLRA","PASUKA","PARWRL","PASMHP","PAVHWJ","PAWRCG")
pt.pairs=c("pair1","pair1","pair2","pair2","pair3","pair3")
pt.mix=data.frame(cbind(pt.ids, pt.pairs), stringsAsFactors = F); colnames(pt.mix)=c("USI_PatientID", "Pair")

results[results$USI_PatientID %in% pt.mix$USI_PatientID,]

head(results); dim(results)
concordant=results[results$Gender==results$predicted_sex,]; dim(concordant)
different=results[results$Gender!=results$predicted_sex,]; dim(different) 
df1=cbind(concordant[c("Gender", "predicted_sex", "Total.DNA..ng.", "QC", "Mix")], rep("Conc", dim(concordant)[1]))
colnames(df1)[6]="Status"
df2=cbind(different[c("Gender", "predicted_sex", "Total.DNA..ng.", "QC", "Mix")], rep("Diff", dim(different)[1]))
colnames(df2)[6]="Status"
df=data.frame(rbind(df1, df2))

ddf=melt(df)

pdf(paste0("Sample_swap/", today,".Boxplots_TotalDNA_QC.pdf"), paper="a4", onefile=TRUE, useDingbats=FALSE)
p=ggplot(ddf, aes(x=Status, y=value, color=QC)) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Total DNA (ng)")+ 
  ggtitle(paste0("Colored by QC results; Pass (n=577) vs Fail (n=73) total DNA amount t-test p-value = ", round(t.test(ddf[ddf$QC=="Pass",]$value, ddf[ddf$QC=="Fail",]$value, alternative = "greater")$p.value, digits = 3)))
p + geom_text(aes(label = rownames(df)),
              size = 3.5)
dev.off()
