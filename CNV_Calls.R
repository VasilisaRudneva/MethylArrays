setwd(dir="/myproject/Analysis/CNV_Calling/")

amp_formula="mean(res[,5])+sd(res[,5])"
del_formula="mean(res[,5])-sd(res[,5])"

today=paste(strsplit(date(), " ")[[1]][c(2,4,6)], collapse = "-")

print(today)
print(paste("Amp: ", amp_formula, " Del: ", del_formula, sep=""))
options(warn=-1)

suppressMessages(library(minfi)) 
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19)) 
suppressMessages(library(IlluminaHumanMethylation450kmanifest)) 
suppressMessages(library(minfiData))
suppressMessages(library(stringr))
suppressMessages(library("conumee"))

chrOrder=c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p",
           "10q", "11p", "11q","12p", "12q", "13p", "13q", "14p", "14q","15p", "15q", "16p", "16q", "17p",
           "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q", "Xp", "Xq", "Yp", "Yq")
chrOrder=paste("chr", chrOrder, sep="")

data(RGcontrolSetEx)
data(exclude_regions)
data(detail_regions)

targets=read.table("/myproject/Targets.txt", header=T, stringsAsFactors = F); 
dim(targets); colnames(targets)[5]="Basename"; colnames(targets)[4]="Sample_Name"
targets$Basename=as.character(targets$Basename)
# Not COG-0332 ==> exclude:
exclude_COG_0332=c("201985320234_R08C01")
targets=targets[!targets$Basename %in% exclude_COG_0332, ]
targetsCOG=targets
dim(targets)
load("/myproject/DATA/RData/gSetCOG.RData")
rgSet=rgSetCOG

RGcontrolSetEx=read.metharray.exp("/data/methylation/Normal_Refs/")
MsetControls <- preprocessRaw(RGcontrolSetEx)
anno <- CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions)
#query.data=CNV.load(preprocessRaw(rgSet))
mSetRaw <- preprocessRaw(rgSet)

probes=featureNames(MsetControls)[featureNames(MsetControls) %in% featureNames(mSetRaw)]

mSetRaw=mSetRaw[probes,]
query.data=CNV.load(mSetRaw)

all_samples=colnames(getM(mSetRaw))
MsetControls <- MsetControls[probes,]
controls.data <- CNV.load(MsetControls)

load("cytoBand.RData")

# Temporary from the internet
Mset <- mapToGenome(mSetRaw)
anno@probes <- subsetByOverlaps(anno@probes, granges(Mset))

controls.data
query.data
anno

output_table=targetsCOG
output_table$CNV=rep(NA, dim(targetsCOG)[1])
output_table$Drivers=rep(NA, dim(targetsCOG)[1])
colnames(output_table)[4]="PID"
head(output_table)

Drivers=c("MYCN","GLI2","MYC","PTCH1","PTEN","SUFU","TP53")
Chromosomes=paste("chr", 1:22, sep="")

#for (a2.filename in metadataInfants[metadataInfants$PID %in% filenames, "METH_450K"])
index=1
for (a2.filename in names(query.data))
{
  print(paste0(index, ": ", a2.filename))
  
  filename=a2.filename
  
  x <- CNV.fit(query.data[a2.filename], controls.data, anno)

  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
  #
  
  #CNV.detailplot_wrap(x)
  out=""
  out.detail=""
  res=CNV.write(x, what="bins")
  
  ext_low=res[which(res[,5]< -0.4),]
  ext_high=res[which(res[,5]> 0.4),]
    
  qc_out=vector()
    
  e=ext_low
    #print("Extreme low...")
  for (c in paste("chr", seq(1:22), sep="")){
      #print(paste(c, " : ", mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]), sep=""))
      qc_out=c(qc_out, mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]))
      
  }
  e=ext_high
    #print("Extreme high...")
  for (c in paste("chr", seq(1:22), sep="")){
      #print(paste(c, " : ", mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]), sep=""))
      qc_out=c(qc_out, mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]))
  }
    
    if (length(qc_out[qc_out> 0]) <= 15){
      
      amp_cutoff=eval(parse(text=amp_formula))
      del_cutoff=eval(parse(text=del_formula))
      
      for (i in 1:22)
      {
        
        chr=paste("chr", i, sep="")
        p.start=cyto.mappings[i,1];p.end=cyto.mappings[i,2]
        q.start=cyto.mappings[i,3];q.end=cyto.mappings[i,4]
        
        this_res=res[res$Chromosome==chr,]
        med.p=median(as.numeric(this_res[(this_res$Start>=p.start)&(this_res$End<=p.end),][,5]))
        med.q=median(as.numeric(this_res[(this_res$Start>=q.start)&(this_res$End<q.end),][,5]))
        
        if(is.na(med.p)){out=paste(out, chr, "p=NA; ", sep="")} else {
          if (med.p > amp_cutoff) {out=paste(out, chr, "p=AMP; ", sep="")} else {
            if (med.p < del_cutoff) {out=paste(out, chr, "p=DEL; ", sep="")} else {out=paste(out, chr, "p=0; ", sep="")}
          }
        }
        if(is.na(med.q)){out=paste(out, chr, "q=NA; ", sep="")} else {
          if (med.q > amp_cutoff) {out=paste(out, chr, "q=AMP; ", sep="")} else {
            if (med.q < del_cutoff) {out=paste(out, chr, "q=DEL; ", sep="")} else {out=paste(out, chr, "q=0; ", sep="")}
          }
        }
      }
      
      Meth_CNV=unlist(lapply(unlist(strsplit(out, split = "; ")), function(x) strsplit(x,"=")[[1]][2])); Meth_CNV[Meth_CNV=="NA"]=NA
      names(Meth_CNV)=chrOrder[1:44]; Meth_CNV=Meth_CNV[which(!is.na(Meth_CNV))]; Meth_CNV=Meth_CNV[which(Meth_CNV!=0)]
      if (length(Meth_CNV)==0){out2_a="No"; out2_d="No"} else {
        Meth_CNV_A=Meth_CNV[which(Meth_CNV=="AMP")]; Meth_CNV_D=Meth_CNV[which(Meth_CNV=="DEL")]
        out2_a=""; for (c in sort(names(Meth_CNV_A))){out2_a=paste(out2_a, c, ";", sep="")}
        out2_d=""; for (c in sort(names(Meth_CNV_D))){out2_d=paste(out2_d, c, ";", sep="")}
      }
      
      ### Copy number for the regions of interest
      res.detail=CNV.write(x, what="detail")
      for (g in res.detail$name) {
        if (g %in% c("MYC", "MYCN", "GLI2")){
          if (res.detail[res.detail$name==g,]$value >= 0.4 ){out.detail=paste(out.detail, g, "=AMP; ", sep="")}
          else {out.detail=paste(out.detail, g, "=0; ", sep="")}
        } else {
          if (res.detail[res.detail$name==g,]$value <= -0.4 ){out.detail=paste(out.detail, g, "=DEL; ", sep="")}
          else {out.detail=paste(out.detail, g, "=0; ", sep="")}
        }
      }
      
      Driver_CNV=unlist(lapply(unlist(strsplit(out.detail, split = "; ")), function(x) strsplit(x,"=")[[1]][2])); Driver_CNV[Driver_CNV=="NA"]=NA
      names(Driver_CNV)=Drivers; Driver_CNV=Driver_CNV[which(!is.na(Driver_CNV))]; Driver_CNV=Driver_CNV[which(Driver_CNV!=0)]
      if (length(Driver_CNV)==0){out1_a="No"; out1_d="No"} else {
        Driver_CNV_A=Driver_CNV[which(Driver_CNV=="AMP")]; Driver_CNV_D=Driver_CNV[which(Driver_CNV=="DEL")]
        out1_a=""; for (c in sort(names(Driver_CNV_A))){out1_a=paste(out1_a, c, ";", sep="")}
        out1_d=""; for (c in sort(names(Driver_CNV_D))){out1_d=paste(out1_d, c, ";", sep="")}
      }
      
      pdf(paste(index, "_", filename,  ".", "CNV.", today, ".pdf", sep=""), width = 17.5, height = 10, onefile=TRUE, useDingbats=FALSE)
      CNV.genomeplot(x)
      legend("topleft",legend = paste("CNV:\nAMP: ", out2_a, "\nDEL: ", out2_d,"\nDrivers:\nAMP: ", out1_a, " DEL: ", out1_d, sep=""), cex=0.75, x.intersp=1, y.intersp=7)
      dev.off()
      
      output_table[output_table$PID==filename,]$CNV=out
      output_table[output_table$PID==filename,]$Drivers=out.detail
      
    } else {
      
    pdf(paste("Excluded/", index, "_", filename, ".CNV.", today,".pdf", sep=""), width = 17.5, height = 10, onefile=TRUE, useDingbats=FALSE)
    CNV.genomeplot(x)
    dev.off()
    }
  
  seg_info <- CNV.write(x, what="segments")
  cnv_segs <- data.frame(Chromosome=seg_info$chrom, Start=seg_info$loc.start,End=seg_info$loc.end, Value=seg_info$seg.mean)
  write.table(cnv_segs, paste0(index, "_", "Segs_", filename,".CNV_segments.",today,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  index=index+1
}


write.table(x = output_table, file = paste("CNV_Calls.", today, ".txt", sep=""), sep="\t", col.names = T, row.names = T, quote=F)
save(output_table, file="/Users/vrudneva/Documents/COG_0331-0332/Analysis/CNV_calls.RData")

tab=output_table

chrOrder=c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p",
           "10q", "11p", "11q","12p", "12q", "13p", "13q", "14p", "14q","15p", "15q", "16p", "16q", "17p",
           "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q")
chrOrder=paste("chr", chrOrder, sep="")
Drivers=c("MYCN","GLI2","MYC","PTCH1","PTEN","SUFU","TP53")

tab$CNV=gsub(x = tab$CNV, pattern = " ", replacement = "")
CNV2=lapply(tab[!is.na(tab$CNV),]$CNV,function(x) data.frame(names=unlist(lapply(unlist(strsplit(as.character(x), split = ";")), function(x) strsplit(x,"=")[[1]][1])), values=unlist(lapply(unlist(strsplit(as.character(x), split = ";")), function(x) strsplit(x,"=")[[1]][2]))))
names(CNV2)=tab[!is.na(tab$CNV),]$PID

m <- cbind(as.character(tab$PID), matrix(NA, ncol = length(chrOrder)+length(Drivers)+2, nrow = dim(tab)[1])); m <- data.frame(m, stringsAsFactors = F)
colnames(m)=c("PID", Drivers, chrOrder, "Drivers", "CNV")

for (s in names(CNV2)){
  tmp=CNV2[s][[1]]; rownames(tmp)=tmp[,1]
  m[m$PID==s, "Drivers"]=as.character(tab[tab$PID==s,"Drivers"])
  m[m$PID==s, "CNV"]=as.character(tab[tab$PID==s,"CNV"])
  for (c in chrOrder){m[m$PID==s,c]=as.character(tmp[c,"values"])}
}

Dr=lapply(tab[!is.na(tab$Drivers),]$Drivers,function(x) data.frame(names=unlist(lapply(unlist(strsplit(as.character(x), split = ";")), function(x) strsplit(x,"=")[[1]][1])), values=unlist(lapply(unlist(strsplit(as.character(x), split = ";")), function(x) strsplit(x,"=")[[1]][2]))))
names(Dr)=tab[!is.na(tab$Drivers),]$PID

for (s in names(Dr)){
  print(s)
  tmp=Dr[s][[1]]; rownames(tmp)=tmp[,1]
  rownames(tmp)=gsub(" ", "", rownames(tmp))
  for (d in Drivers){m[m$PID==s,d]=as.character(tmp[d,"values"])}
}
head(m)

write.table(x = m, file = paste("CNV_Calls.", today, ".parsed.txt", sep=""), sep="\t", col.names = T, row.names = T, quote=F)
