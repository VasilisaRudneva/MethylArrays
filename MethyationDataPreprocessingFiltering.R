### Methylation data preprocessing and t-SNE
suppressMessages(library("minfi")) 
suppressMessages(library("limma")) 
suppressMessages(library("minfiData"))
suppressMessages(library("stringr"))
suppressMessages(library("Rtsne"))
suppressMessages(library("weights"))

## Preprocessing
PreprocessFilterMethylationData <- function(rgSet, material) {

mSetRaw <- preprocessRaw(rgSet)
all_samples=sampleNames(rgSet)

Meth=log2(minfi::getMeth(mSetRaw))
Unmeth=log2(minfi::getUnmeth(mSetRaw))

data <- data.frame(id=colnames(Meth)); data$vector1 <- material[match(data$id, material$id),]
dim(data); head(data)

Meth=Meth[,as.vector(na.omit(data$vector1$id))]
Unmeth=Unmeth[,as.vector(na.omit(data$vector1$id))]

if (length(unique(data$vector1$mat))>1){
  
  Meth2=removeBatchEffect(Meth, batch = as.vector(data$vector1$mat))
  Unmeth2=removeBatchEffect(Unmeth, batch = as.vector(data$vector1$mat))
  
  Meth=2^Meth2
  Unmeth=2^Unmeth2
}

mSetRaw=MethylSet(Meth, Unmeth)

# Filtering
# load appropriate library
if (annotation(rgSet)[[1]] == "IlluminaHumanMethylationEPIC"){
  suppressMessages(library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19"))
  suppressMessages(library("IlluminaHumanMethylationEPICmanifest")) 
  annThisArrayType = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19); head(annThisArrayType)
} else {
  suppressMessages(library("IlluminaHumanMethylation450kanno.ilmn12.hg19")) 
  suppressMessages(library("IlluminaHumanMethylation450kmanifest")) 
  annThisArrayType = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19); head(annThisArrayType)
}

# remove probes on the sex chromosomes
keep <- !(featureNames(mSetRaw) %in% annThisArrayType$Name[annThisArrayType$chr %in% c("chrX","chrY")])
table(keep)
mSetRaw <- mSetRaw[keep,]; dim(mSetRaw)
# exclude cross reactive probes
xReactiveProbes <- read.csv(file="/Users/rudneva/scripts/MethylationProbes_to_filter/450k/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE) 
keep <- !(featureNames(mSetRaw) %in% xReactiveProbes$TargetID)
table(keep)
mSetRaw <- mSetRaw[keep,]
dim(mSetRaw)
# remove probes with any SNP
snps=getSnpInfo(rgSet)
remove=rownames(snps[!is.na(snps$SBE_rs) | !is.na(snps$CpG_rs),])
keep <- !(featureNames(mSetRaw) %in% remove)
table(keep)
mSetRaw <- mSetRaw[keep,]
dim(mSetRaw)
# remove probes that are not uniquely mapped to the hg19 genome
# BOWTIE2 multi-mapped
multi.map <- read.csv('/Users/rudneva/scripts/MethylationProbes_to_filter/450k/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', head = F, as.is = T)
multi.map.probes <- as.character(multi.map$V1)
keep <- !(featureNames(mSetRaw) %in% multi.map.probes)
table(keep)
mSetRaw <- mSetRaw[keep,]
dim(mSetRaw)

### Oct-2018 Additional probes for EPIC arrays. Reference: https://github.com/sirselim/illumina450k_filtering
if (annotation(rgSet)[[1]] == "IlluminaHumanMethylationEPIC"){
  suppressMessages(library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19"))
  suppressMessages(library("IlluminaHumanMethylationEPICmanifest")) 
  # probes from Pidsley 2016 (EPIC)
  epic.cross1 <- read.csv('/Users/rudneva/scripts/MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM1_ESM.csv', head = T)
  epic.variants1 <- read.csv('/Users/rudneva/scripts/MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM4_ESM.csv', head = T)
  epic.variants2 <- read.csv('/Users/rudneva/scripts/MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM5_ESM.csv', head = T)
  epic.variants3 <- read.csv('/Users/rudneva/scripts/MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM6_ESM.csv', head = T)
  # additional filter probes
  epic.add.probes <- c(as.character(epic.cross1$X), as.character(epic.variants1$PROBE), as.character(epic.variants2$PROBE), 
                     as.character(epic.variants3$PROBE))
  # final list of unique probes
  epic.add.probes <- unique(epic.add.probes)

  keep <- !(featureNames(mSetRaw) %in% epic.add.probes)
  table(keep)
  mSetRaw <- mSetRaw[keep,]
  dim(mSetRaw)
}

bVals <- getBeta(mSetRaw, offset=100); head(bVals[,1:5]); dim(bVals)
names=colnames(bVals)

return(bVals)
}

# Run standard t-SNE
RuntSNE <- function(bVals) {

betas=bVals
betas.sd <- apply(betas, 1, sd)
m.sd <- names(sort(betas.sd, decreasing = TRUE))

## cor
s <- "all"
cor.mXprobes <- list()

# cor weighted, set lowest used probe to 0
p = 0.25
p <- length(which(betas.sd > p))
cor.mXprobes[["all-sd-w-0.25probes"]] <- wtd.cors(betas[m.sd[1:p],], y=NULL, weight=betas.sd[m.sd[1:p]]-betas.sd[m.sd[p]])


## tsne seed (figure out which seed to use)
set.seed(20160501)
Y <- Rtsne(1-cor.mXprobes[["all-sd-w-0.25probes"]], pca = FALSE, verbose = FALSE,
             is_distance = TRUE, theta = 0, max_iter = 2000, dims = 2, perplexity = 10)$Y
rownames(Y) <- rownames(cor.mXprobes[["all-sd-w-0.25probes"]])

return(Y)
}

## plotting function for tsne
plot.tsne <- function(Y, xlab="TSNE 1", ylab="TSNE 2", pch=21, col="black", bg="lightgrey", las=2, centroid=TRUE, ...) {
  Y.range <- apply(Y, 2, range)
  Y.diff <- apply(Y.range, 2, diff)
  Y.center <- apply(Y.range, 2, mean)
  
  set.seed("123")
  if(length(bg) > 1) {
    Y.order <- order(is.element(bg, "lightgrey"), is.element(bg, "white"), decreasing = TRUE)
  } else {
    Y.order <- seq(nrow(Y))
  }
  if(length(bg) > 1) bg <- bg[Y.order]
  if(length(pch) > 1) pch <- pch[Y.order]
  if(length(col) > 1) col <- col[Y.order]
  
  par(mar=c(4, 4, 4, 4))
  plot(Y[Y.order, ], xlim=Y.center[1] + c(-0.5, 0.5)*max(Y.diff), ylim=Y.center[2] + c(-0.5, 0.5)*max(Y.diff),
       xlab=xlab, ylab=ylab, pch=pch, bg=bg, col=col, las=las, ...)
  
  if(length(unique(bg))>1 & centroid) {
    Y.centroid <- cbind(sapply(split(Y[, 1], bg), mean), sapply(split(Y[, 2], bg), mean))
    for(i in seq(nrow(Y))) {
      lines(c(Y.centroid[bg[i], 1], Y[i, 1]), c(Y.centroid[bg[i], 2], Y[i, 2]), col=bg[i])
    }
  }
}
