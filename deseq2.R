
library(DESeq2)
library(org.Mm.eg.db)
library(openxlsx)
library(UpSetR)
library(pheatmap)
library(IHW)
library(ggplot2)
library(ggrepel)
library(data.table)


source("/Users/gandrieux/bitbucket/work/tools/iwanthue.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

toNum <- function(x) return(as.numeric(levels(x))[x])

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Mm.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

entrez2ensembl <- function(entrez)
{
  esbl <- mget(as.character(entrez), org.Mm.egENSEMBL, ifnotfound=NA)
  esbl <- lapply(esbl, function(i) return(i[1]))
  return(unlist(esbl))
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Mm.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

entrez2name <- function(entrez)
{
  gn <- mget(as.character(entrez), org.Mm.egGENENAME, ifnotfound=NA)
  gn <- unlist(lapply(gn, function(i) return(i[1])))
  return(gn)
}


getNbDEG <- function(deseq, pv, fc, adjusted = FALSE)
{
  deseqFC <- deseq$log2FoldChange
  if(class(deseqFC)=="factor") deseqFC <- as.numeric(levels(deseqFC))[deseqFC]
  ifelse(adjusted, deseqPV <- deseq$padj, deseqPV <- deseq$pvalue)
  if(class(deseqPV)=="factor") deseqPV <- as.numeric(levels(deseqPV))[deseqPV]
  
  return(sum(abs(deseqFC)>= fc & deseqPV <= pv, na.rm = TRUE))
}

getDEGMatrix <- function(ddsObj, contrast, rld = NULL)
{
  dds.ann.gene <- as.data.frame(rowData(ddsObj))
  dds.ann.samples <- as.data.frame(colData(ddsObj))
  dds.ann.samples <- dds.ann.samples[order(dds.ann.samples$Group), ]# re-order rows according to group
  
  #dds.counts <- counts(ddsObj, normalized=TRUE)
  #dds.counts <- dds.counts[, match(dds.ann.samples$Sample, colnames(dds.counts))]
  #colnames(dds.counts) <- paste0(colnames(dds.counts), ".count")
  
  resLFC <- lfcShrink(ddsObj, contrast= contrast, cooksCutoff=FALSE, independentFiltering=FALSE)
  
  res.final <- as.data.frame(resLFC)
  res.final <- cbind(dds.ann.gene[, c("entrez", "symbol", "ensembl", "gene.name")], res.final)# add gene annotation
  #res.final <- cbind(res.final, dds.counts)
  
  if(!is.null(rld))
  {
    dds.rld <- rld
    dds.rld <- dds.rld[, match(dds.ann.samples$SAMPLE, colnames(dds.rld))]
    colnames(dds.rld) <- paste0(colnames(dds.rld), ".rLog")
    res.final <- cbind(res.final, dds.rld)
  }
  
  res.final <- res.final[order(res.final$pvalue),]
  
  return(res.final)
}


ggVolcano <- function(deseq, pvalue = 0, genes = NULL, outFile){
  ggmat <- data.frame(X = deseq$log2FoldChange, Y = -log10(deseq$padj), GENE = deseq$symbol)
  ggmat$isSignif <- "no"
  ggmat$isSignif[deseq$padj < pvalue] <- "yes"
  ggmat$toShow <- "no"
  if(!is.null(genes)) ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = isSignif, label = GENE))
  p <- p + geom_point(size = 1, alpha = 0.25)
  p <- p + theme_bw(base_size = 20)
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("-log10 Pvalue")
  p <- p + theme(legend.position="none")
  p <- p + scale_color_manual(values=c(yes = "red2", no = "grey"))
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), force = 2)
  
  pdf(outFile)
  plot(p)
  dev.off()
}



log2df <- function(logM)
{
  reads.tot <- logM$V2[logM$V1 == "                          Number of input reads |"]
  mapped <- logM$V2[logM$V1 == "                   Uniquely mapped reads number |"]
  mapped.pc <- logM$V2[logM$V1 == "                        Uniquely mapped reads % |"]
  multi <- logM$V2[logM$V1 == "        Number of reads mapped to multiple loci |"]
  multi.pc <- logM$V2[logM$V1 == "             % of reads mapped to multiple loci |"]
  unmapped.mismatches <- logM$V2[logM$V1 == "       % of reads unmapped: too many mismatches |"]
  unmapped.short <- logM$V2[logM$V1 == "                 % of reads unmapped: too short |"]
  unmapped.other <- logM$V2[logM$V1 == "                     % of reads unmapped: other |"]
  chimeric <- logM$V2[logM$V1 == "                            % of chimeric reads |"]
  
  return(data.frame(Reads.TOTAL = reads.tot,
                    Mapped = mapped,
                    Mapped.PC = mapped.pc,
                    MultipleLoci = multi,
                    MultipleLoci.PC = multi.pc,
                    UnMapped.mismatches.PC = unmapped.mismatches,
                    UnMapped.tooShort.PC = unmapped.short,
                    UnMapped.other.pc = unmapped.other,
                    Chimeric.PC = chimeric)
  )
  
}


ggLog <- function(logDF, outFile){
  ggmat <- data.frame(SAMPLE = logDF$SAMPLE,
                      Reads.TOTAL = toNum(logDF[, "Reads.TOTAL"]),
                      Mapped.PC = as.numeric(gsub("%", "", logDF[, "Mapped.PC"])),
                      MultipleLoci.PC = as.numeric(gsub("%", "", logDF[, "MultipleLoci.PC"])),
                      UnMapped.mismatches.PC = as.numeric(gsub("%", "", logDF[, "UnMapped.mismatches.PC"])),
                      UnMapped.tooShort.PC = as.numeric(gsub("%", "", logDF[, "UnMapped.tooShort.PC"])),
                      UnMapped.other.pc = as.numeric(gsub("%", "", logDF[, "UnMapped.other.pc"])),
                      Chimeric.PC = as.numeric(gsub("%", "", logDF[, "Chimeric.PC"]))
  )
  ggmat <- melt(ggmat)
  
  p <- ggplot(ggmat, aes(x = SAMPLE, y = value))
  p <- p + geom_bar(stat="identity")
  p <- p + facet_wrap(~ variable, scales = "free_y")
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(outFile, width = 15, height = 10)
  plot(p)
  dev.off()
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


#################################################################
# LOG 2 DATAFRAME
setwd(file.path("~/cluster/master/Duyster/RNA_010719/quant/"))
logFiles <- list.files(pattern = "_Log.final.out$", recursive = TRUE)
#names(logFiles) <- gsub("_Log.final.out", "", logFiles)
names(logFiles) <- unlist(lapply(strsplit(logFiles, split = "/"), function(i) i[1]))

logList <- lapply(logFiles, read.delim, header = FALSE)

logDF <- lapply(logList, log2df)
logDF <- do.call(rbind, logDF)

logDF$SAMPLE <- rownames(logDF)

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count/"))
write.xlsx(logDF, "RNA_Log.final.out.xlsx", row.names= TRUE)

# LOG 2 GGPLOT
logDF$SAMPLE <- factor(logDF$SAMPLE, levels = sort(logDF$SAMPLE))
ggLog(logDF, "RNA_Log.final.out.pdf")



#################################################################
# GET COUNT MATRIX

# Load count data
setwd(file.path("~/cluster/master/Duyster/RNA_010719/quant/"))
countFiles <- list.files(pattern = "_ReadsPerGene.out.tab", recursive = TRUE)
names(countFiles) <- dirname(countFiles)
countList <- lapply(countFiles, read.delim, skip = 4, header = FALSE, stringsAsFactors = FALSE)

# Get Ensembl
ensembl <- countList[[1]]$V1

# Make count table
countList <- lapply(countList, function(i) return(i[,2]))
countMat <- do.call(cbind, countList)
rownames(countMat) <- ensembl

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(countMat), ann$FASTQ),]# RE-ORDER annotation

colnames(countMat) <- ann$SAMPLE

#setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count/"))
#write.table(countMat, "count.txt", sep = "\t", quote = FALSE)

##########################################################################################
# Load count matrix
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
countMat <- read.delim("count.txt")

# Load annotation

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(countMat), ann$SAMPLE),]# RE-ORDER annotation

# Remove outliers
idx <- which(ann$SAMPLE == "S10")
ann <- ann[-idx,]
countMat <- countMat[, -idx]

##################
# BUILD DESEQ DATA
dds <- DESeqDataSetFromMatrix(countData = countMat,
                              colData = ann,
                              design = ~ 1)

# ADD GENE META DATA
ensID <- rownames(countMat)
entrezID <- ensembl2entrez(ensID)
entrezID[is.na(entrezID)] <- "NA"
symbol <- entrez2symbol(entrezID)
geneName <- entrez2name(entrezID)

featureData <- data.frame(entrez = entrezID, symbol = symbol, ensembl = ensID, gene.name = geneName)
mcols(dds) <- DataFrame(mcols(dds), featureData)

# FILTERING
dds <- dds[entrezID != "NA",]# remove NA entrez IDs

keep <- rowSums(counts(dds)) >= 5 #remove low count
dds <- dds[keep,]

# GET COUNT
dds <- estimateSizeFactors(dds)
dds.counts <- counts(dds, normalized=TRUE)

# GET NORMALIZED
rld <- rlogTransformation(dds, blind=TRUE)
#rld <- vst(dds, blind=TRUE)
rld <- assay(rld)

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
write.table(rld, "rlogTransformation_DESeq2.txt", sep = "\t", quote = FALSE)

deseqDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2")

###################################################################################################
# Group wise; all subset together

# SET DESIGN
Group <- ann$GROUP
Group <- factor(Group)
colData(dds)$Group  <- Group

design(dds) <- formula(~ + Group)

# DIFFERENTIAL ANALYSIS
dds <- DESeq(dds, test = "Wald", fitType = "mean")

# CONTRAST
contMat <- matrix(c("OSM", "CTL"),
                  ncol = 2, byrow = TRUE)
contMat <- cbind("Group", contMat)					 
contName <- apply(contMat, 1, function(i) paste(i[2], i[3], sep = "-"))

# RUN PIPELINE
resList <- apply(contMat, 1, function(i)
  getDEGMatrix(dds, i, rld = NULL)
)
names(resList) <- contName	

# SAVE
setwd(deseqDir)
lapply(1:length(resList), function(i)
  write.xlsx(resList[[i]], paste0(names(resList)[i], "_DESeq2.xlsx"), row.names = FALSE)
)

###################################################################################################
# NB DEG

# Load limma outputs
setwd(deseqDir)
limmaFiles <- list.files(pattern = "_DESeq2.xlsx")
names(limmaFiles) <- gsub("_DESeq2.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

# NUMBER OF DEG
fcCutoff <- 0
pvCutoff <- 0.05
nbList <- lapply(limmaList, getNbDEG, pv = pvCutoff, fc = fcCutoff, adjusted = FALSE)
nbList.adj <- lapply(limmaList, getNbDEG, pv = pvCutoff, fc = fcCutoff, adjusted = TRUE)

nbMat <- data.frame(DESeq2 = names(limmaList), NbDEG = unlist(nbList), NbDEG.adj = unlist(nbList.adj))
write.xlsx(nbMat, paste0("nbDEG_pv", pvCutoff, "_fc", fcCutoff, ".xlsx"), row.names = FALSE)



###################################################################################################
# VOLCANO PLOTS

# Load limma outputs
setwd(deseqDir)
limmaFiles <- list.files(pattern = "_DESeq2.xlsx")
names(limmaFiles) <- gsub("_DESeq2.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

lapply(1:length(limmaList), function(i)
  ggVolcano(limmaList[[i]], pvalue = 0.05, genes = head(limmaList[[i]]$symbol, 25),
            outFile = paste0(names(limmaList)[i], "_volcano.pdf"))
)



