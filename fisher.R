
library(ggplot2)
library(genefilter)
library(data.table)
library(pheatmap)
library(GO.db)
library(org.Mm.eg.db)
library(openxlsx)
library(viridis)

source("/Users/gandrieux/bitbucket/work/tools/venn_script.r")
source("/Users/gandrieux/bitbucket/work/tools/iwanthue.r")
source("/Users/gandrieux/bitbucket/work/tools/hyperG_parallel.R")

# Consensus Path DB
load("/Volumes/Home/Geoffroy/database/consensus/Consensus_unique08.mouse.RData") # load cons2 annotation list
# Misgdb 
load("/Volumes/Home/Geoffroy/database/msigdb_mouse/msigdb_v6.2.mouse.RData")

dbList <- list("Reactome" = reactomedb, "KEGG" = keggdb,
               "goBP" = goBP, "goMF" = goMF, "goCC" = goCC,
               "Hallmark" = hallmarkdb, "CGP" = cgpdb,
               "Consensus" = consdb)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Mm.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}


toNum <- function(x) return(as.numeric(levels(x))[x])


saveFisherDetails <- function(fisherList, outFile)
{
  fisherList.size <- sapply(fisherList, nrow)
  fisherList <- fisherList[fisherList.size != 0]
  
  fh_info <- data.frame(Sheet = paste("sheet_", seq(1, length(fisherList)), sep = ""),
                        Group = names(fisherList))
  
  toXLSX <- list("info" = fh_info)
  toXLSX <- c(toXLSX, fisherList)
  names(toXLSX)[2:length(toXLSX)] <- paste0("sheet_", seq(1, length(fisherList)))
  
  write.xlsx(toXLSX, outFile, row.names = FALSE)	
}

FisherFromVenn <- function(venn, db, org.library, univ = NULL)
{
  if(is.null(univ)) univ <- venn$entrez
  groups <- unique(venn$Group)
  
  vList <- lapply(groups, function(i)
    hyperG(venn$entrez[venn$Group == i], geneSets = db, universe = univ, org.library = org.library, cutoff = 1, mincount = 2)
  )
  names(vList) <- groups
  return(vList)
}


my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}

getHeatmapMatrix <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i$Term))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(gs$Term, trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"adj.P.Val")		
    else(m[idx, i] <- as.numeric(gs$"p-value"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

fisherHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat <- gsMat[, c(idx.up, idx.down)]
  
  # add annotation
  ann.col <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
  rownames(ann.col) <- colnames(gsMat)
  myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"))
  
  #gsMat.min <- rowMin(gsMat)
  #gsMat <- gsMat[order(gsMat.min),, drop = FALSE]	
  #if(nrow(gsMat) > nb) gsMat <- gsMat[1:nb, ]
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  if(nrow(idxMat) > nb) idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  gsMat <- -log10(gsMat)	
  gsMat[gsMat > 5] <- 5 # set the limit to 5
  
  paletteLength <- 10
  myColor <- magma(paletteLength)
  #myColor <- viridis(paletteLength)
  #myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- c(seq(0, max(gsMat), length.out=ceiling(paletteLength)))
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
           cluster_cols = FALSE, cluster_rows = doClust,
           annotation_col = ann.col, annotation_colors = myColor.ann,
           cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8
           #, width = 10, height = 10
  )
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

####################
# FISHER FROM DESEQ2

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/"))
deseqFiles <- list.files(pattern = "_DESeq2.xlsx")
names(deseqFiles) <- gsub("_DESeq2.xlsx", "", deseqFiles)
deseqList <- lapply(deseqFiles, read.xlsx, sheet = 1)

deseqList <- lapply(deseqList, function(i){
  i$padj[is.na(i$padj)] <- 1
  return(i)
})

genes <- deseqList[[1]]$entrez
genes <- unique(genes[!is.na(genes)])

pvCutoff <- 0.05
fcCutoff <- 0
deseqList.UP <- lapply(deseqList, function(i) i$entrez[i$padj < pvCutoff & i$log2FoldChange > fcCutoff])
deseqList.DOWN <- lapply(deseqList, function(i) i$entrez[i$padj < pvCutoff & i$log2FoldChange < -fcCutoff])
deseqList.DEG <- lapply(deseqList, function(i) i$entrez[i$padj < pvCutoff & abs(i$log2FoldChange) > fcCutoff])

# Fisher
fisherDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/gsea/Fisher/")
setwd(fisherDir)


lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  dir.create(mydb.name, showWarnings = FALSE)
  
  
  fh.UP <- lapply(deseqList.UP, hyperG, geneSets = mydb, universe = genes, org.library = "org.Mm.eg.db", cutoff = 1, mincount = 2)
  fh.DOWN <- lapply(deseqList.DOWN, hyperG, geneSets = mydb, universe = genes, org.library = "org.Mm.eg.db", cutoff = 1, mincount = 2)
  fh.DEG <- lapply(deseqList.DEG, hyperG, geneSets = mydb, universe = genes, org.library = "org.Mm.eg.db", cutoff = 1, mincount = 2)
  
  # Save
  for(j in 1:length(fh.UP))
  {
    setwd(fisherDir)
    #dir.create(names(fh.UP)[j], showWarnings = FALSE)
    #setwd(names(fh.UP)[j])
    fhName <- file.path(mydb.name,
                        paste(names(fh.UP)[j], "_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx", sep = ""))
    write.xlsx(list(UP = fh.UP[[j]], DOWN = fh.DOWN[[j]], DEG = fh.DEG[[j]]), fhName, row.names = FALSE)	
  }
})


#######################################
# HEATMAPS

setwd(fisherDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # biopsy1 vs. biopsy2 vs. biopsy3 vs. HD
  fhFiles <- c(file.path(mydb.name, paste0("OSM-CTL_pv0.05_fc0_", mydb.name, "_Fisher.xlsx")))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_pv0.05_fc0_", mydb.name, "_Fisher.xlsx"), "", fhFiles)
  fhList <- do.call(c, fhList)
  fhList <- fhList[-grep("DEG", names(fhList))]
  
  fisherHeatmap(fhList, file.path(mydb.name,
                                  paste0("OSM-CTL_", mydb.name, "_Fisher_heatmap.pdf")), nb = 20, adjusted = FALSE)
  
})

#######################################
# OVERALL BARPLOT

# Load all databases for a given comparison
setwd(fisherDir)
nb <- 50
comp <- c("OSM-CTL")
dbList.sub <- dbList[-7]# remove CGP

lapply(comp, function(i){
  fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv0.05_fc0_", j, "_Fisher.xlsx")))
  names(fhFiles) <- paste0(i, "_", names(fhFiles))
  fhList.up <- lapply(fhFiles, read.xlsx, sheet = "UP")
  fhList.down <- lapply(fhFiles, read.xlsx, sheet = "DOWN")

  # UP
  fh.up <- do.call(rbind, fhList.up)
  fh.up <- fh.up[order(fh.up$"p-value"), ]
  
  ggmat.up <- fh.up[1:nb, c("Term", "p-value", "Count")]
  ggmat.up$Score <- -log10(ggmat.up$"p-value")
  ggmat.up$Term <- substring(ggmat.up$Term, 1, 60)
  ggmat.up$Term <- factor(ggmat.up$Term, levels = rev(ggmat.up$Term[order(ggmat.up$"p-value")]))
  
  p <- ggplot(data=ggmat.up, aes(x=Term, y=Score, fill = Count))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("-log10 pvalue")
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "UP_overall_barplot.pdf"), width = 8, height = 10)
  
  # DOWN
  fh.down <- do.call(rbind, fhList.down)
  fh.down <- fh.down[order(fh.down$"p-value"), ]
  
  ggmat.down <- fh.down[1:nb, c("Term", "p-value", "Count")]
  ggmat.down$Score <- -log10(ggmat.down$"p-value")
  ggmat.down$Term <- substring(ggmat.down$Term, 1, 60)
  ggmat.down$Term <- factor(ggmat.down$Term, levels = rev(ggmat.down$Term[order(ggmat.down$"p-value")]))
  
  p <- ggplot(data=ggmat.down, aes(x=Term, y=Score, fill = Count))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("-log10 pvalue")
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "DOWN_overall_barplot.pdf"), width = 8, height = 10)
  
})






