
library(gage)
library(GO.db)
library(org.Mm.eg.db)
library(openxlsx)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(ggplot2)

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

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Mm.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Mm.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}


gagePipeline <- function(test, ref, db, up = NULL, down = NULL)
{
  gageMat <- cbind(ref, test)
  gos <- gage(gageMat, gsets = db, ref = seq(1:ncol(ref)),
              same.dir = T, compare = "unpaired", rank = T,
              saaTest = gs.tTest, set.size = c(5, 500))
  significant.groups <- sigGeneSet(gos, cutoff = 1, qpval = c("q.val"))
  
  # Greater
  greater <- significant.groups$greater
  greater.entrez <- lapply(rownames(greater), function(i) return(intersect(db[[i]], up)))
  greater.symbol <- lapply(greater.entrez, entrez2symbol)
  greater.nb <- lapply(greater.entrez, length)
  
  greater.entrez <- unlist(lapply(greater.entrez, toString))
  greater.symbol <- unlist(lapply(greater.symbol, toString))
  greater.nb <- unlist(greater.nb)
  greater <- cbind(greater, nb = greater.nb, entrez = greater.entrez, symbol = greater.symbol)
  
  # Less
  less <- significant.groups$less
  less.entrez <- lapply(rownames(less), function(i) return(intersect(db[[i]], down)))
  less.symbol <- lapply(less.entrez, entrez2symbol)
  less.nb <- lapply(less.entrez, length)
  
  less.entrez <- unlist(lapply(less.entrez, toString))
  less.symbol <- unlist(lapply(less.symbol, toString))
  less.nb <- unlist(less.nb)
  less <- cbind(less, nb = less.nb, entrez = less.entrez, symbol = less.symbol)
  
  return(list(greater = greater, less = less))
}

gagePipelineUni <- function(test, ref, db, deg = NULL)
{
  gageMat <- cbind(ref, test)
  gos <- gage(gageMat, gsets = db, ref = seq(1:ncol(ref)),
              same.dir = F, compare = "unpaired", rank = T,
              saaTest = gs.tTest, set.size = c(5, 500))
  significant.groups <- sigGeneSet(gos, cutoff = 1, qpval = c("q.val"))
  
  # Greater
  greater <- significant.groups$greater
  greater.entrez <- lapply(rownames(greater), function(i) return(intersect(db[[i]], deg)))
  greater.symbol <- lapply(greater.entrez, entrez2symbol)
  greater.nb <- lapply(greater.entrez, length)
  
  greater.entrez <- unlist(lapply(greater.entrez, toString))
  greater.symbol <- unlist(lapply(greater.symbol, toString))
  greater.nb <- unlist(greater.nb)
  greater <- cbind(greater, nb = greater.nb, entrez = greater.entrez, symbol = greater.symbol)
  
  return(list(greater = greater))
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
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(gs[,1], trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"q.val")		
    else(m[idx, i] <- as.numeric(gs$"p.val"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

gageHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE)
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
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  idxMat <- idxMat[1:nb, ]
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
           cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
           width = 12, height = 12)
}


gageComplexHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat <- gsMat[, c(idx.up, idx.down)]
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  gsMat <- -log10(gsMat)	
  gsMat[gsMat > 5] <- 5 # set the limit to 5
  
  #rownames(gsMat) <- gsub("^HALLMARK_", "", rownames(gsMat))
  
  # Column annotation
  ha.df <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
  ha <- HeatmapAnnotation(df = ha.df, col = list(Sign = c("UP" = "red", "DOWN" = "blue")), show_legend = FALSE)
  
  
  # Split UP / DOWN gene-sets
  mySplit <- c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down)))
  mySplit <- factor(mySplit, levels = c("UP", "DOWN"))
  
  # size of the cells
  myheight <- 0.5*nrow(gsMat)
  mywidth <- 0.5*ncol(gsMat)
  
  ht <- Heatmap(gsMat,
                col = magma(10),
                heatmap_legend_param = list(title = "",
                                            legend_height = unit(4, "cm")),
                clustering_distance_rows = "euclidean",
                #cluster_columns = dend,
                show_column_names = TRUE,
                top_annotation = ha,
                column_split = mySplit,
                cluster_columns = FALSE,
                rect_gp = gpar(col = "grey", lwd = 1),
                row_names_gp = gpar(fontsize = 10),
                row_names_side = "left", #row_dend_side = "right",
                show_row_dend = FALSE,
                width = unit(mywidth, "cm"), height = unit(myheight, "cm")
  )
  
  pdf(outFile, width = 12)
  draw(ht, padding = unit(c(2, 150, 2, 2), "mm"), row_sub_title_side = "right")
  dev.off()
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# CREATE A CUSTOM GENE-SET
#setwd(file.path("/Volumes/Home/Geoffroy/Bengsch/RNA_190619/doc/Bertram_email_310619"))
#geneMat <- read.xlsx("ESG_summary_human_transcriptomic only.xlsx", sheet = 1)
#geneList <- apply(geneMat, 2, list)
#geneList <- lapply(geneList, function(i) i[[1]])
#geneList <- lapply(geneList, function(i) i[!is.na(i)])
#geneList <- lapply(geneList, symbol2entrez)
#dbList <-list(CUSTOM = geneList)

# EXPRESSION
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
expMat <- read.delim("rlogTransformation_DESeq2.txt")

expMat.entrez <- ensembl2entrez(rownames(expMat))
expMat.entrez <- expMat.entrez[!is.na(expMat.entrez)]
expMat.entrez <- expMat.entrez[!duplicated(expMat.entrez)]

expMat <- expMat[names(expMat.entrez), ]
rownames(expMat) <- expMat.entrez


# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(expMat), ann$SAMPLE), ]

Group <- ann$GROUP
cond <- Group


# CONTRAST
contMat <- matrix(c("OSM", "CTL"),
  ncol = 2, byrow = TRUE)
contMat <- cbind("Group", contMat)					 
contName <- apply(contMat, 1, function(i) paste(i[2], i[3], sep = "-"))


deseqFolder <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/")
gageFolder <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/gsea/gage/")

######################
lapply(1:length(dbList), function(i){
  
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  lapply(1:nrow(contMat), function(j){
    testCond <- contMat[j, 2]
    refCond <- contMat[j, 3]
    
    deseqMat <- read.xlsx(file.path(deseqFolder, paste0(contName[j], "_DESeq2.xlsx")), sheet=1)
    deseqFC <- deseqMat$log2FoldChange
    deseqUP <- deseqMat$entrez[deseqFC > 0]
    deseqDOWN <- deseqMat$entrez[deseqFC < 0]
    
    testMat <- expMat[,cond==testCond]
    refMat <- expMat[,cond==refCond, drop = FALSE]
    gage.res <- gagePipeline(testMat, refMat, mydb, deseqUP, deseqDOWN)
    
    dir.create(file.path(gageFolder, mydb.name), showWarnings = FALSE)
    
    write.xlsx(list(UP = gage.res$greater, DOWN = gage.res$less),
               file.path(gageFolder, mydb.name, paste(contName[j], "_", mydb.name, "_gage.xlsx", sep = "")), row.names = TRUE)
    
  })
  
})

##########
# HEATMAPS

setwd(gageFolder)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]

  fhFiles <- c(file.path(mydb.name, paste0("OSM-CTL_", mydb.name, "_gage.xlsx")))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_", mydb.name, "_gage.xlsx"), "", fhFiles)
  fhList <- do.call(c, fhList)
  
  gageHeatmap(fhList, file.path(mydb.name, paste0("OSM-CTL_", mydb.name, "_gage_heatmap.pdf")), nb = 20, adjusted = FALSE)
  
})


#######################################
# OVERALL BARPLOT

# Load all databases for a given comparison
setwd(gageFolder)
nb <- 50
comp <- c("OSM-CTL")
dbList.sub <- dbList[-7]# remove CGP



lapply(comp, function(i){
  fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_", j, "_gage.xlsx")))
  names(fhFiles) <- paste0(i, "_", names(fhFiles))
  fhList.up <- lapply(fhFiles, read.xlsx, sheet = "UP")
  fhList.down <- lapply(fhFiles, read.xlsx, sheet = "DOWN")
  
  # UP
  fh.up <- do.call(rbind, fhList.up)
  fh.up$"p.val" <- as.numeric(fh.up$"p.val")
  fh.up$"set.size" <- as.numeric(fh.up$"set.size")
  fh.up <- fh.up[order(fh.up$"p.val"), ]
  colnames(fh.up)[1] <- "Term"
  
  ggmat.up <- fh.up[1:nb, c("Term", "p.val", "set.size")]
  ggmat.up$Score <- -log10(ggmat.up$"p.val")
  ggmat.up$Term <- substring(ggmat.up$Term, 1, 60)
  ggmat.up$Term <- factor(ggmat.up$Term, levels = rev(ggmat.up$Term[order(ggmat.up$"p.val")]))
  
  p <- ggplot(data=ggmat.up, aes(x=Term, y=Score, fill = set.size))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("-log10 pvalue")
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "UP_overall_barplot.pdf"), width = 8, height = 10)
  
  # down
  fh.down <- do.call(rbind, fhList.down)
  fh.down$"p.val" <- as.numeric(fh.down$"p.val")
  fh.down$"set.size" <- as.numeric(fh.down$"set.size")
  fh.down <- fh.down[order(fh.down$"p.val"), ]
  colnames(fh.down)[1] <- "Term"
  
  ggmat.down <- fh.down[1:nb, c("Term", "p.val", "set.size")]
  ggmat.down$Score <- -log10(ggmat.down$"p.val")
  ggmat.down$Term <- substring(ggmat.down$Term, 1, 60)
  ggmat.down$Term <- factor(ggmat.down$Term, levels = rev(ggmat.down$Term[order(ggmat.down$"p.val")]))
  
  p <- ggplot(data=ggmat.down, aes(x=Term, y=Score, fill = set.size))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("-log10 pvalue")
 
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "DOWN_overall_barplot.pdf"), width = 8, height = 10)
  
})



