
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(openxlsx)
library(ComplexHeatmap)
library(viridis)
library(GeneAnswers)
library(circlize)
library(measurements)

# Consensus Path DB
load("/Volumes/Home/Geoffroy/database/consensus/Consensus_unique08.mouse.RData") # load cons2 annotation list
# Misgdb 
load("/Volumes/Home/Geoffroy/database/msigdb_mouse/msigdb_v6.2.mouse.RData")

dbList <- list("Reactome" = reactomedb, "KEGG" = keggdb,
               "goBP" = goBP, "goMF" = goMF, "goCC" = goCC,
               "Hallmark" = hallmarkdb, "CGP" = cgpdb,
               "Consensus" = consdb)


source("~/bitbucket/work/tools/iwanthue.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getGeneSetsFromKeyword <- function(keyword, dbList){
  return(lapply(dbList, function(i){
    sort(grep(keyword, names(i), ignore.case = TRUE, value = TRUE))
  }))
} 

getGeneSetsFromGene <- function(gene, dbList, isEID = TRUE){
  if(!isEID) gene <- symbol2entrez(gene)
  return(lapply(dbList, function(i){
    isInside <- sapply(i, function(y) gene %in% y)
    return(sort(names(i)[isInside]))
  }))
} 

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

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Mm.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

getGenesDESeq2 <- function(deseq, fcCt, pvCt, dsign, adjusted)
{
  deseqFC <- deseq$log2FoldChange
  ifelse(adjusted, deseqPV <- deseq$padj, deseqPV <- deseq$pvalue)
  deseqPV[is.na(deseqPV)] <- 1
  
  deseqUP <- as.character(deseq$entrez)[deseqFC > fcCt & deseqPV <= pvCt]
  deseqUP <- deseqUP[!deseqUP=="NA"]
  deseqDOWN <- as.character(deseq$entrez)[deseqFC < -fcCt & deseqPV <= pvCt]
  deseqDOWN <- deseqDOWN[!deseqDOWN=="NA"]
  if(dsign == "up") return(deseqUP)
  else if(dsign == "down") return(deseqDOWN)
  else(return(c(deseqUP, deseqDOWN)))
}


scaledAnnoHeatmapLabel <- function(mat, ann.row, ann.col, outFile)
{
  paletteLength <- 25
  myMax <- ceiling(max(abs(mat)))
  
  myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  myBreaks <- myBreaks[myBreaks != 0]
  #myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
  
  #myColor.ann <- list(Group = c(Dc = "black", Tc = "grey", WEHI = "red"))
  
  if(nrow(mat)<2) doCluster <- FALSE
  else(doCluster <- TRUE)
  
  pheatmap(mat, color = myColor, breaks = myBreaks, filename = outFile,
           annotation_row = ann.row, annotation_col = ann.col,
           #annotation_colors = myColor.ann,
           cluster_cols = FALSE, cluster_rows = doCluster, show_rownames = TRUE,
           cellwidth = 12, cellheight = 12
           #height = 10, width = 10
  )
  
}	

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

complexHeatmapLabel <- function(mat, outName, annCol = NULL, colCol = NULL, fc = NULL, pch = NULL, colFct = NULL){
  myLimit <- max(abs(mat), na.rm = TRUE)
  
  if(is.null(colFct)){
    col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("dodgerblue", "white", "firebrick1"))}else col_fun <- colFct
    
    if(!is.null(annCol) & is.null(colCol)){
      if(ncol(annCol) != 1){
        colCol <- apply(annCol, 2, getColor)
      }else{
        colCol <- list(getColor(annCol[,1]))
        names(colCol) <- colnames(annCol)[1]
        }
    
    }  
    
    if(!is.null(fc)){
      fcLimit <- floor(max(abs(fc)))
      #colRow_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("springgreen", "black", "red3"))
      colRow_fun <- colorRamp2(c(-fcLimit, 0, fcLimit), c("dodgerblue", "white", "firebrick1"))
      
      lgd_fc = Legend(title = "fc", col = colRow_fun, at = c(-fcLimit, 0, fcLimit), 
                      labels = c(-fcLimit, 0, fcLimit))
    }
    
    if(is.null(pch)) pch <- rep(NA, nrow(mat))  
    
    # size of the cells
    myheight <- 0.4*nrow(mat)
    mywidth <- 0.4*ncol(mat)
    
    ha <- Heatmap(mat, name = "mat", col = col_fun, na_col = "black",
                  #clustering_distance_rows = robust_dist, row_dend_reorder = TRUE,
                  cluster_columns = FALSE,
                  row_dend_side = "left",
                  #cluster_rows = FALSE,
                  row_names_gp = gpar(fontsize = 12),
                  top_annotation = HeatmapAnnotation(df = annCol, col = colCol),
                  right_annotation = rowAnnotation(fc = anno_simple(fc, col = colRow_fun, pch = pch)),
                  width = unit(mywidth, "cm"), height = unit(myheight, "cm")
    )
    
    # save pdf
    pdfHeight <- conv_unit(as.numeric(ComplexHeatmap:::height(draw(ha))), "mm", "inch")+1
    pdfWidth <- conv_unit(as.numeric(ComplexHeatmap:::width(draw(ha))), "mm", "inch")+1
    
    pdf(outName, width = pdfWidth, height = pdfHeight)
    draw(ha, annotation_legend_list = list(lgd_fc))
    dev.off()
    
}

complexHeatmapNOLabel <- function(mat, outName, annCol = NULL, colCol = NULL, fc = NULL, pch = NULL, colFct = NULL){
  myLimit <- max(abs(mat), na.rm = TRUE)
  
  if(is.null(colFct)){
    col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("dodgerblue", "white", "firebrick1"))}else col_fun <- colFct
    
    if(!is.null(annCol) & is.null(colCol)){
      if(ncol(annCol) != 1){
        colCol <- apply(annCol, 2, getColor)
      }else{
        colCol <- list(getColor(annCol[,1]))
        names(colCol) <- colnames(annCol)[1]
      }
      
    }  
    
    if(!is.null(fc)){
      fcLimit <- floor(max(abs(fc)))
      #colRow_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("springgreen", "black", "red3"))
      colRow_fun <- colorRamp2(c(-fcLimit, 0, fcLimit), c("dodgerblue", "white", "firebrick1"))
      
      lgd_fc = Legend(title = "fc", col = colRow_fun, at = c(-fcLimit, 0, fcLimit), 
                      labels = c(-fcLimit, 0, fcLimit))
    }
    
    # size of the cells
    mywidth <- 0.4*ncol(mat)
    
    if(!(is.null(fc))){
      ha <- Heatmap(mat, name = "mat", col = col_fun, na_col = "black",
                    #clustering_distance_rows = robust_dist, row_dend_reorder = TRUE,
                    cluster_columns = FALSE, row_dend_side = "left",
                    show_row_names = FALSE,
                    top_annotation = HeatmapAnnotation(df = annCol, col = colCol),
                    right_annotation = rowAnnotation(fc = anno_simple(fc, col = colRow_fun)),
                    width = unit(mywidth, "cm"), heatmap_height = unit(9, "cm")
      )
      
      # save pdf
      pdfHeight <- conv_unit(as.numeric(ComplexHeatmap:::height(draw(ha))), "mm", "inch")+1
      pdfWidth <- conv_unit(as.numeric(ComplexHeatmap:::width(draw(ha))), "mm", "inch")+1
      
      pdf(outName, width = pdfWidth, height = pdfHeight)
      draw(ha, annotation_legend_list = list(lgd_fc))
      dev.off()
      
      
    } else{
      ha <- Heatmap(mat, name = "mat", col = col_fun, na_col = "black",
                    #clustering_distance_rows = robust_dist, row_dend_reorder = TRUE,
                    cluster_columns = FALSE, row_dend_side = "left",
                    show_row_names = FALSE,
                    top_annotation = HeatmapAnnotation(df = annCol, col = colCol),
                    width = unit(mywidth, "cm"), heatmap_height = unit(9, "cm")
      )
      
      # save pdf
      pdfHeight <- conv_unit(as.numeric(ComplexHeatmap:::height(draw(ha))), "mm", "inch")+1
      pdfWidth <- conv_unit(as.numeric(ComplexHeatmap:::width(draw(ha))), "mm", "inch")+1
      
      pdf(outName, width = pdfWidth, height = pdfHeight)
      draw(ha)
      dev.off()
    }
    
    
    
    
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

###################################################################
# BUILD CUSTOM GENE SETS (Michael Rassner 30.08.19)

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc/"))
ex.all <- read.xlsx("T_cell_exhaustion.xlsx", sheet = "ALL")
ex.up <- read.xlsx("T_cell_exhaustion.xlsx", sheet = "UP")
ex.down <- read.xlsx("T_cell_exhaustion.xlsx", sheet = "DOWN")

# Convert to list
ex.all <- lapply(as.list(ex.all), function(i) i[!is.na(i)])
ex.all <- lapply(ex.all, function(i) unique(gsub(" ", "", i)))
ex.all <- lapply(ex.all, alias2Symbol, species = "Mm")
ex.all <- lapply(ex.all, symbol2entrez)

ex.up <- lapply(as.list(ex.up), function(i) i[!is.na(i)])
ex.up <- lapply(ex.up, function(i) unique(gsub(" ", "", i)))
ex.up <- lapply(ex.up, alias2Symbol, species = "Mm")
ex.up <- lapply(ex.up, symbol2entrez)

ex.down <- lapply(as.list(ex.down), function(i) i[!is.na(i)])
ex.down <- lapply(ex.down, function(i) unique(gsub(" ", "", i)))
ex.down <- lapply(ex.down, alias2Symbol, species = "Mm")
ex.down <- lapply(ex.down, symbol2entrez)

dbList <- list(exhaustion_all = ex.all,
               exhaustion_up = ex.up,
               exhaustion_down = ex.down)



################################################################################
# TOP XX DEG (DESEQ2)

# LOAD EXPRESSION
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
expMat <- read.delim("rlogTransformation_DESeq2.txt")

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(expMat), ann$SAMPLE), ]

# Remove outliers
#idx <- which(ann$SAMPLE == "S10")
#ann <- ann[-idx,]
#countMat <- countMat[, -idx]

# RE-ORDER ANNOTATION
ann <- ann[order(ann$GROUP, ann$SAMPLE), ]
expMat <- expMat[, ann$SAMPLE]

cond <- ann$GROUP


# CONTRAST
contMat <- matrix(c("OSM", "CTL"),
                  ncol = 2, byrow = TRUE)
contMat <- cbind("Group", contMat)					 
contName <- apply(contMat, 1, function(i) paste(i[2], i[3], sep = "-"))

deseqDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/")
heatmapDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/complexHeatmaps/")



# Select top XX genes per sample
nb <- 50
fcCutoff <- 1
pvCutoff <- 0.001


lapply(1:nrow(contMat), function(i){
  
  testCond <- contMat[i,2]
  refCond <- contMat[i,3]
  
  # load deseq
  deseqMat <- read.xlsx(file.path(deseqDir, paste0(contName[i], "_DESeq2.xlsx")))
  
  # Select significanct and re-order according to the foldchange
  deseqMat <- deseqMat[deseqMat$padj < pvCutoff & abs(deseqMat$log2FoldChange) > fcCutoff, ]
  deseqMat <- na.omit(deseqMat)
  deseqMat <- deseqMat[order(-abs(deseqMat$log2FoldChange)), ]
  
  # get DEG
  deg <- deseqMat$ensembl[deseqMat$padj < pvCutoff & abs(deseqMat$log2FoldChange) > fcCutoff]
  deg <- deg[!is.na(deg)]
  
  # get columns
  mysamples <- c(ann$SAMPLE[cond == refCond], ann$SAMPLE[cond == testCond])
  
  expMat.sub <- expMat[deg, mysamples]
  
  # ensembl to symbol
  ensembl <- rownames(expMat.sub)
  symbols <- deseqMat$symbol[match(ensembl, deseqMat$ensembl)]
  symbols[is.na(symbols)] <- ensembl[is.na(symbols)]
  symbols[duplicated(symbols)] <- ensembl[duplicated(symbols)]
  rownames(expMat.sub) <- symbols
  
  # heatmap COLUMN annotation
  phAnn.col <- ann[match(colnames(expMat.sub), ann$SAMPLE), c("GROUP"), drop = FALSE]
  rownames(phAnn.col) <- colnames(expMat.sub)
  
  # heatmap ROW annotation
  fc <- deseqMat$log2FoldChange[match(deg, deseqMat$ensembl)]
  pch <- ifelse(ensembl %in% deg, "*", "")
  
  if(nrow(expMat.sub) < nb){
    complexHeatmapLabel(t(scale(t(expMat.sub))),
                        outName = file.path(heatmapDir, paste0(contName[i], "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf")),
                        annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
  }else{
    complexHeatmapLabel(t(scale(t(expMat.sub[1:nb,]))),
                        outName = file.path(heatmapDir, paste0(contName[i], "_pv", pvCutoff, "_fc", fcCutoff, "_top", nb, "_heatmap.pdf")),
                        annCol = phAnn.col, colCol = NULL, fc = fc[1:nb], pch = pch[1:nb], colFct = NULL)
    
    
    complexHeatmapNOLabel(t(scale(t(expMat.sub))),
                          outName = file.path(heatmapDir, paste0(contName[i], "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf")),
                          annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
  }
  
})


#################################################################
# SELECTED GENE-SETS

# BUILD GENE LISTS
getGeneSetsFromKeyword("cytokine", dbList)
getGeneSetsFromKeyword("OXID", dbList)
getGeneSetsFromKeyword("glycol", dbList)
getGeneSetsFromKeyword("lymp", dbList)
getGeneSetsFromKeyword("differentiation", dbList)
getGeneSetsFromKeyword("inflam", dbList)


getGeneSetsFromGene("Havcr2", dbList, isEID = FALSE)
getGeneSetsFromGene("Pdcd1", dbList, isEID = FALSE)
getGeneSetsFromGene("Ctla4", dbList, isEID = FALSE)
getGeneSetsFromGene("Lag3", dbList, isEID = FALSE)

a <- list(unlist(getGeneSetsFromGene("Havcr2", dbList, isEID = FALSE)),
          unlist(getGeneSetsFromGene("Pdcd1", dbList, isEID = FALSE)),
          unlist(getGeneSetsFromGene("Ctla4", dbList, isEID = FALSE)),
          unlist(getGeneSetsFromGene("Lag3", dbList, isEID = FALSE)))


# cytokine
myGS <- getGeneSetsFromKeyword("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", dbList)
myGS.length <- sapply(myGS, length)
myGS <- myGS[myGS.length != 0]
dbList.sub <- lapply(1:length(myGS), function(i){
    db <- dbList[[names(myGS)[i]]][myGS[[i]]]
    return(db)
})
names(dbList.sub) <- names(myGS)


# LOAD EXPRESSION
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
expMat <- read.delim("rlogTransformation_DESeq2.txt")

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(expMat), ann$SAMPLE), ]

# Remove outliers
#idx <- which(ann$SAMPLE == "S10")
#ann <- ann[-idx,]
#countMat <- countMat[, -idx]

# RE-ORDER ANNOTATION
ann <- ann[order(ann$GROUP, ann$SAMPLE), ]
expMat <- expMat[, ann$SAMPLE]

cond <- ann$GROUP


# LOAD DESEQ2 AND ORDER EXPRESSION ACCORDING TO PVALUE
deseqDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/")
heatmapDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/complexHeatmaps/")


# load deseq
deseqMat <- read.xlsx(file.path(deseqDir, "OSM-CTL_DESeq2.xlsx"))

# re-order according to the pvalue
#expMat <- expMat[match(deseqMat$ensembl, rownames(expMat)), ]

# ensembl to symbol
ensembl <- rownames(expMat)
symbols <- deseqMat$symbol[match(ensembl, deseqMat$ensembl)]
symbols[is.na(symbols)] <- ensembl[is.na(symbols)]
symbols[duplicated(symbols)] <- ensembl[duplicated(symbols)]
rownames(expMat) <- symbols

# get DEG
fcCutoff <- 0
pvCutoff <- 0.05
deg <- deseqMat$symbol[deseqMat$padj < pvCutoff & abs(deseqMat$log2FoldChange) > fcCutoff]
deg <- deg[!is.na(deg)]

# Select top XX genes per sample
nb <- 50


lapply(1:length(dbList.sub), function(i){
  
  lapply(1:length(dbList.sub[[i]]), function(j){
    mygenes <- dbList.sub[[i]][[j]]
    mygenes.symbol <- entrez2symbol(mygenes)
    mygenes.symbol <- intersect(rownames(expMat), mygenes.symbol)
    
    if(length(mygenes.symbol) != 0){
      expMat.sub <- expMat[mygenes.symbol, ]
      
      # heatmap COLUMN annotation
      phAnn.col <- ann[match(colnames(expMat.sub), ann$SAMPLE), c("GROUP"), drop = FALSE]
      rownames(phAnn.col) <- colnames(expMat.sub)
      
      # heatmap ROW annotation
      fc <- deseqMat$log2FoldChange[match(mygenes.symbol, deseqMat$symbol)]
      pch <- ifelse(mygenes.symbol %in% deg, "*", "")
      
      new_order <- order(deseqMat$pvalue[match(mygenes.symbol, deseqMat$symbol)])
      expMat.sub <- expMat.sub[new_order, ]
      fc <- fc[new_order]
      pch <- pch[new_order]
      
      outDir <- file.path(heatmapDir, names(dbList.sub)[[i]])
      dir.create(outDir, showWarnings = FALSE)
      outName <- file.path(outDir, paste0(names(dbList.sub[[i]])[j]))
      
      print(outName)
      print(dim(expMat.sub))
      if(nrow(expMat.sub) < nb & nrow(expMat.sub) > 2){
        complexHeatmapLabel(t(scale(t(expMat.sub))),
                            outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf"),
                            annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
      }else if(nrow(expMat.sub) > 2){
        complexHeatmapLabel(t(scale(t(expMat.sub[1:nb,]))),
                            outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_top", nb, "_heatmap.pdf"),
                            annCol = phAnn.col, colCol = NULL, fc = fc[1:nb], pch = pch[1:nb], colFct = NULL)
        
        
        #complexHeatmapNOLabel(t(scale(t(expMat.sub))),
        #                     outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf"),
        #                      annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
      }
      
    }
    
  })
  
})


##############################################################################
# SELECTED GENES

# LOAD TRANSCRIPTION FACTOR LIST
tfMat <- read.delim(file.path("/Volumes/Home/Geoffroy/database/transcription_factors/mouse/AnimalTFDB3.0_Mus_musculus_TF.txt"))
tf <- unique(tfMat$Symbol)


# LOAD EXPRESSION
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count"))
expMat <- read.delim("rlogTransformation_DESeq2.txt")

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann <- ann[match(colnames(expMat), ann$SAMPLE), ]

# Remove outliers
#idx <- which(ann$SAMPLE == "S10")
#ann <- ann[-idx,]
#countMat <- countMat[, -idx]

# RE-ORDER ANNOTATION
ann <- ann[order(ann$GROUP, ann$SAMPLE), ]
expMat <- expMat[, ann$SAMPLE]

cond <- ann$GROUP


# LOAD DESEQ2 AND ORDER EXPRESSION ACCORDING TO PVALUE
deseqDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/")
heatmapDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/complexHeatmaps/")


# load deseq
deseqMat <- read.xlsx(file.path(deseqDir, "OSM-CTL_DESeq2.xlsx"))

# re-order according to the pvalue
#expMat <- expMat[match(deseqMat$ensembl, rownames(expMat)), ]

# ensembl to symbol
ensembl <- rownames(expMat)
symbols <- deseqMat$symbol[match(ensembl, deseqMat$ensembl)]
symbols[is.na(symbols)] <- ensembl[is.na(symbols)]
symbols[duplicated(symbols)] <- ensembl[duplicated(symbols)]
rownames(expMat) <- symbols

# get DEG
fcCutoff <- 0
pvCutoff <- 0.05
deg <- deseqMat$symbol[deseqMat$padj < pvCutoff & abs(deseqMat$log2FoldChange) > fcCutoff]
deg <- deg[!is.na(deg)]

# Select top XX genes per sample
nb <- 100

mygenes.symbol <- tf
mygenes.symbol <- intersect(rownames(expMat), mygenes.symbol)

if(length(mygenes.symbol) != 0){
  expMat.sub <- expMat[mygenes.symbol, ]
  
  # heatmap COLUMN annotation
  phAnn.col <- ann[match(colnames(expMat.sub), ann$SAMPLE), c("GROUP"), drop = FALSE]
  rownames(phAnn.col) <- colnames(expMat.sub)
  
  # heatmap ROW annotation
  fc <- deseqMat$log2FoldChange[match(mygenes.symbol, deseqMat$symbol)]
  pch <- ifelse(mygenes.symbol %in% deg, "*", NA)
  
  new_order <- order(deseqMat$pvalue[match(mygenes.symbol, deseqMat$symbol)])
  expMat.sub <- expMat.sub[new_order, ]
  fc <- fc[new_order]
  pch <- pch[new_order]
  
  outDir <- file.path(heatmapDir, "TF")
  dir.create(outDir, showWarnings = FALSE)
  outName <- file.path(outDir, "AnimalTFDB3.0_Mouse")
  
  print(outName)
  if(nrow(expMat.sub) < nb){
    complexHeatmapLabel(t(scale(t(expMat.sub))),
                        outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf"),
                        annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
  }else{
    complexHeatmapLabel(t(scale(t(expMat.sub[1:nb,]))),
                        outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_top", nb, "_heatmap.pdf"),
                        annCol = phAnn.col, colCol = NULL, fc = fc[1:nb], pch = pch[1:nb], colFct = NULL)
    
    
    #complexHeatmapNOLabel(t(scale(t(expMat.sub))),
    #                     outName = paste0(outName, "_pv", pvCutoff, "_fc", fcCutoff, "_heatmap.pdf"),
    #                      annCol = phAnn.col, colCol = NULL, fc = fc, pch = pch, colFct = NULL)
  }
  
}


