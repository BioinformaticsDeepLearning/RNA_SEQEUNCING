
library(org.Mm.eg.db)
library(openxlsx)
library(pheatmap)
library(GO.db)
library(fgsea)
library(viridis)
library(ggplot2)
library(limma)

source("/Users/gandrieux/bitbucket/work/tools/fgsea_fct.r")
source("/Users/gandrieux/bitbucket/work/tools/hyperG.R")


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
  return(entrez)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




getHeatmapMatrix <- function(gsea_list)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    m[idx, i] <- as.numeric(gs$"NES")
  }
  m[is.na(m)] <- 0	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

getHeatmapMatrix.pv <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"pval")		
    else(m[idx, i] <- as.numeric(gs$"padj"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

fgseaHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE)
{
  gsMat <- getHeatmapMatrix.pv(fh_list, adjusted)
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat <- gsMat[, c(idx.up, idx.down)]
  
  # add annotation
  ann.col <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
  rownames(ann.col) <- colnames(gsMat)
  myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"))
  
  #gsMat <- getHeatmapMatrix(fh_list)	
  #colnames(gsMat) <- names(fh_list)	
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  gsMat <- -log10(gsMat)	
  gsMat[gsMat > 5] <- 5 # set the limit to 5
  
  paletteLength <- 10
  myColor <- magma(paletteLength)
  myBreaks <- c(seq(0, max(gsMat), length.out=ceiling(paletteLength)))
  
  #paletteLength <- 25
  #myMax <- ceiling(max(abs(gsMat)))
  
  #myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  #myBreaks <- myBreaks[myBreaks != 0]
  #myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
           cluster_cols = FALSE, cluster_rows = doClust,
           annotation_col = ann.col, annotation_colors = myColor.ann,
           cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8
           #width = 9, height = 9
  )
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





###################################################################
# FSGEA BASED ON FC

# Load FC
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/deseq2/"))
#limmaFiles <- list.files(pattern = "_DESeq2.xlsx")
limmaFiles <- c("OSM-CTL_DESeq2.xlsx")
names(limmaFiles) <- gsub("_DESeq2.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)
limmaList <- lapply(limmaList, function(i) i[order(i$entrez), ])

# Foldchange matrix
fcMat <- do.call(cbind, lapply(limmaList, function(i) i$log2FoldChange))
rownames(fcMat) <- limmaList[[1]]$entrez

# Rank genes based on z-score
rankedGenesList <- lapply(1:ncol(fcMat), function(i) {
  zs <- fcMat[,i]
  names(zs) <- rownames(fcMat)
  zs <- sort(zs)
})
names(rankedGenesList) <- colnames(fcMat)	


# intersect db and fcmat
dbList <- lapply(dbList, function(i){
  current <- lapply(i, intersect, y = rownames(fcMat))
  current.length <- sapply(current, length)
  return(current[current.length >= 5])
})


# FGSEA
fgseaDir <- file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/gsea/fgsea/")
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  plotDir <- fgseaDir   
  
  dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(plotDir, mydb.name))   
  
  #########################
  # Perform fgsea analysis
  
  fgseaResList <- lapply(rankedGenesList,
                         fgsea, pathways = mydb, minSize=15, maxSize=500, nperm=10000
  )
  fgseaResList <- lapply(fgseaResList, function(i) i[order(i$pval),])	
  
  # Add symbols
  fgseaResList <- lapply(fgseaResList, as.data.frame)
  
  fgseaResList <- lapply(fgseaResList, function(k){
    k.symbol <- lapply(k$leadingEdge, entrez2symbol)
    k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
    k.symbol <- lapply(k.symbol, toString)
    k$leadingEdge.symbol <- unlist(k.symbol)
    k
  })
  
  # Save fgsea output   
  lapply(1:length(fgseaResList), function(i)
    write.xlsx(fgseaResList[[i]],
               paste(names(fgseaResList)[i], "_", mydb.name, "_fgsea.xlsx", sep = ""),
               row.names = FALSE)
  )	
  
  # Plot significant pathways
  #gsSignif <- lapply(fgseaResList, function(i)
  #  i$pathway[1:5])
  #gsSignif <- sort(unique(unlist(gsSignif)))
  
  #gsNames <- gsSignif
  #for(gs in gsNames)
  #{
  #  gs.name <- gsub("\\/", "-", gs)
  #  
  #  p1 <- plotEnrichment(mydb[[gs]],
  #                       rankedGenesList[[1]]) + labs(title=names(rankedGenesList)[1])
  #  
  #  pdf(paste(gs.name, ".pdf", sep = ""), width = 5, height = 5)
  #  plot(p1)
  #  dev.off()                
  #  
  #}
  
  fgseaResList <- lapply(fgseaResList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j)
  })
  
  # Heatmaps
  
  # Divide UP and DOWN
  #fgseaResList.UP <- lapply(fgseaResList, function(j) j[j$NES > 0, ])
  #names(fgseaResList.UP) <- paste0(names(fgseaResList), ".UP")
  #fgseaResList.DOWN <- lapply(fgseaResList, function(j) j[j$NES < 0, ])
  #names(fgseaResList.DOWN) <- paste0(names(fgseaResList), ".DOWN")
  #fgseaResList <- c(fgseaResList.UP, fgseaResList.DOWN)
  #fgseaHeatmap(fgseaResList, paste0(mydb.name, "_fgsea_heatmap.pdf"), 10, adjusted = TRUE)
  
})



# PLOT RELEVENT GENE-SETS

# goBP
plotDir <- fgseaDir  
mydb <- dbList[["exhaustion_down"]]
mydb.name <- "exhaustion_down"

dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
setwd(file.path(plotDir, mydb.name))   

#gsNames <- c("GO_LYMPHOCYTE_ACTIVATION")
gsNames <- names(mydb)
for(gs in gsNames)
{
  gs.name <- gsub("\\/", "-", gs)
  
  p1 <- plotEnrichment(mydb[[gs]],
                       rankedGenesList[[1]]) + labs(title=names(rankedGenesList)[1])

  pdf(paste(gs.name, ".pdf", sep = ""), width = 5, height = 5)
  plot(p1)
  dev.off()                
  
}


#######################################
# OVERALL BARPLOT

# Load all databases for a given comparison
setwd(fgseaDir)
nb <- 50
comp <- c("OSM-CTL")
dbList.sub <- dbList[-7]# remove CGP



lapply(comp, function(i){
  fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_", j, "_fgsea.xlsx")))
  names(fhFiles) <- paste0(i, "_", names(fhFiles))
  fhList <- lapply(fhFiles, read.xlsx, sheet = 1)

  fhList.up <- lapply(fhList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j[j$NES > 0, ])
  }) 
  
  fhList.down <- lapply(fhList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j[j$NES < 0, ])
  }) 
  
  # UP
  fh.up <- do.call(rbind, fhList.up)
  fh.up <- fh.up[order(-fh.up$NES), ]
  
  ggmat.up <- fh.up[1:nb, c("pathway", "pval", "size", "NES")]
  ggmat.up$Score <- -log10(ggmat.up$"pval")
  ggmat.up$pathway <- substring(ggmat.up$pathway, 1, 60)
  ggmat.up$pathway <- factor(ggmat.up$pathway, levels = ggmat.up$pathway[order(ggmat.up$NES)])
  
  p <- ggplot(data=ggmat.up, aes(x=pathway, y=NES, fill = size))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("Normalized Enrichment Score")
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "UP_overall_barplot.pdf"), width = 8, height = 10)
  
  # DOWN
  fh.down <- do.call(rbind, fhList.down)
  fh.down <- fh.down[order(fh.down$NES), ]
  
  ggmat.down <- fh.down[1:nb, c("pathway", "pval", "size", "NES")]
  ggmat.down$Score <- -log10(ggmat.down$"pval")
  ggmat.down$pathway <- substring(ggmat.down$pathway, 1, 60)
  ggmat.down$pathway <- factor(ggmat.down$pathway, levels = rev(ggmat.down$pathway[order(ggmat.down$NES)]))
  
  p <- ggplot(data=ggmat.down, aes(x=pathway, y=NES, fill = size))
  p <- p +  geom_bar(stat="identity")
  p <- p + coord_flip()
  p <- p + scale_fill_viridis_c()
  p <- p + theme_bw()
  p <- p + xlab("") + ylab("Normalized Enrichment Score")
  
  ggsave(plot = p, filename = paste0(i, "_top", nb, "DOWN_overall_barplot.pdf"), width = 8, height = 10)
  
})







