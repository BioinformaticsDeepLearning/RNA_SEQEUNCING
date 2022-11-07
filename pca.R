library(ggplot2)
library(labdsv)
library(limma)
library(openxlsx)
library(plyr)
library(ggrepel)
library(ggforce)
library(org.Mm.eg.db)
library(raster)

source("/Users/gandrieux/bitbucket/work/tools/tools.r")

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


ggPCAplot <- function(pcaObj, groupColor, groupShape = NULL, text = FALSE, polygon = FALSE, loadings = FALSE, output)
{
  pcaSum <- summary(pcaObj)
  pcaMatrix <- data.frame(sample = rownames(pcaObj$score),
                          pcaObj$score, GroupColor = groupColor)	
  if(!is.null(groupShape)) pcaMatrix$GroupShape = groupShape
  
  xlim <- ceiling(max(abs(pcaMatrix$PC1)) + 0.15 * max(abs(pcaMatrix$PC1)))
  ylim <- ceiling(max(abs(pcaMatrix$PC2)) + 0.15 * max(abs(pcaMatrix$PC2)))
  p <- ggplot(pcaMatrix)
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
                 axis.line.y = element_line(color="black", size = 0.5))		
  if(!is.null(groupShape)) p <- p + geom_point(aes(PC1, PC2, col = GroupColor, shape = GroupShape), size = 4, alpha = 0.5)
  else p <- p + geom_point(aes(PC1, PC2, col = GroupColor), size = 4, alpha = 0.5)
  #p <- p + scale_colour_gradientn(colours = terrain.colors(10))
  p <- p + scale_x_continuous(limits = c(-xlim, xlim))
  p <- p + scale_y_continuous(limits = c(-ylim, ylim))
  p <- p + theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
  p <- p + xlab(paste("PC1 (",round(100*pcaSum[2,1],2),"%)",sep="")) +  ylab(paste("PC2 (",round(100*pcaSum[2,2],2),"%)",sep=""))
  #p <- p + geom_text(aes(PC1, PC2, label = sample, colour = GroupColor), vjust = -0.5, size = 4)
  
  if(text){
    p <- p + geom_text_repel(aes(PC1, PC2, label=sample), cex = 3,
                             segment.size = 0.1, col = "black", min.segment.length = unit(0, 'lines'))
  }
  
  if(polygon){
    #getting the convex hull of each unique point set
    df <- pcaMatrix
    find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
    hulls <- ddply(df, "GroupColor", find_hull)
    p <- p + geom_polygon(data = hulls, alpha = 0.5, aes(PC1, PC2, fill = GroupColor))
  }
  
  if(loadings){
    # Set new limits
    abslim <- max(xlim, ylim)
    p <- p + scale_x_continuous(limits = c(-abslim, abslim))
    p <- p + scale_y_continuous(limits = c(-abslim, abslim))
    
    # Select top loadings
    ldg <- pcaObj$loadings
    aldg <- abs(ldg)
    aldg <- sweep(aldg, 2, colSums(aldg), "/")
    
    if(nrow(aldg) > 10)
    {
      pc1.top <- rownames(aldg[order(aldg[,1], decreasing = TRUE),])[1:5]
      pc2.top <- rownames(aldg[order(aldg[,2], decreasing = TRUE),])[1:5] 
      ldg <- ldg[unique(c(pc1.top, pc2.top)), ]			
    }
    
    # change coordinates	
    circleLimit <- abslim - (0.20 * abslim)
    ldg.max <- max(abs(ldg[,1:2]))
    ldg[,1] <- (circleLimit * ldg[,1]) / ldg.max
    ldg[,2] <- (circleLimit * ldg[,2]) / ldg.max
    
    ldg.gg <- data.frame(x=rep(0, nrow(ldg)), y=rep(0, nrow(ldg)),
                         vx=ldg[,1], vy=ldg[,2], label = rownames(ldg))	
    
    p <- p + geom_segment(data=ldg.gg, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                          arrow=arrow(unit(0.05, "inches"), type = "closed", angle = 30), size=0.5, color="black", alpha = 0.5)
    p <- p + geom_text_repel(data=ldg.gg, aes(vx, vy, label = label), size = 3)	
    
    # draw circle
    circleRadius <- apply(ldg[,1:2], 1, function(i)
      pointDistance(c(0, 0), i, lonlat = FALSE)
    )
    circleRadius <- max(circleRadius)	
    circles <- data.frame(
      x0 = 0,
      y0 =  0,
      r = circleRadius
    )
    p <- p + geom_circle(aes(x0=x0, y0=y0, r=r), data=circles) + coord_fixed()
    p <- p + annotate("text", x = circleRadius, y = circleRadius,
                      label = paste("r == ", ldg.max),
                      parse = TRUE)
  }
  
  pdf(output)
  plot(p)
  dev.off()	
}


ggLoadings <- function(pcaObj, output)
{
  pcaSum <- summary(pcaObj)
  
  ldg <- pcaObj$loadings
  aldg <- abs(ldg)
  aldg <- sweep(aldg, 2, colSums(aldg), "/")
  
  if(nrow(aldg) > 10)
  {
    pc1.top <- rownames(aldg[order(aldg[,1], decreasing = TRUE),])[1:5]
    pc2.top <- rownames(aldg[order(aldg[,2], decreasing = TRUE),])[1:5] 
    ldg <- ldg[unique(c(pc1.top, pc2.top)), ]			
  }
  
  ldg.gg <- data.frame(x=rep(0, nrow(ldg)), y=rep(0, nrow(ldg)),
                       vx=ldg[,1], vy=ldg[,2], label = rownames(ldg))
  
  xlim <- max(abs(ldg.gg$vx)) + 0.15 * max(abs(ldg.gg$vx))
  ylim <- max(abs(ldg.gg$vy)) + 0.15 * max(abs(ldg.gg$vy))
  abslim <- max(xlim, ylim)
  
  
  circleRadius <- apply(ldg[,1:2], 1, function(i)
    pointDistance(c(0, 0), i, lonlat = FALSE)
  )
  circleRadius <- max(circleRadius)	
  
  circles <- data.frame(
    x0 = 0,
    y0 =  0,
    r = circleRadius
  )
  
  p <- ggplot(ldg.gg)
  p <- p + geom_segment(mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                        arrow=arrow(unit(0.05, "inches"), type = "closed", angle = 30), size=0.5, color="black")
  p <- p + geom_text_repel(data=ldg.gg, aes(vx, vy, label = label), size = 3,
                           segment.color = "grey")	
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.line.x = element_line(color="black", size = 0.5),
                 axis.line.y = element_line(color="black", size = 0.5))		
  p <- p + scale_x_continuous(limits = c(-abslim, abslim))
  p <- p + scale_y_continuous(limits = c(-abslim, abslim))
  p <- p + xlab(paste("PC1 (",round(100*pcaSum[2,1],2),"%)",sep="")) +  ylab(paste("PC2 (",round(100*pcaSum[2,2],2),"%)",sep=""))
  
  p <- p + geom_circle(aes(x0=x0, y0=y0, r=r), data=circles) + coord_fixed()
  
  pdf(output)
  plot(p)
  dev.off()
}




############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


##########################################################################################

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/count/"))
expr <- read.delim("rlogTransformation_DESeq2.txt")# DESeq2

# ANNOTATIONS
setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/doc"))
ann.sample <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann.sample <- ann.sample[match(colnames(expr), ann.sample$SAMPLE), ]

mat <- expr
group.color <- ann.sample$GROUP
group.shape <- NULL

# PCA
p5 <- pca(t(mat),dim=4,cor=F)

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/pca/"))
ggPCAplot(p5,
          groupColor = group.color,
          groupShape = group.shape,
          text = TRUE,
          polygon = FALSE,
          loadings = FALSE,
          output = "all_RNA_category_PCA.pdf")

pdf("all_RNA_category_hclust.pdf")
plot(hclust(dist(t(mat))))
dev.off()

#########
# MAD 10%
mat <- expr
nbMAD <- round(nrow(mat) / 10)
mad <- apply(mat, 1, mad)
mat <- mat[order(mad, decreasing = TRUE)[1:nbMAD],]

# PCA
p5 <- pca(t(mat),dim=4,cor=F)

setwd(file.path("/Volumes/Home/Geoffroy/Duyster/RNA_010719/pca/"))
ggPCAplot(p5,
          groupColor = group.color,
          groupShape = group.shape,
          text = TRUE,
          polygon = FALSE,
          loadings = FALSE,
          output = "all_MAD10PC_RNA_category_PCA.pdf")

pdf("all_MAD10PC_RNA_category_hclust.pdf")
plot(hclust(dist(t(mat))))
dev.off()

