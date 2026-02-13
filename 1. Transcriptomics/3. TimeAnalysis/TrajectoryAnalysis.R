# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(igraph)
library(tidyverse)
library(ggpubr)

# Capitalize first letter
firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/3. TimeAnalysis"))

# Load data:

# GO annotation
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))

# Gene annotation
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"geneAnnotation.RData"))

# Expression matrix
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))

# Sample info
load(paste0(homeDir,"/SampleInfo.RData"))
sampleInfo$Tissue[sampleInfo$Tissue == "Cell"] <- "iPSC"

# Forebrain development
GOid <- "GO:0030900"

# Retrieve genes
geneList <- GOgenes[[GOid]]

# filter gx matrix
gxMatrix_fil <- gxMatrix_norm[rownames(gxMatrix_norm) %in% geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% geneList],]
#gxMatrix_fil <- gxMatrix_norm

# Function for calculating Z score
calZ <- function(x, rowM = rowMeans(x), rowSD = apply(x,1,sd)){
    (x-rowM)/(rowSD)
}


# select Dorsal or ventral
tissue <- "Ventral"

# Select relevant samples
sampleInfo_IC <- sampleInfo[(sampleInfo$Tissue %in% c(tissue, "iPSC")) &
                                (sampleInfo$Group == "IC"),]
sampleInfo_RTT <- sampleInfo[(sampleInfo$Tissue %in% c(tissue, "iPSC")) &
                                 (sampleInfo$Group == "RTT"),]


gxMatrix_IC <- gxMatrix_fil[,sampleInfo_IC$SampleID]
gxMatrix_IC <- gxMatrix_IC[apply(gxMatrix_IC,1,sd) > 0,]
gxMatrix_IC_Z <- calZ(gxMatrix_IC)


gxMatrix_RTT <- gxMatrix_fil[rownames(gxMatrix_IC), sampleInfo_RTT$SampleID]
gxMatrix_RTT_Z <- calZ(gxMatrix_RTT,
                     rowM = rowMeans(gxMatrix_IC),
                     rowSD = apply(gxMatrix_IC,1,sd))

# Run PCA
pca <- prcomp(t(gxMatrix_IC_Z), 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = FALSE,
              scale = FALSE)


# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)

# Combine sample info in PCA scores
plotPCA <- as.data.frame(pca$x)

# Add sample information to dataframe
if (all(rownames(plotPCA) == sampleInfo_IC$SampleID)){
    plotPCA <- cbind(plotPCA, sampleInfo_IC)
} else{
    print("Samples in wrong order")
}

# Fit spline
spline_fit <- smooth.spline(plotPCA$PC1, plotPCA$PC2, spar = 0.5)

# Project RTT samples onto PCA plot
projectPCA <- as.data.frame(t(gxMatrix_RTT_Z) %*% pca$rotation)

# Add sample information to dataframe
if (all(rownames(projectPCA) == sampleInfo_RTT$SampleID)){
    projectPCA <- cbind(projectPCA, sampleInfo_RTT)
}

plotPCA <- rbind(plotPCA, projectPCA)


# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")

plotSpline <- data.frame(predict(spline_fit, seq(min(plotPCA$PC1),
                                                 max(plotPCA$PC1),0.01)))



p <- ggplot() +
    geom_point(data = plotPCA, 
               aes(x = PC1, y = PC2, color = Colour), size = 3) +
    geom_line(data = plotSpline,
              aes(x = x, y = y),
              color = "#084594") +
    xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
    ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
    labs(color = NULL) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 18,
                                    face = "bold"),
          panel.grid.major = element_line(color = "grey",
                                          linewidth = 0.1,
                                          linetype = 1),
          panel.grid.minor= element_line(color = "lightgrey",
                                         linewidth = 0.05,
                                         linetype = 1),          
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"))


p

ggsave(filename = paste0("PseudoTime_",tissue,".png"), p, width = 7, height = 4)






plotSpline$dist <- NA
plotSpline$dist[1] <- 0
for (i in 2:nrow(plotSpline)){
    plotSpline$dist[i] <- plotSpline$dist[i-1] + 
        sqrt(((plotSpline$x[i]-plotSpline$x[i-1])^2 + 
                  (plotSpline$y[i]-plotSpline$y[i-1])^2))
    
}


pseudoTime <- function(x, y, plotSpline){
    distance <- sqrt(((plotSpline$x-x)^2 + (plotSpline$y-y)^2))
    return(plotSpline$dist[which.min(distance)])
}

plotPCA$pseudoTime <- NA
for (i in 1:nrow(plotPCA)){
    plotPCA$pseudoTime[i] <- pseudoTime(x = plotPCA$PC1[i], 
                                        y = plotPCA$PC2[i], plotSpline)
}


p <- ggplot() +
    geom_point(data = plotPCA,
               aes(x = as.numeric(str_remove(Time, "D")), 
                   y = pseudoTime, color = Group)) +
    geom_line(data = plotPCA,
              aes(x = as.numeric(str_remove(Time, "D")), 
                  y = pseudoTime, color = Group,
                  group = paste0(Group, Replicate))) +
    scale_color_manual(values = c("#2171B5", "#CB181D")) +
    ylab("Pseudotime") +
    xlab("Time (days)") +
    labs(color = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 18,
                                    face = "bold"),
          panel.grid.major = element_line(color = "grey",
                                          linewidth = 0.1,
                                          linetype = 1),
          panel.grid.minor= element_line(color = "lightgrey",
                                         linewidth = 0.05,
                                         linetype = 1),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"))

ggsave(filename = paste0("TimevsPseudoTime_",tissue,".png"), p, width = 7, height = 4)


t.test(plotPCA$pseudoTime[(plotPCA$Time == "D0")&
                              (plotPCA$Group == "RTT")],
       plotPCA$pseudoTime[(plotPCA$Time == "D0")&
                              (plotPCA$Group == "IC")])

t.test(plotPCA$pseudoTime[(plotPCA$Time == "D13")&
                              (plotPCA$Group == "RTT")],
       plotPCA$pseudoTime[(plotPCA$Time == "D13")&
                              (plotPCA$Group == "IC")])

t.test(plotPCA$pseudoTime[(plotPCA$Time == "D40")&
                              (plotPCA$Group == "RTT")],
       plotPCA$pseudoTime[(plotPCA$Time == "D40")&
                              (plotPCA$Group == "IC")])


t.test(plotPCA$pseudoTime[(plotPCA$Time == "D75")&
                              (plotPCA$Group == "RTT")],
       plotPCA$pseudoTime[(plotPCA$Time == "D75")&
                              (plotPCA$Group == "IC")])


