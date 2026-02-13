# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(stringr)
library(edgeR)
library(patchwork)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"                                      # Raw count data
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/",
            "geneAnnotation.RData"))                                       
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

###############################################################################

# Finding differentially expressed X-linked genes

###############################################################################

# Threshold for significance
pvalue_thres <- 0.05
logFC_thres <- 1
p_type <- "PValue"

# Make empty lists
up_genes <- list()
down_genes <- list()
both_genes <- list()
all_genes <- list()

# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

up_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[1]] <- rownames(output)

# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

up_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[2]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_iPSC.RData"))
output <- top.table_iPSC

up_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[3]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

up_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[4]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

up_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[5]] <- rownames(output)

# Load statistics
accessID <- "GSE117511"
load(file = paste0(accessID,"/top.table.RData"))
output <- top.table

up_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[6]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_CC.RData"))
output <- top.table_CC

up_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[7]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_TC.RData"))
output <- top.table_TC

up_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[8]] <- rownames(output)


validationGenes <- table(unlist(both_genes))
validationGenes <- names(validationGenes)[validationGenes>=4]
intersect(validationGenes, geneAnnotation$GeneName[geneAnnotation$chromosome_name == "X"])

###############################################################################

# Make plot

###############################################################################

# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron
PLP2 <- output["PLP2",]


# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])


# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_iPSC.RData"))
output <- top.table_iPSC
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])

# Load statistics
accessID <- "GSE117511"
load(file = paste0(accessID,"/top.table.RData"))
output <- top.table
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_CC.RData"))
output <- top.table_CC
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_TC.RData"))
output <- top.table_TC
PLP2 <- rbind.data.frame(PLP2,
                         output["PLP2",])


PLP2$Name <- factor(c("Neuron", "NPC", 
               "iPSC", "NPC", "Neuron",
               "Interneuron",
               "Cingulate\ncortex",
               "Temporal\ncortex"),
               levels = c("iPSC", "NPC", "Neuron",
                          "Interneuron",
                          "Cingulate\ncortex",
                          "Temporal\ncortex"))

PLP2$Dataset <- factor(c(rep("GSE123753",2),
                  rep("GSE107399",3),
                  rep("GSE117511",1),
                  rep("GSE128380",2)),
                  levels = c("GSE107399", "GSE123753",
                             "GSE117511", "GSE128380"))

p <- ggplot() +
    geom_bar(data = PLP2,
             aes(x = Name, y = -log10(PValue), fill = logFC), 
             color = "black", stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    facet_grid(cols = vars(Dataset), scale = "free", space = "free") +
    scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                         limit = c(-1.5,1.5), oob = scales::squish) +
    xlab(NULL) +
    ylab(expression(-log[10]~"P value")) +
    labs(fill = expression(log[2]~"FC")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(p, file = "PLP2_validation.png", 
       width = 7, height = 4)

