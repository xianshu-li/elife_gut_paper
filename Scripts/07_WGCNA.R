########################################################################
#clean.workplace
########################################################################

rm(list = ls())
# clear all plots
dev.off()
# clear console
cat("\014")

########################################################################
#library
########################################################################
library(WGCNA)
library(edgeR)
library(limma)
library(FactoMineR)
library(plotly)
library(cowplot)
library(RColorBrewer)
library(pals)
library(preprocessCore)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(scales)
library(here)
library(dplyr)
library(tidyr)
library(parallel)
library(reshape2)
library(stringr)
library(Biobase)
library(ggfortify)
`%ni%` <- Negate(`%in%`)

clusterWGCNA <- function(dataset, powers, analysis){ #dataset = datExpr0 analysis = diet
  print("running pickSoftThreshold")
  set.seed(20000001)
  sft <- pickSoftThreshold(t(dataset), 
                           powerVector = powers, 
                           verbose = 5, 
                           networkType = "signed hybrid", 
                           blockSize = 25000, 
                           corFnc = "bicor",corOptions = list(quick = 0.5, use = 'pairwise.complete.obs'),
                           moreNetworkConcepts = T)
  
  pdf(paste0("./Data/WGCNA/", analysis, "/09_module_generation_pickSoftThreshold.pdf"), width = 10, height = 10)
  cex1 = 0.9;
  
  par(mfrow = c(1,2));
  plot(sft$fitIndices[,1], sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  dev.off()    
  
  fit.sequence <- sft$fitIndices[,"SFT.R.sq"]
  names(fit.sequence) <- sft$fitIndices[,"Power"]
  fit.sequence[1] <- 0 #never take the first power
  bestpower <-  as.numeric(ifelse(max(fit.sequence)> 0.85, names(fit.sequence)[min(which(fit.sequence > 0.85))], names(fit.sequence)[which.max(fit.sequence)]))
  minmodsize <- 30
  
  set.seed(20000001)
  net <- blockwiseModules(t(dataset), power = bestpower, networkType = "signed hybrid",
                          TOMType = "signed", minModuleSize = minmodsize,
                          reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                          numericLabels = TRUE, pamStage = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, minKMEtoStay = 0.2, minCoreKME = 0.3, minCoreKMESize = 5,
                          saveTOMFileBase = paste0("./Data/WGCNA/", analysis), useBranchEigennodeDissim = FALSE,
                          verbose = 1, maxBlockSize = 25000, corType = "bicor", quickCor = 1,
                          nThreads = 0)
  
  mergedColors = labels2colors(net$colors)
  
  
  rownames(net$MEs) <- colnames(dataset)
  ME <- net$MEs
  ME.colors <- labels2colors(as.numeric(gsub("ME", "", colnames(ME))))
  METree <- hclust(dist(t(ME)), method = "average");
  
  net$METree <- METree
  net$ME.colors <- ME.colors
  net$sft <- sft
  clusteringResult <- data.frame(cluster = as.numeric(net$colors), gene = rownames(dataset))
  net$clusteringResult <- clusteringResult
  net$analysis <- analysis
  net$dataset <- dataset
  net$powers  <- powers
  net$bestpower <- bestpower
  saveRDS(net,paste0("./Data/WGCNA/", analysis, "_net.RDS"))
}



Diet_expression <- c("CD", "HFD")

load("./Data/Proximal_colon_gene_expression.RData", verbose = T)
  
#load annotated gene expression data
data$id         <- paste0(data$strain, "_", data$diet)
tmp             <- as.data.frame(data[, c("id", "gene_symbol", "value")])
tmp             <- reshape2::dcast(tmp, gene_symbol~id, value.var = "value")
tmp             <- tmp[!is.na(tmp$gene_symbol),]
rownames(tmp)   <- tmp$gene_symbol
tmp_new         <- 	tmp[, !grepl("C57BL/6J|DBA/2J", colnames(tmp))]
exp             <- edgeR::DGEList(tmp_new[, -1])
exp$samples$strain <- gsub("_cd|_hfd|_CD|_HFD", "", rownames(exp$samples))
exp$samples$diet   <- gsub(".*_", "", rownames(exp$samples))
exp$samples$strain_diet <- rownames(exp$samples)
  
#normalization
exp <- calcNormFactors(exp, method = "TMM")
load(paste0("./Result/Differential_expression.RData"))
tmp <- tt$`HFD - CD`
exp$counts  <- exp$counts[tmp$gene_id, ]
  
exp$design <- model.matrix(~0 + diet + strain, data = exp$samples)
colnames(exp$design) <- make.names(colnames(exp$design))
exp$voom   <- voom(exp, design = exp$design, plot = T)
  
lapply(Diet_expression, function(diet){ #diet = Diet_expression[2]
    
print(diet)
file_out     <- paste0("./Data/WGCNA/", diet)
    
    if(!dir.exists(file_out)){
      dir.create(file_out, recursive = T)
    }
    
    data.Expr <- exp$counts[, grepl(diet, colnames(exp$counts))]

    data.Expr  <- as.data.frame(data.Expr)
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    clusterWGCNA(dataset = data.Expr, powers = powers, analysis = diet)
    
    return(NULL)
    
})

