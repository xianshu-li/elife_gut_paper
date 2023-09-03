########################################################################
#clean.workplace
########################################################################

rm(list = ls())

########################################################################
#library
########################################################################

# install.packages("dendextend")
# install.packages("circlize")

library(msigdbr)
library(cowplot)
library(plotly)
library(data.table)
library(RColorBrewer)
library(clusterProfiler)
library(parallel)
library(enrichplot)
library(ggplot2)
library(MASS)
library(plyr) # For reformatting data
library(forcats) # For the "fct_reorder" function wich reorders factors on a numeric variable
library(scales) # For "oob =squish" Makes out-of-bounds values not considered as NA in the color scale
library(htmlwidgets)
library(Hmisc)
library(ggrepel)
library(dendextend)
library(circlize)
`%ni%` <- Negate(`%in%`)

##########################Figure 1D##############################

genesets_mouse     <- read.table("./Data/custom_disease_signature_mouse.txt", header = T, sep = "\t", stringsAsFactors = F)
genesets_human     <- read.table("./Data/custome_disease_signature_human.txt", header = T, sep = "\t", stringsAsFactors = F)

genesets           <- rbind(genesets_mouse, genesets_human)

load(paste0("./Result/Differential_expression.RData"), verbose = T)
  #GSEA analysis
  allGSEA <- lapply(tt, function(x){ # x= tt[[1]]
    geneList <- x$logFC
    names(geneList) <- x$gene_id
    geneList <- sort(geneList, decreasing = T)
    z <- genesets[, c("condition", "GeneName")]
    colnames(z) <- c("ont", "gene")
    set.seed(30072021)
    GSEA_result <- GSEA(geneList, TERM2GENE = z, nPerm = 100000, minGSSize = 0, maxGSSize = 5000, pvalueCutoff = 1, verbose = F, seed = F)
    return(GSEA_result)
  })
  names(allGSEA) <- names(tt)
  
  #merge result
  allgsea_result <- lapply(allGSEA,  function(x){ # x= allGSEA[[1]]
    result = x@result
    return(result)
  })
  names(allgsea_result) <- names(allGSEA) 

  result <- allgsea_result$`HFD - CD`


Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))
result$stars   <- cut(result$p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
result$direction <- ifelse(grepl("up", result$Description), "Up", "Down")
result$disease   <- gsub("_down|_up", "", result$ID)
result           <- result[!grepl("T2D", result$ID) & !grepl("Ileum", result$tissue), ]

result$direction <- factor(result$direction, levels = c("Up", "Down"))

g <- ggplot(result, aes(x = direction, y = disease)) + 
  geom_point(aes(size = -log10(p.adjust), col = NES )) +
  scale_size(range = c(3,8)) +
  geom_text(data = subset(result, signficant = "TRUE"), aes(label=stars), color="black", size=5, vjust = 0.8, hjust = 0.5)  +
  scale_colour_gradientn(colours= Heatmap_palette , limits=c(-2.5,2.5) ,na.value="gray87", oob = squish) +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 13)

  ) 
