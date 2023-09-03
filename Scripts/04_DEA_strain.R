
########################################################################
#clean.workplace
########################################################################

rm(list = ls())

########################################################################
#library
########################################################################
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
library(data.table)
####################################################################
######data processing 
####################################################################

`%ni%` <- Negate(`%in%`)

#read metadata
load(paste0("./Data/Proximal_colon_gene_expression.RData"), verbose = T)

data$id      <- paste0(data$strain, "_", data$diet)

# check if there are strains only have one diet, and remove them

strain_cd    <- unique(data$strain[data$diet == "CD"])
strain_hfd   <- unique(data$strain[data$diet == "HFD"])
    
strains      <- intersect(strain_cd, strain_hfd)

#Remove strains with only one type of diet
genedata_res       <- data[data$strain %in% strains, ]
genedata_res$split <- paste0(genedata_res$strain, "_", genedata_res$gene_symbol)
lt                 <- split(genedata_res, genedata_res$split)

FC                 <- do.call(rbind, mclapply(lt, function(x){ #x = lt[[1]]
  
     tmp           <- log2(x$value[x$diet == "HFD"]) - log2(x$value[x$diet == "CD"])
     res           <- data.frame(strain = unique(x$strain), gene_symbol = unique(x$gene_symbol), FC_value = tmp)
     res
     
}, mc.cores = 20))

#genesets

Genesets <-list( 
  Reactome = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
  Hallmarks = msigdbr(species = "Mus musculus", category = "H"),
  GOBP = msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
  GOCC = msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"),
  Custom_genesets = read.csv("./Data/Custom_genesets.csv")
  #GOMF = msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF")
)

load(paste0("./Result/Differential_expression.RData"))

FC           <- FC[FC$gene_symbol %in% tt$`HFD - CD`$gene_id, ]

tt           <- split(FC, FC$strain)

save(tt, file="./Result/Differential_expression_strain.RData")

#GSEA analysis 
allGSEA <- lapply(tt, function(x){ # x= tt[[1]]
  print(unique(x$strain))
  geneList <- x$FC_value
  names(geneList) <- x$gene_symbol
  geneList <- sort(geneList, decreasing = T)
  lapply(Genesets, function(z){ #z = Genesets[[1]]
    print(unique(z$gs_subcat))
    z <- z[, c("gs_name", "gene_symbol")]
    colnames(z) <- c("ont", "gene")
    set.seed(30072021)
    GSEA_result <- GSEA(geneList, TERM2GENE = z, nPerm = 100000, minGSSize = 0, maxGSSize = 5000, pvalueCutoff = 1, verbose = F, seed = F)
    return(GSEA_result)
  })
})
names(allGSEA) <- names(tt)

#merge result
allgsea_result <- lapply(allGSEA,  function(x){ # x= allGSEA[[1]]
  result_1 <- mapply(z = x, y=names(x),function(z, y){ # z= x[[1]]
    result = z@result
    result$category = y
    return(result)
  }, SIMPLIFY = FALSE)
  result_2 <- do.call(rbind, result_1)
})
names(allgsea_result) <- names(allGSEA) 

#Figure 1—figure supplement 1D

plot_lt <- do.call(rbind, mapply(comparison = allgsea_result, name = names(allgsea_result), function(comparison, name){ # comparison <- allgsea_result[[1]]
  
  result   <- comparison
  result$core_enrichmentNgenes <- sapply(strsplit(result$core_enrichment, split = "/", fixed = T), length)
  result$core_enrichmentgenes <- unlist(sapply(strsplit(result$core_enrichment, split = "/", fixed = T), function(x){paste(sort(x), collapse = " ")}))
  result$core_enrichmentgenes <- gsub('(.{1,60})(\\s|$)', '\\1\n', result$core_enrichmentgenes)
  result$signficant <- result$qvalues <= 0.05 & abs(result$NES) >1
  result$label <- gsub("\\(.*\\)", "", gsub("_", " ", result$Description))
  result$label <- tolower(result$label)
  result$label <- gsub("reactome ", "", gsub("^kegg.*? ", "", gsub("^hallmark.*? ", "", gsub("^go.*? ", "", result$label))))
  result$label <- capitalize(result$label)
  result$gene_ratio <- result$core_enrichmentNgenes / result$setSize
  result$comparison <- name
  
  result <- result[order(result$ID), ]
  
  return(result)
}, SIMPLIFY = FALSE))

plot_lt$stars   <- cut(plot_lt$p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
plot_lt$tissue  <- tissue
plot_lt


Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))


plot_lt_colon   <- plot_lt


#####Sort by significance
plot_lt_colon         <- plot_lt_colon[order(plot_lt_colon$NES, decreasing = T),]
# Hierarchical clustering on phenotypes1
data_dcast              <- reshape2::dcast(plot_lt_colon, comparison~label, value.var="NES")
data_dcast[is.na(data_dcast)] <- 0 # For the purpose of clustering, we consider "NA" equivalent to no correlation
rownames(data_dcast) <- data_dcast$comparison
hc <- hclust(dist(t(data_dcast[,-1])))
rowInd <- hclust(dist(data_dcast[,-1]))$order
colInd <- hclust(dist(t(data_dcast[,-1])))$order


#plot_lt$label   <- factor(plot_lt$label, levels = colnames(data_dcast)[-1][colInd])
plot_lt$comparison <- factor(plot_lt$comparison, levels = data_dcast$comparison[rowInd])

g <- ggplot(plot_lt, aes(x = comparison, y = label)) + 
  geom_point(aes(size = -log10(p.adjust), col = NES )) +
  scale_size(range = c(2,10)) +
  geom_text(data = subset(plot_lt, signficant = "TRUE"), aes(label=stars), color="black", size=5, vjust = 0.8, hjust = 0.5)  +
  #facet_grid(rows = vars(tissue) , space="free", scales="free") +
  scale_colour_gradientn(colours= Heatmap_palette , limits=c(-2.5,2.5) ,na.value="gray87", oob = squish) +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle= 90, hjust=1, vjust=0.5,size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 13)
    #panel.border = element_blank()
  ) 

