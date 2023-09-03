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
library(plyr)
library(forcats)
library(scales)
library(htmlwidgets)
library(Hmisc)
library(ggrepel)
library(dendextend)
library(circlize)
library(ggdendro)

`%ni%` <- Negate(`%in%`)

Zscore_for      <- function(x){
  z <- (x - mean(x)) / sd(x)
  return(z)
}

genesets_mouse     <- read.table("./Data/custom_disease_signature_mouse.txt", header = T, sep = "\t", stringsAsFactors = F)
genesets_human     <- read.table("./Data/custome_disease_signature_human.txt", header = T, sep = "\t", stringsAsFactors = F)

genesets           <- rbind(genesets_mouse, genesets_human)

load("./Result/Differential_expression_strain.RData")
#GSEA analysis
allGSEA <- lapply(tt, function(x){ # x= tt[[1]]
  geneList <- x$FC_value
  names(geneList) <- x$gene_symbol
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

plot_lt <- do.call(rbind, mapply(comparison = allgsea_result, name = names(allgsea_result), function(comparison, name){ # comparison <- allgsea_result[[1]]
  comparison$ID     <- rownames(comparison)
  result <- comparison
  result$core_enrichmentNgenes <- sapply(strsplit(result$core_enrichment, split = "/", fixed = T), length)
  result$core_enrichmentgenes <- unlist(sapply(strsplit(result$core_enrichment, split = "/", fixed = T), function(x){paste(sort(x), collapse = " ")}))
  result$core_enrichmentgenes <- gsub('(.{1,60})(\\s|$)', '\\1\n', result$core_enrichmentgenes)
  result$signficant <- result$p.adjust <= 0.05 & abs(result$NES) >1
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
plot_lt$species <- ifelse(grepl("day", plot_lt$ID), "mouse", "human")
save(plot_lt, file = "./Data/GSEA_IBD_strain.RData")
#load("./Data/GSEA_IBD_strain.RData", verbose = T)

IBD_colon_mouse <- plot_lt[plot_lt$species == "mouse", ]
IBD_colon_mouse <- IBD_colon_mouse[!grepl("day2|day4", IBD_colon_mouse$ID), ]
IBD_colon_mouse$Category<- ifelse(grepl("up", IBD_colon_mouse$ID), "Up","Down")

Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))

#####Sort by significance
plot_lt_colon_mouse         <- IBD_colon_mouse[order(IBD_colon_mouse$NES, decreasing = T),]

#only use mouse data to do the cluster again
which.col    <- c("ID", "NES", "p.adjust","comparison")

IBD_colon_mouse_cluster          <- IBD_colon_mouse[, which.col]
IBD_colon_mouse_cluster$sig      <- IBD_colon_mouse_cluster$NES * -log10(IBD_colon_mouse_cluster$p.adjust)

IBD_colon_mouse_cluster          <- reshape2::dcast(IBD_colon_mouse_cluster, ID ~ comparison, value.var = "sig")
rownames(IBD_colon_mouse_cluster) <- IBD_colon_mouse_cluster$ID

data_dcast              <- IBD_colon_mouse_cluster
data_dcast[is.na(data_dcast)] <- 0 # For the purpose of clustering, we consider "NA" equivalent to no correlation
rownames(data_dcast) <- data_dcast$ID
hc <- hclust(dist(t(data_dcast[,-1]), method = "manhattan"))
rowInd <- hclust(dist(data_dcast[,-1], method = "manhattan"))$order
colInd <- hclust(dist(t(data_dcast[,-1]), method = "manhattan"))$order

tree <- hclust(d = dist(x = as.matrix(t(IBD_colon_mouse_cluster[, colnames(IBD_colon_mouse_cluster) %ni% "ID"])), method = "manhattan"))
cl_members <- cutree(tree = tree, k = 3)


# Hierarchical clustering dendrogram
hc <- as.dendrogram(tree) 

# Colors
hc <- hc %>%
  color_branches(k = 3) %>%
  color_labels(k = 3)

# Circular dendrogram
circlize_dendrogram(hc,
                    labels_track_height = NA,
                    dend_track_height = 0.5)

clu <- ggdendrogram(hc, theme_dendro = F) +
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
    legend.position = "none",
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 13)
    #panel.border = element_blank()
  ) 

clu <- clu + scale_x_discrete(expand = c(0.01,0), breaks = c(colnames(data_dcast)[-1][colInd]))



strain        <- as.data.frame(cl_members)
strain$strain <- rownames(strain)
strain$direction <- strain$cl_members
strain$direction[strain$cl_members %in% c("1")] <- "Down"
strain$direction[strain$cl_members == "3"] <- "Up"
strain$direction[strain$cl_members %in% c("2")] <- "Other"

write.table(strain, "./Data/cluster_IBD_colon_mouse_new.txt", row.names = F, quote = F, sep = "\t")

plot_tmp            <- IBD_colon_mouse
plot_tmp$comparison <- factor(plot_tmp$comparison, levels = colnames(data_dcast)[-1][colInd])
plot_tmp$label   <- factor(plot_tmp$label, levels = rev(c(unique(plot_tmp$label[grepl(" up", plot_tmp$label)]), unique(plot_tmp$label[grepl("down", plot_tmp$label)]))))
g1 <- ggplot(plot_tmp, aes(x = comparison, y = label)) + 
  geom_point(aes(size = -log10(p.adjust), col = NES )) +
  scale_size(range = c(2,10)) +
  geom_text(data = subset(plot_tmp, signficant = "TRUE"), aes(label=stars), color="black", size=5, vjust = 0.8, hjust = 0.5)  +
  #facet_grid(rows = vars(Category) , space="free", scales="free") +
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
    legend.position = "bottom",
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 13)
    #panel.border = element_blank()
  ) 

IBD_colon_human          <- plot_lt[plot_lt$species == "human", ]

IBD_colon_human$Category <- ifelse(grepl("Up", IBD_colon_human$ID), "Up","Down")

plot_tmp                 <- IBD_colon_human
Heatmap_palette          <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))


plot_tmp$comparison      <- factor(plot_tmp$comparison, levels = colnames(data_dcast)[-1][colInd])
plot_tmp$label           <- factor(plot_tmp$label, levels = rev(c(unique(plot_tmp$label[grepl("\\.up", plot_tmp$label)]), unique(plot_tmp$label[grepl("\\.down", plot_tmp$label)]))))
g <- ggplot(plot_tmp, aes(x = comparison, y = label)) + 
  geom_point(aes(size = -log10(p.adjust), col = NES )) +
  scale_size(range = c(2,10)) +
  geom_text(data = subset(plot_tmp, signficant = "TRUE"), aes(label=stars), color="black", size=5, vjust = 0.8, hjust = 0.5)  +
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
    legend.position = "bottom",
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 13)
  ) 

p   <- plot_grid(clu, g1, g, nrow = 3, ncol = 1, align = "hv", rel_heights = c(0.8, 1, 1.5))

