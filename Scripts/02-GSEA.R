
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
`%ni%` <- Negate(`%in%`)

##########################Read Genesets from msigdbr R package##############################

Genesets <- list(
  Reactome = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
  Hallmarks = msigdbr(species = "Mus musculus", category = "H"),
  GOBP = msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
  GOCC = msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"),
  Custom_genesets = read.csv("./Data/Custom_genesets.csv")
)


##########################GSEA analysis ##############################

load(paste0("./Result/Differential_expression.RData"))
tmp <- tt$`HFD - CD`

#GSEA analysis
allGSEA <- lapply(tt, function(x){ # x= tt[[1]]
  print(unique(x$contrast))
  geneList <- x$logFC
  names(geneList) <- x$gene_id
  geneList <- sort(geneList, decreasing = T)
  lapply(Genesets, function(z){ #z = Genesets[[1]]
    print(unique(z$gs_subcat))
    z <- z[, c("gs_name", "gene_symbol")]
    colnames(z) <- c("ont", "gene")
    set.seed(30072021)
    GSEA_result <- GSEA(geneList, TERM2GENE = z, nPerm = 100000, minGSSize = 30, maxGSSize = 5000, pvalueCutoff = 1, verbose = F, seed = F)
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

#save GSEA data to a Rdata file
save(allGSEA, allgsea_result, file =paste0("./Result/GSEA.RData"))

#########################Plot for GSEA (Figure 1C)##############################

load(paste0("./Result/GSEA.RData"), verbose = T)

colon_result       <- allgsea_result$`HFD - CD`

genesets <- openxlsx::read.xlsx("./Data/Selected_genesets_Colon.xlsx")
result   <- merge(result, genesets[, c("ID", "facet_category", "Label")], by = "ID")
result$core_enrichmentNgenes <- sapply(strsplit(result$core_enrichment, split = "/", fixed = T), length)
result$core_enrichmentgenes <- unlist(sapply(strsplit(result$core_enrichment, split = "/", fixed = T), function(x){paste(sort(x), collapse = " ")}))
result$core_enrichmentgenes <- gsub('(.{1,60})(\\s|$)', '\\1\n', result$core_enrichmentgenes)
result$signficant <- result$p.adjust <= 0.05 & abs(result$NES) >1
result$gene_ratio <- result$core_enrichmentNgenes / result$setSize
result$comparison <- "HFD vs CD"
result            <- result[order(result$ID), ]

Heatmap_palette <- c(rev(brewer.pal(7,"Blues")),"white",brewer.pal(7,"Reds"))

result$category_function <- factor(result$facet_category, level = c("Cell cycle", "Inflammation", "Mitochondria","Stress","Keratins"))

result$stars <- cut(result$p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

orders <- unique(result[order(result$NES, decreasing = T), "Label"])

result$Labels <- factor(result$Label, levels = rev(orders), ordered = T)

g <- ggplot(result, aes(x = NES, y = Labels)) + 
  geom_point(aes(size = -log10(p.adjust), col = NES )) +
  scale_size(range = c(3,8)) +
  geom_text(data = subset(result, signficant = "TRUE"), aes(label=stars), color="black", size=5, vjust = 0.7, hjust = 0.5)  +
  facet_grid(rows = vars(category_function) , space="free", scales="free") +
  scale_colour_gradientn(colours= Heatmap_palette , limits=c(-3,3) ,na.value="gray87", oob = squish) +
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
    strip.text = element_text(size = 13),
    legend.position = "bottom"
    #panel.border = element_blank()
  ) 

