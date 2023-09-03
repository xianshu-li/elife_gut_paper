
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
library(Biobase)
library(ggfortify)
`%ni%` <- Negate(`%in%`)

#Z score
Zscore_for      <- function(x){
  z <- (x - mean(x)) / sd(x)
  return(z)
}

#read gene expression data
load("./Data/Proximal_colon_gene_expression.RData", verbose = T)
  
data$id         <- paste0(data$strain, "_", data$diet)

genedata_res    <- data[, c("id", "gene_symbol", "value")]
genedata_res    <- reshape2::dcast(genedata_res, gene_symbol~id)
genedata_res    <- genedata_res[!is.na(genedata_res$gene_symbol),]
genedata_res    <- genedata_res[, !grepl("C57BL/6J|DBA/2J", colnames(genedata_res))]
rownames(genedata_res)   <- genedata_res$gene_symbol
exp             <- edgeR::DGEList(genedata_res[, colnames(genedata_res) %ni% "gene_symbol"])

exp$samples$strain <- gsub("_cd|_hfd|_CD|_HFD", "", rownames(exp$samples))
exp$samples$diet   <- gsub(".*_", "", rownames(exp$samples))
exp$samples$strain_diet <- rownames(exp$samples)

#normalization
exp <- calcNormFactors(exp, method = "TMM")

#filtering
palmieri_medians <- rowMedians(exp$counts)

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

man_threshold <- 5.5


hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

idx_man_threshold <- apply(exp$counts, 1,
                           function(x){
                             sum(x > man_threshold) >= min(length(exp$samples$strain[exp$samples$diet == "CD"]), length(exp$samples$strain[exp$samples$diet == "HFD"]))})
table(idx_man_threshold)
exp$counts  <- exp$counts[idx_man_threshold, ]

plot_tmp        <- as.data.frame(exp$counts)
plot_tmp_zscore <- as.data.frame(apply(plot_tmp, 1, Zscore_for))
plot_tmp_zscore <- as.matrix(t(plot_tmp_zscore))
type <- gsub(".*_", "", colnames(plot_tmp_zscore))

col = list(Diets = c("CD" = "#6495ED", "HFD" = "#DE3163"))
ha = HeatmapAnnotation(
  df = data.frame(Diets = type),
  annotation_height = unit(4, "mm"),
  col = col
)

#Figure 1—figure supplement 1C
pdf("./Plots/heatmap_cluster.pdf", width = 8, height = 8, useDingbats = FALSE)
ht  <- Heatmap(plot_tmp_zscore, name = "expression", km = 5, top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE) 
draw(ht)
dev.off()


exp$design <- model.matrix(~0 + diet + strain, data = exp$samples)
colnames(exp$design) <- make.names(colnames(exp$design))
exp$voom =voom(exp, design = exp$design, plot = T)

###PCA plot (Figure 1—figure supplement 1A)#######

exp$pca <- PCA(X = t(exp$voom$E), graph = F)
exp$pca$ind$coord <- cbind(exp$pca$ind$coord, exp$samples)
exp$pca$ind$coord <-  exp$pca$ind$coord[,unique(colnames(exp$pca$ind$coord))]
pca_tmp <- exp$pca$ind$coord

g <- ggplot(data = exp$pca$ind$coord, aes(x= Dim.1, y = Dim.2, col = diet)) + 
  geom_point(size = 3) + 
  #scale_colour_gradientn(colours= Heatmap_palette , limits=c(-2.5,2.5) ,na.value="gray87", oob = squish) +
  #geom_text(aes(label=strain_diet), color="black", size=3, vjust = 0.8, hjust = 0.5)  +
  stat_ellipse(level =0.9) +
  scale_color_manual(values = c("CD" = "#6495ED", "HFD" = "#DE3163")) + 
  xlab(label = paste0("PC1 (", signif(exp$pca$eig[1,2], 2), "%)")) +
  ylab(label = paste0("PC2 (", signif(exp$pca$eig[2,2], 2), "%)")) +
  theme_bw() +
  theme(axis.text.x = element_text(size =10, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11),
        strip.background = element_rect( fill = "white"),
        panel.grid = element_blank(),
        #panel.border = element_blank()
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
  ) 

strains          <- intersect(exp$samples$strain[exp$samples$diet == "CD"], exp$samples$strain[exp$samples$diet == "HFD"])
pca_tmp          <- pca_tmp[pca_tmp$strain %in% strains, ]
pca_tmp_order    <- pca_tmp[pca_tmp$diet == "HFD", ]
pca_tmp_order    <- pca_tmp_order[order(pca_tmp_order$Dim.1, decreasing = T), ]
pca_tmp$strain   <- factor(pca_tmp$strain, levels = unique(pca_tmp_order$strain))

#bar plot (Figure 1—figure supplement 1B)
g <- ggplot(data = pca_tmp, aes(x= strain, y = Dim.1, fill = diet)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values = c("CD" = "#6495ED", "HFD" = "#DE3163")) + 
  xlab(NULL) +
  ylab(label = "PC1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =9, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.ticks = element_blank(),
        legend.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11),
        strip.background = element_rect( fill = "white"),
        panel.grid = element_blank(),
        #panel.border = element_blank()
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
  ) 


#create contrast
level1.comp <- c("dietHFD - dietCD") 
contrastsm  <- makeContrasts(contrasts = level1.comp, levels =exp$design)
colnames(contrastsm) <- gsub("diet", "", colnames(contrastsm))

fit         <- lmFit(exp$counts, design = exp$design)
fitc        <- contrasts.fit(fit, contrasts = contrastsm)
fitc        <- eBayes(fitc)
decideT     <- decideTests(fitc, lfc = 0.5, p.value = 0.05, adjust.method = "BH")

summary(decideT)


tt <- lapply(colnames(contrastsm), function(x){
  out <- topTable(fitc, coef = x, n = Inf)
  out$contrast <- x
  dftoadd <- data.frame(gene_id = as.character(rownames(out)))
  out <- cbind(dftoadd, out)
  return(out)
})
names(tt) <- colnames(contrastsm)

save(tt, file="./Result/Differential_expression.RData")

################################################################################
## Volcano plots (Figure 1B)
################################################################################

library(cowplot)
library(ggrepel)

load(paste0("./Result/Differential_expression.RData"), verbose = T)
  # plot all volcano plots
  names(tt) <- c("HFD vs CD")
  vp <- mapply(x=tt, y=names(tt), function(x,y){ # x =tt$`HFD vs CD`  y = names(tt)[1]
    file <- gsub(" ", "", gsub("\\+", "_", y))
    data <- x %>%
      mutate( is_annotate = ifelse(((logFC > 0.5| logFC < -0.5) & adj.P.Val < 0.05), "yes", "no")) %>%
      mutate( is_highlight = ifelse(((logFC > 0.5| logFC < -0.5) & adj.P.Val < 0.05), "yes", "no"))
    
    sig     <- data[data$is_highlight == "yes", ]
    data$label <- data$gene_id
    data$label[data$is_annotate=="no"] <- NA
    data$label[grepl("^Gm|Rik$", data$label)] <- NA
    data$sort <- abs(data$logFC) * -log10(data$adj.P.Val)
    data <- data[order(data$sort, decreasing = T),]
    data$logPValue <- -log10(data$adj.P.Val)
    labels <- unique(data$label[!is.na(data$label)])

    if (length(labels) > 0) {
      value <- 0
      data1 <- data[data$logFC > 0.5,]
      data2 <- data[data$logFC < -0.5,]
      label1 <- as.character(data1$label[1:10])
      label2 <- as.character(data2$label[1:10])
      labels <- c(label1, label2)
      labels <- labels[!is.na(labels)]
      data$label[!data$label %in% labels] <- NA
    }
    
    g <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), label= label)) +
      geom_hex(data = subset(data, is_highlight=="no"), bins = 70) +
      geom_point(data=subset(data, is_highlight=="yes" & logFC >0), color= "#DE3163", size=2) +
      geom_point(data=subset(data, is_highlight=="yes" & logFC <0), color= "#6495ED", size=2) +
      scale_fill_gradientn(colors = (pals::brewer.greys(200)[100:200])) +
      geom_text_repel(size=3, parse = F) +
      geom_hline(yintercept = -log10(0.05), linetype='dotted')+
      geom_vline(xintercept = c(-0.5,0.5), linetype='dotted') +
      #ggtitle(label = y) + 
      theme_bw() + 
      theme(axis.text.x = element_text(size =10, color = "black"),
            plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
            axis.title = element_text(size = 13),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.ticks = element_blank(),
            legend.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 11),
            strip.background = element_rect( fill = "white"),
            panel.grid = element_blank(),
            #panel.border = element_blank()
            plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
            legend.position = "none"
      )
  })

  
  