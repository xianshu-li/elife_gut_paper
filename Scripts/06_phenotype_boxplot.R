########################################################################
#clean.workplace
########################################################################

rm(list = ls())

########################################################################
#library
########################################################################

library(cowplot)
library(plotly)
library(data.table)
library(RColorBrewer)
library(parallel)
library(ggplot2)
library(plyr) # For reformatting data
library(forcats) # For the "fct_reorder" function wich reorders factors on a numeric variable
library(scales) # For "oob =squish" Makes out-of-bounds values not considered as NA in the color scale
library(htmlwidgets)
library(Hmisc)
library(ggrepel)
library(ggpubr)
library(effsize)

`%ni%` <- Negate(`%in%`)

#read phenotypes
phenotypic_combine       <- read.table("./Data/phenotypes.txt", header = T, stringsAsFactors = F)

phenotype_name           <- c("Cytokine_IL15_.pg.mL.", "Cytokine_IL10_.pg.mL.")
strain_info              <- c("Strain_id", "diet", "strain")

cluster                  <- read.table("./Data/cluster_IBD_colon_mouse_new.txt", header = T)


lapply(phenotype_name, function(x){
  #x= "Cytokine_IL15_.pg.mL."
  data                 <- phenotypic_combine[, c("Strain_id", "strain", "diet", x)]
  colnames(data)       <- c("Strain_id", "strain", "diet", "value")

  data                 <- merge(data, cluster, by = "strain")
  data                 <- data[data$direction %in% c("Up", "Down", "Other"), ]
  data                 <- data[!is.na(data$value), ]
  
  g <- ggplot(data, aes(x=direction, y=value, fill=direction)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_boxplot(alpha = 0.5, outlier.color = "white") + 
    geom_jitter(aes(colour = direction), size = 0.5, alpha=0.5) +
    #stat_summary(fun.data=data_summary) +
    scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    scale_color_manual(values = c("Down" = "#4daf4a", "Up" = "#e41a1c", "Other" = "#3366cc")) +
    scale_fill_manual(values = c("Down" = "#4daf4a", "Up" = "#e41a1c", "Other" = "#3366cc")) +
    labs(subtitle= paste0(x, " (HFD)"),
         x="",
         y="") +
    stat_compare_means(comparisons = list(c("Down", "Up")), label = "p.signif",position = "identity",method = "t.test") +
    stat_compare_means(comparisons = list(c("Down", "Other")), label = "p.signif",position = "identity",method = "t.test") +
    stat_compare_means(comparisons = list(c("Up", "Other")), label = "p.signif",position = "identity",method = "t.test") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      strip.background = element_rect( fill = "white"),
      strip.text = element_text(size = 13)
      #panel.border = element_blank()
    ) 
})

