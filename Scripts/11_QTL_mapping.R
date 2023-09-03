########################################################################
#clean.workplace
########################################################################
rm(list = ls())

########################################################################
#Library and Function
########################################################################
library(parallel)
library(reshape2)
library(clusterProfiler)
library(dplyr)
library(data.table)
library(tidyverse)
library(sva)
library(plotly)
library(readxl)
library(GenomicRanges)
library(foreach)
library(car)
library(preprocessCore)
library(qtl2)
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)

`%ni%` <- Negate(`%in%`)
########################################################################
#Data input and output paths
########################################################################
geno_raw <- as.data.frame(openxlsx::read.xlsx(paste0("./Data/fasted_genotype.xlsx")))
colnames(geno_raw) <- gsub("C57BL/6J","C57BL.6J",  gsub( "DBA/2J", "DBA.2J",colnames(geno_raw)))
geno_raw <- geno_raw[geno_raw$Locus %ni% c("UNC5348732", "rs6377403"),]

Diet_expression  <- c("HFD","CD")

lapply(Diet_expression, function(d){ #d= Diet_expression[2]
    print(d)
  data <- readRDS(paste0("./Data/WGCNA/", d, "_net.RDS"))
  data$analysis <- d
  
  allsamples <- colnames(data$dataset)
  MEs2       <- data.matrix(data$MEs)
  MEs2       <- MEs2[, colnames(MEs2) %ni% "ME0"]
  #remove outliers
    pheno <- colnames(MEs2)
    
    tmp_lt <- lapply(pheno, function(x){# x= colnames(MEs2)[1]
     # print(x)
      data_tmp <- data.frame(MEs2[,x])
      colnames(data_tmp) <- x
      data_tmp$id <- rownames(MEs2)
      data_tmp <- reshape2::melt(data_tmp, id.var = "id")
      data_tmp$value <- gsub("dead|outlier", NA, gsub("\\,", ".", data_tmp$value))
      data_tmp <- data_tmp[!is.na(data_tmp$value),]
      data_tmp$value <- as.numeric(as.character(data_tmp$value))
      IQR <- IQR(data_tmp$value)
      data_out <- data_tmp[(data_tmp$value > (quantile(data_tmp$value, probs = c(0.75)) + 1.5*IQR))|(data_tmp$value < (quantile(data_tmp$value, probs = c(0.25)) - 1.5*IQR)),]
      data_tmp <- data_tmp[data_tmp$id %ni% data_out$id,  ]
      return(data_tmp)
      })
    
    tmp_lt    <- do.call(rbind, tmp_lt)
    phenotype <- reshape2::dcast(tmp_lt[,c("id", "variable", "value")], id~variable)
    rownames(phenotype) <- phenotype$id
  
#data normalization
  normalized_data <- lapply(pheno, function(x){# x= colnames(MEs2)[1]
   # print(x)

    pheno.data = phenotype[, x]
    names(pheno.data) <- rownames(phenotype)
    # test if there is any zero or negative values in the phenotype data
    if(length(which(pheno.data <= 0)) != 0){
      # if non-positive value exists, use quantile transformation
      #message(paste0("Phenotype has non-positive values, quantile transformation will be applied!"))
      quantNorm = function(x){
        qnorm(rank(x, na.last = "keep",ties.method = "average")/(length(x)+1))
      }
      pheno.data.norm <- quantNorm(pheno.data)
    }else{
      # if all values are positive, boxcox transformation could be used
      #message(paste0("Phenotype has only positive values, Boxcox transformation will be applied!"))
      # scale all data into data that centers around 1, by dividing the average
      if (length(which(pheno.data >= 0)) != 0) {
        pheno.center <- pheno.data / mean(pheno.data, na.rm = TRUE)
        # run the box-cox transformation
        bc <- boxcox(pheno.center ~ 1, lambda = seq(-2, 2, 0.25), plotit=FALSE)
        trans <- bc$x[which.max(bc$y)]
        if(trans == 0){
          pheno.data.norm <- log(pheno.center)
        }else{
          pheno.data.norm <- (pheno.center^trans - 1)/trans
        }
      }else{
        pheno.data.norm <- NA
      }
      return(pheno.data.norm)
    }
  })
  
  normalized_data <- as.data.frame(do.call(cbind, normalized_data))
  colnames(normalized_data) <- pheno
  
    if(!dir.exists(paste0("./result/QTL_result/", d, "/Rawdata"))){
      dir.create(paste0("./result/QTL_result/", d, "/Rawdata"), recursive = T)
    }
    #######################################################################################
    #
    # eQTL preprocessing for rqtl2 packages
    #
    #######################################################################################
    
    system(paste0("cp ./Data/QTLmapping.yaml ", "./result/QTL_result/", d, "/Rawdata"))
  
    normalized_data$strain  <- gsub("_.*", "", rownames(normalized_data))
    genedata_res         <- normalized_data
    strain               <- unique(genedata_res$strain)
    geno                 <- geno_raw[,!colnames(geno_raw)%in%c("Chr","Mb_mm10","Mb_mm9","cM_BXD")]
    colnames(geno)       <- gsub( "C57BL.6J", "C57BL/6J", colnames(geno))
    colnames(geno)       <- gsub( "DBA.2J", "DBA/2J", colnames(geno))
    rownames(geno)       <- geno$Locus
    geno                 <- geno[, (colnames(geno) %in% strain)]
    geno$Locus           <- rownames(geno)
    setdiff(colnames(geno), strain)
    setdiff( strain, colnames(geno))
    geno <- reshape2::dcast(
      reshape2::melt(geno,"Locus"),
      variable ~ Locus, 
      value.var = "value"
    )
    colnames(geno)[colnames(geno)=="variable"] <- "id"
    geno$id <- as.character(geno$id)
    
    geno  <- cbind(id = geno$id, geno[,(colnames(geno) %in% geno_raw$Locus)])
    write.csv(
      x = geno,
      file = paste0("./result/QTL_result/", d, "/Rawdata/geno.csv"),
      quote = FALSE,
      na = "",
      row.names = FALSE
    )
    #gmap
    gmap <- geno_raw
    gmap <- gmap[,colnames(gmap) %in% c("Chr","Locus","cM_BXD")]
    colnames(gmap) <- c("chr","marker","pos")
    c <- 1:19
    c <- as.character(c)
    c <- c(c,"X")
    gmap_Order <- lapply(c, function(chr_ID){
      df <- gmap[gmap$chr==chr_ID,]
      df <- df[order(as.numeric(as.character(df$pos)),decreasing = FALSE),]
    })
    gmap_Order <- do.call(rbind, gmap_Order)
    gmap_Order <- gmap_Order[,c("marker","chr","pos")]
    
    write.csv(
      x = gmap_Order,
      file = paste0("./result/QTL_result/", d,"/Rawdata/gmap.csv"),
      quote = FALSE,
      na = "",
      row.names = FALSE
    )
    
    # covar
    covar <- data.frame(
      id = as.character(geno$id),
      sex = "m",
      cross_direction = "(BxD)x(BxD)"
    )
    
    write.csv(
      x = covar,
      file = paste0("./result/QTL_result/", d,"/Rawdata/covar.csv"),
      quote = FALSE,
      na = "",
      row.names = FALSE
    )
    
    #Expression Trait normal distribution
    rownames(genedata_res) <- normalized_data$strain
    genedata_res           <- genedata_res[, colnames(genedata_res) %ni% "strain"]
    pheno2                 <- as.data.frame(cbind(Strain = rownames(genedata_res),genedata_res))
    colnames(pheno2)[colnames(pheno2)=="Strain"] <- "id"
    pheno2 <- pheno2[pheno2$id %in% geno$id,]
    
    write.csv(
      x = pheno2,
      file = paste0("./result/QTL_result/", d,"/Rawdata/pheno.csv"),
      quote = FALSE, 
      row.names = FALSE
    )

    
    # Compute tissue specific and diet LMM with kinship
    cross2 <- read_cross2(file = paste0("./result/QTL_result/", d, "/Rawdata/QTLmapping.yaml"))
    map    <- insert_pseudomarkers(cross2$gmap, step=1)
    pr     <- calc_genoprob(cross = cross2, map = map, error_prob = 0.002, cores = 4)
    apr    <- genoprob_to_alleleprob(pr)
    grid   <- calc_grid(map = map, step=1)
    pr_grid <- probs_to_grid(pr, grid)
    kinship_grid <- calc_kinship(pr_grid, "loco")
    Xcovar <- get_x_covar(cross2)
    
    
    print("Compute LOD score using the LMM algorithm") 
    out_kinship <- scan1(
      genoprobs = pr_grid,
      pheno = cross2$pheno,
      kinship = kinship_grid,
      cores = 10
    )
    if(!dir.exists(paste0("./result/QTL_result/", d, "/Result"))){
      dir.create(paste0("./result/QTL_result/", d, "/Result"), recursive = T)
    }
 
    save(
      out_kinship,
      file = paste0( "./result/QTL_result/", d,"/Result/out_kinship_qtl2.RData")
    )
    
    print("Start Compute permutations") 
    print(Sys.time())
    
    operm <- scan1perm(pr_grid, 
                       cross2$pheno, 
                       kinship = kinship_grid, 
                       cores = 50,
                       n_perm = 10000,
    )
    print("finish Compute permutations") 
    print(Sys.time())
    save(
      operm,
      file = paste0(file_out, "/result/", d, "/Result/permutation_test.RData")
    )
    
    peaks <- find_peaks(out_kinship, map = map, threshold = summary(operm,alpha = 0.05), prob=0.95, peakdrop=5)
    if(nrow(peaks) > 0){
      peaks$analysis <- data$analysis
    }
    write.table(peaks, paste0(file_out, "/result/", d, "/Result/significant_qtl2_gene_position_0.05.txt"), quote=F, sep='\t', row.names=F)
    
    peaks <- find_peaks(out_kinship, map = map, threshold = summary(operm,alpha = 0.1), prob=0.95, peakdrop=5)
    if(nrow(peaks) > 0){
      peaks$analysis <- data$analysis
    }
    write.table(peaks, paste0(file_out, "/result/", d, "/Result/significant_qtl2_gene_position_0.1.txt"), quote=F, sep='\t', row.names=F)
    
    return(NULL)
})
})
