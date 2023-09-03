clusteredData <- list.files("./Data/WGCNA/", pattern = "_net\\.RDS$", full.names = T, recursive = T)
netList <- lapply(clusteredData, readRDS)
names(netList) <- gsub("\\/", "_", gsub("./Data/WGCNA/|_net.RDS", "", clusteredData))

netList <- mapply(x= netList,y = names(netList), function(x,y){
  x$analysis <- y
  x
}, SIMPLIFY = F)


library(clusterProfiler)
library(org.Mm.eg.db)
library(parallel)
library(msigdbr)

Genesets <-list(
  Hallmarks = msigdbr(species = "Mus musculus", category = "H")
)

#for Hallmark
mclapply(X = netList, function(net){ #net = netList[[4]]
  universe <- as.character(net$clusteringResult$gene)
  gene_clusters <- lapply(split(net$clusteringResult, net$clusteringResult$cluster), function(x){as.character(x$gene)})
  gene_clusters[["0"]] <- NULL
  
  Exp.Modules.En <-
    mapply(mod = names(gene_clusters), function(mod){ # mod = names(gene_clusters)[1]
      # mod = "c1_80"
      res1 <- mapply(geneset = names(Genesets), function(geneset){ # geneset = names(Genesets)[1]
        print(mod)
        print(geneset)
        
        if(length(gene_clusters[mod]) > 0){
          
          
          genes.names.tmp <- unlist(gene_clusters[mod])
          res             <- clusterProfiler::enricher(gene = genes.names.tmp,
                                                       TERM2GENE = Genesets[[geneset]][, c("gs_name", "gene_symbol")],
                                                       universe = universe,
                                                       minGSSize = 30,
                                                       maxGSSize = 800)
          return(res)
        }
      }, SIMPLIFY = FALSE)
      return(res1)
    }, SIMPLIFY = FALSE)
  

    saveRDS(Exp.Modules.En, paste0("./Data/WGCNA/Hallmark.RDS"))

  return(NULL)
}, mc.cores = 3)



load("./Data/gene_signitures_inflammtion_celltype.RData", verbose = T)
Genesets <- genesets_tt

mclapply(X = netList, function(net){ #net = netList[[1]]
  universe <- as.character(net$clusteringResult$gene)
  gene_clusters <- lapply(split(net$clusteringResult, net$clusteringResult$cluster), function(x){as.character(x$gene)})
  gene_clusters[["0"]] <- NULL
  
  Exp.Modules.En <-
    mapply(mod = names(gene_clusters), function(mod){ # mod = names(gene_clusters)[1]
      # mod = "c1_80"
      res1 <- mapply(geneset = names(Genesets), function(geneset){ # geneset = names(Genesets)[1]
        print(mod)
        print(geneset)
        
        if(length(gene_clusters[mod]) > 0){
          
          
          genes.names.tmp <- unlist(gene_clusters[mod])
          res             <- clusterProfiler::enricher(gene = genes.names.tmp,
                                                       TERM2GENE = Genesets[[geneset]][, c("match", "gene_symbol_mmusculus")],
                                                       #TERM2GENE = Genesets[[geneset]][, c("gs_name", "gene_symbol")],
                                                       universe = universe,
                                                       minGSSize = 20,
                                                       maxGSSize = 800)
          return(res)
        }
      }, SIMPLIFY = FALSE)
      return(res1)
    }, SIMPLIFY = FALSE)
  
    saveRDS(Exp.Modules.En, paste0("/home/li3/common/Users/xiaoxu/BXD_gut_project/WGCNA/Celltype_infla.RDS"))

  return(NULL)
}, mc.cores = 3)

### For IBD signature
genesets_mouse     <- read.table("./Data/custom_disease_signature_mouse.txt", header = T, sep = "\t", stringsAsFactors = F)
genesets_human     <- read.table("./Data/custome_disease_signature_human.txt", header = T, sep = "\t", stringsAsFactors = F)

genesets           <- rbind(genesets_mouse, genesets_human)

mclapply(X = netList, function(net){ #net = netList[[1]]
  universe <- as.character(net$clusteringResult$gene)
  gene_clusters <- lapply(split(net$clusteringResult, net$clusteringResult$cluster), function(x){as.character(x$gene)})
  gene_clusters[["0"]] <- NULL
  
  Exp.Modules.En <-
    mapply(mod = names(gene_clusters), function(mod){ # mod = names(gene_clusters)[1]
      # mod = "c1_80"
      print(mod)
      
      if(length(gene_clusters[mod]) > 0){
        
        
        genes.names.tmp <- unlist(gene_clusters[mod])
        res             <- clusterProfiler::enricher(gene = genes.names.tmp,
                                                     TERM2GENE = genesets[, c("condition", "GeneName")],
                                                     universe = universe,
                                                     minGSSize = 10,
                                                     maxGSSize = 800)
        return(res)
      }
      
    }, SIMPLIFY = FALSE)
  
    saveRDS(Exp.Modules.En, paste0("./Data/WGCNA/IBD.RDS"))
    
  return(NULL)
}, mc.cores = 1)


##############using CD modules as custom genesets##############

net <- readRDS(paste0("./Data/WGCNA/CD_net.RDS"))
  genesets <- net$clusteringResult
  genesets$MEcluster <- paste0("ME",genesets$cluster)
  genesets <- genesets[!grepl("ME0", genesets$MEcluster),]
  
  net_HFD <- readRDS(paste0("./Data/WGCNA/HFD_net.RDS"))
  universe <- as.character(net_HFD$clusteringResult$gene)
  gene_clusters <- lapply(split(net_HFD$clusteringResult, net_HFD$clusteringResult$cluster), function(x){as.character(x$gene)})
  gene_clusters[["0"]] <- NULL
  Exp.Modules.En <-
    mapply(mod = names(gene_clusters), function(mod){ # mod = names(gene_clusters)[1]
      # mod = "c1_80"
      print(mod)
      
      if(length(gene_clusters[mod]) > 0){
        
        
        genes.names.tmp <- unlist(gene_clusters[mod])
        res             <- clusterProfiler::enricher(gene = genes.names.tmp,
                                                     TERM2GENE = genesets[, c("MEcluster", "gene")],
                                                     universe = universe,
                                                     minGSSize = 10,
                                                     maxGSSize = 3000)
        return(res)
      }
      
    }, SIMPLIFY = FALSE)
  
  saveRDS(Exp.Modules.En, paste0("./Data/WGCNA/overlap_CD.RDS"))



