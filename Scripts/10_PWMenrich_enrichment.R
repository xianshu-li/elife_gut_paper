clusteredData <- list.files("./Data/WGCNA/", pattern = "_net\\.RDS$", full.names = T, recursive = T)
netList <- lapply(clusteredData, readRDS)
names(netList) <- gsub("\\/", "_", gsub("./Data/WGCNA/|_net.RDS", "", clusteredData))

netList <- mapply(x= netList,y = names(netList), function(x,y){
  
  x$analysis <- y
  x
}, SIMPLIFY = F)

library(clusterProfiler)
library(org.Mm.eg.db)
library(PWMEnrich)
library(MotifDb)

bg <- readRDS("./Data/WGCNA/PWMenrich/pwm.hocomoco.logn.background.RDS")

proms.gr     <- readRDS("./Data/WGCNA/PWMenrich/proms.gr.RDS")
allproms.seq <- readRDS("./Data/WGCNA/PWMenrich/allproms.seq.RDS")


registerCoresPWMEnrich(30)

lapply(netList, function(net){ #net = netList[[6]]
  
  gene_clusters <- lapply(split(net$clusteringResult, net$clusteringResult$cluster), function(x){as.character(x$gene)})
  gene_clusters[["0"]] <- NULL
  
  gene_clusters_PWMenrich <- lapply(gene_clusters, function(x){ #x= gene_clusters[[25]]
    promsMatching <- match(x, names(allproms.seq))
    promsMatching <- promsMatching[!is.na(promsMatching)]
    curProms <- allproms.seq[promsMatching]
    out <- motifEnrichment(curProms, bg)
    out
  })
  gene_clusters_PWMenrich_matrix <- list(matrix = do.call(cbind, lapply(gene_clusters_PWMenrich, function(x){x$group.bg})),
                                         analysis = net$analysis)
  
  saveRDS(gene_clusters_PWMenrich_matrix, paste0("./Data/WGCNA/PWMenrich/", net$analysis, "_PWMenrichMatrix.RDS"))
  saveRDS(gene_clusters_PWMenrich, paste0("./Data/WGCNA/PWMenrich/", net$analysis, "_PWMenrichAllres.RDS"))
  return(NULL)
})
