
########################################################################
#clean.workplace
########################################################################

rm(list = ls())
# clear all plots
dev.off()
# clear console
cat("\014")

library(PWMEnrich)
library(PWMEnrich.Mmusculus.background)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(MotifDb)
library(clusterProfiler)

genome <- BSgenome.Mmusculus.UCSC.mm10
standard.chr <- paste0("chr", c(1:19, "X", "Y"))

genes.gr <- keepSeqlevels(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), standard.chr, pruning.mode = "coarse")

# txdb <- AnnotationDbi::loadDb("./Data/JANSSENCOMMON/txdb.ref")
# seqlevelsStyle(txdb) <- "Ensembl"
# genes.gr <- keepSeqlevels(genes(txdb), standard.chr, pruning.mode = "coarse")

proms.gr <- promoters(genes.gr, upstream = 5000, downstream = 0)
gnconv <- bitr(proms.gr$gene_id, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
allproms.seq <- getSeq(genome,proms.gr)
allproms.seq <- allproms.seq[-which(grepl("N", allproms.seq))]
names(allproms.seq) <- gnconv$SYMBOL[match(names(allproms.seq), gnconv$ENTREZID)]

proms.gr$SYMBOL <- gnconv$SYMBOL[match(proms.gr$gene_id, gnconv$ENTREZID)]
pwm <- toPWM(as.list(subset (MotifDb, organism=='Mmusculus' & dataSource == "HOCOMOCOv10")))

registerCoresPWMEnrich(30)
useBigMemoryPWMEnrich(FALSE)

bg <- makePWMLognBackground(allproms.seq, pwm) # time-consuming

#download.file("http://hocomoco10.autosome.ru/final_bundle/MOUSE/mono/HOCOMOCOv10_annotation_MOUSE_mono.tsv", destfile = "./Data/JANSSENCOMMON/PWMenrich/HOCOMOCOv10_annotation_MOUSE_mono.tsv")

saveRDS(proms.gr, file = "./Data/WGCNA/PWMenrich/proms.gr.RDS")
saveRDS(allproms.seq, file = "./Data/WGCNA/PWMenrich/allproms.seq.RDS")
saveRDS(bg, file = "./Data/WGCNA/PWMenrich/pwm.hocomoco.logn.background.RDS")
