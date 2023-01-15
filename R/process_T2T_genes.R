#d =read.table("T2T_genes.tsv",header=F,sep="\t")
#head(d)
#d=d[which(d$V30 !="N/A"),]
## protein coding only
#d=d[which(d$V18=="protein_coding"),]
#d =d[,c(1,2,3,6,13)]
#library(plyr)
#genes_T2T = ddply(d, .(V1,V13, V6), summarise,  start=min(V2), end=max(V3))
#genes_T2T =genes_T2T[,c(1,4,5,3,2)]
#
# names(genes_T2T) = c("chr", "start", "end", "strand", "gene")
#saveRDS(genes_T2T, "gene_coord_T2T.rds")
#
#
