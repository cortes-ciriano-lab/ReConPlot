
#HG38
karyotype_data = readRDS("chr_info_hg38.rds")
names(karyotype_data)[1] = "chr"
karyotype_data$color[karyotype_data$gieStain == "gneg"] = "white"
karyotype_data$color[karyotype_data$gieStain == "gpos25"] = "grey75"
karyotype_data$color[karyotype_data$gieStain == "gpos50"] = "grey50"
karyotype_data$color[karyotype_data$gieStain == "gpos75"] = "grey25"
karyotype_data$color[karyotype_data$gieStain == "gpos100"] = "grey0"
karyotype_data$color[karyotype_data$gieStain == "acen"] = "red"
gene_coord = read.table("Homo_sapiens.GRCh38.93.gene.coord_strand_name.bed", header=F, sep="\t"); head(gene_coord)
names(gene_coord) = c("chr", "start", "end", "strand", "gene")

#HG19
karyotype_data_hg19 = readRDS("chr_info_hg19.rds")
names(karyotype_data_hg19)[1] = "chr"
# karyotype_data_hg19$chr =paste0("chr",karyotype_data_hg19$chr)
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "gneg"] = "white"
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "gpos25"] = "grey75"
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "gpos50"] = "grey50"
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "gpos75"] = "grey25"
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "gpos100"] = "grey0"
karyotype_data_hg19$color[karyotype_data_hg19$gieStain == "acen"] = "red"
gene_coord_hg19 = read.table("Homo_sapiens.GRCh37.87.genes.bed", header=F, sep="\t")
names(gene_coord_hg19) = c("chr", "start", "end", "gene")
gene_coord_hg19$chr <- paste0("chr", gene_coord_hg19$chr)

#T2T
#cytobands_T2T.tsv
library(tidyverse)
karyotype_data_T2T = read_tsv("cytobands_T2T.tsv")
names(karyotype_data_T2T)[1] = "chr"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "gneg"] = "white"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "gpos25"] = "grey75"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "gpos50"] = "grey50"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "gpos75"] = "grey25"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "gpos100"] = "grey0"
karyotype_data_T2T$color[karyotype_data_T2T$gieStain == "acen"] = "red"
saveRDS(karyotype_data_T2T,file="cytobands_T2T.rds")

gene_coord_T2T_raw = read_tsv("T2T_genes.tsv")
gene_coord_T2T_f=gene_coord_T2T_raw[which(gene_coord_T2T_raw$V30 !="N/A"),]
gene_coord_T2T_f=gene_coord_T2T_f[which(gene_coord_T2T_f$V18=="protein_coding"),]
gene_coord_T2T_ff=gene_coord_T2T_f[,c(1,2,3,13)]
names(gene_coord_T2T_ff) = c("chr", "start", "end", "gene")
gene_coord_T2T <- gene_coord_T2T_ff %>% 
  group_by(chr, gene) %>% 
  summarise(start=min(start), end=max(end)) %>% 
  ungroup() %>% 
  select(chr, start, end, gene)
rm(gene_coord_T2T_raw, gene_coord_T2T_f, gene_coord_T2T_ff)

#mm39
karyotype_data_mm39 = read_tsv("cytobands_mm39.tsv")
names(karyotype_data_mm39)[1:3] = c("chr", "start", "end")
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "gneg"] = "white"
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "gpos25"] = "grey75"
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "gpos50"] = "grey50"
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "gpos75"] = "grey25"
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "gpos100"] = "grey0"
karyotype_data_mm39$color[karyotype_data_mm39$gieStain == "acen"] = "red"
gene_coord_mm39_raw = read_tsv("genes_mm39.tsv")
gene_coord_mm39_f=gene_coord_mm39_raw[which(gene_coord_mm39_raw$V23=="protein_coding"),]
# gene_coord_mm39_f=gene_coord_mm39_raw[which(gene_coord_mm39_raw$V30 !="N/A"),]
gene_coord_mm39_ff=gene_coord_mm39_f[,c(1,2,3,18)]
names(gene_coord_mm39_ff) = c("chr", "start", "end", "gene")
gene_coord_mm39 <- gene_coord_mm39_ff %>% 
  group_by(chr, gene) %>% 
  summarise(start=min(start), end=max(end)) %>% 
  ungroup() %>% 
  select(chr, start, end, gene)
rm(gene_coord_mm39_raw, gene_coord_mm39_f, gene_coord_mm39_ff)

#mm10
karyotype_data_mm10 = read_tsv("cytobands_mm10.tsv")
names(karyotype_data_mm10)[1:3] = c("chr", "start", "end")
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "gneg"] = "white"
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "gpos33"] = "grey75"
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "gpos66"] = "grey50"
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "gpos75"] = "grey25"
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "gpos100"] = "grey0"
karyotype_data_mm10$color[karyotype_data_mm10$gieStain == "acen"] = "red"
gene_coord_mm10_raw = read_tsv("genes_mm10.tsv")
gene_coord_mm10_ff <- gene_coord_mm10_raw %>%
  select(-name)
names(gene_coord_mm10_ff) = c("chr", "start", "end", "gene")
gene_coord_mm10 <- gene_coord_mm10_ff %>% 
  group_by(chr, gene) %>% 
  summarise(start=min(start), end=max(end)) %>% 
  ungroup() %>% 
  select(chr, start, end, gene)
rm(gene_coord_mm10_raw, gene_coord_mm10_ff)

save(list=ls(), file="../R/sysdata.rda")
