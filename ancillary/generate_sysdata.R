
karyotype_data = readRDS("chr_info_hg38.rds")
names(karyotype_data)[1] = "chr"
 karyotype_data$chr = gsub("chr","",karyotype_data$chr)
karyotype_data$color[karyotype_data$gieStain == "gneg"] = "white"
karyotype_data$color[karyotype_data$gieStain == "gpos25"] = "grey75"
karyotype_data$color[karyotype_data$gieStain == "gpos50"] = "grey50"
karyotype_data$color[karyotype_data$gieStain == "gpos75"] = "grey25"
karyotype_data$color[karyotype_data$gieStain == "gpos100"] = "grey0"
karyotype_data$color[karyotype_data$gieStain == "acen"] = "red"


gene_coord = read.table("Homo_sapiens.GRCh38.93.gene.coord_strand_name.bed", header=F, sep="\t"); head(gene_coord) 
names(gene_coord) = c("chr", "start", "end", "strand", "gene")

save(list=ls(), file="sysdata.rda")
