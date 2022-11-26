
karyotype_data = readRDS("chr_info_hg38.rds")
names(karyotype_data)[1] = "chr"
 karyotype_data$chr = gsub("chr","",karyotype_data$chr)
karyotype_data$color[karyotype_data$gieStain == "gneg"] = "white"
karyotype_data$color[karyotype_data$gieStain == "gpos25"] = "grey75"
karyotype_data$color[karyotype_data$gieStain == "gpos50"] = "grey50"
karyotype_data$color[karyotype_data$gieStain == "gpos75"] = "grey25"
karyotype_data$color[karyotype_data$gieStain == "gpos100"] = "grey0"
karyotype_data$color[karyotype_data$gieStain == "acen"] = "red"


# T2T
#cytobands_T2T.tsv
karyotype_data = read.table("cytobands_T2T.tsv", header=T, sep="\t")
names(karyotype_data)[1] = "chr"
#chr1        0  1735965 p36.33   gneg
#names(karyotype_data) = c("chr", "start","end","chr_location","gieStain")
 karyotype_data$chr = gsub("chr","",karyotype_data$chr)
karyotype_data$color[karyotype_data$gieStain == "gneg"] = "white"
karyotype_data$color[karyotype_data$gieStain == "gpos25"] = "grey75"
karyotype_data$color[karyotype_data$gieStain == "gpos50"] = "grey50"
karyotype_data$color[karyotype_data$gieStain == "gpos75"] = "grey25"
karyotype_data$color[karyotype_data$gieStain == "gpos100"] = "grey0"
karyotype_data$color[karyotype_data$gieStain == "acen"] = "red"
saveRDS(karyotype_data,file="cytobands_T2T.rds")


gene_coord = read.table("Homo_sapiens.GRCh38.93.gene.coord_strand_name.bed", header=F, sep="\t"); head(gene_coord) 
names(gene_coord) = c("chr", "start", "end", "strand", "gene")

save(list=ls(), file="sysdata.rda")
