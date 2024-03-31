#--------------------------------------------------------------------------------------------------------------------------
# formatter1000 <- function(x, step){
#   return(x/step)
# }

scaler <- function(step) {
  function(y) seq(0, ceiling(max(y)), by = step)
}

#-------------------------------------------------
# load dependencies
#-------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(scales)))
options(scipen=999)



#' Test structural variation data
#'
#' A dataset containing an example data frame with 
#' 538 SVs. The variables are as follows:
#'
#' \itemize{
#'   \item chr1. Chromosome 1 of the SV
#'   \item pos1. Position 1 of the SV
#'   \item chr2. Chromosome 2 of the SV
#'   \item pos2. Position 2 of the SV
#'   \item strands. Orientation of the SV (++, --, +- or -+)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sv_data
#' @usage data(sv_data)
#' @format A data frame with 538 rows and 5 variables
NULL

#' Test copy number data
#'
#' A dataset containing an example data frame with 
#' 2223 CN segments. The variables are as follows:
#'
#' \itemize{
#'   \item chr. Chromosome 
#'   \item start. Start
#'   \item end. End 
#'   \item copyNumber. Total copy number
#'   \item minorAlleleCopyNumber.  Minor allele copy number
#' }
#'
#' @docType data
#' @keywords datasets
#' @name cn_data
#' @usage data(cn_data)
#' @format A data frame with 2223 rows and 5 variables
NULL

#' Test single nucleotide variation data
#'
#' A dataset containing an example data frame with 
#' 100 SNVs to test the custom_annotation module. The variables are as follows:
#'
#' \itemize{
#'   \item chr. Chromosome 
#'   \item pos. Position
#'   \item y. Y coordinate (in this case, variant allele frequency) 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name snv_data
#' @usage data(snv_data)
#' @format A data frame with 100 rows and 3 variables
NULL


#' Function to plot genomic rearrangements (allele-specific copy number profiles and structural variants)
#'
#' This function generates publication-quality plots for copy number and structural variation profiles, which are particularly useful in the context of cancer genome analysis projects.
#' @import ggplot2
#' @import cowplot
#' @import dplyr
#' @import scales
#' @param sv Dataframe with SV information. Required, no default value.
#' @param cnv Dataframe with copy number information. Required, no default value.
#' @param title Title of the plot. Defaults to "".
#' @param genes Vector of genes (HUGO gene names) to be shown on the plot. Defaults to NULL.
#' @param chr_selection Chromosomes, start and end positions to be depicted. Defaults to NULL.
#' @param scaling_cn_SVs Relative dimension of the panel representing the SVs with respect to the copy number profile. Defaults to 1/6.
#' @param scale_separation_SV_type_labels Separatation (denoted as a fraction of the y axis XXX) between the SV labels/legend. Defaults to 1/18.
#' @param pos_SVtype_description Position for the SV labels/legend on the y axis. Defaults to 1000000.
#' @param window extra spacing on the x axis around the leftmost and rightmost breakpoints detected in the chromosomes selected unless a specific start and end positions for the plot are input.                      
#' @param xscale Scale for the x axis. Defaults to 10^6 to put the x axis in Mbp.
#' @param percentage_increase_y_axis Relative percentage to increase the distance of the y axis. XX
#' @param legend_SV_types Whether to show the legend for the SV types or not. Defaults to TRUE.
#' @param size_chr_labels Size of the chromsosome labels. Defaults to 7pt (use 5-7pt for publication-ready figures).
#' @param size_title Size of the plot title. Defaults to 7pt (use 5-7pt for publication-ready figures).
#' @param size_text Size of the plot text. Defaults to 5pt (use 5-7pt for publication-ready figures).
#' @param size_gene_label Size of the gene label text. Defaults to 1.5.
#' @param colour_band1 Colour of the first background horizontal stripe. Defaults to "grey86".
#' @param colour_band2 Colour of the second background horizontal stripe. Defaults to "antiquewhite1".
#' @param colour_DEL Colour of the arcs representing deletions (DEL). Defaults to "orange".
#' @param colour_h2hINV Colour of the arcs representing head-to-head inversions (h2hINV). Defaults to "forestgreen".
#' @param colour_DUP Colour of the arcs representing duplications (DUP). Defaults to "darkblue".
#' @param colour_t2tINV Colour of the arcs representing tail-to-tail inversions (t2tINV). Defaults to "black".
#' @param colour_TRA Colour of the arcs representing translocations (TRA). Defaults to "darkgray".
#' @param colour_INS Colour of the lines representing insertions (INS). Defaults to "darkgoldenrod4".
#' @param colour_SBE Colour of the arcs representing single breakends (SBE). Defaults to "bisque3".
#' @param color_minor_cn Colour for the horizontal bars representing the minor copy number values. Defaults to "#8491B4B2".
#' @param upper_limit_karyotype Upper limit in y-axis for karyotype ideogram. Defaults to -0.2 .
#' @param karyotype_rel_size Size of the karyotype ideogram, relative to the CN y axis. Defaults to .2.
#' @param custom_annotation Dataframe with custom annotation. Must contain columns "chr", "pos" and "y". Defaults to NULL.
#' @param ann_dot_col Colour of the dots in annotation plot. Defaults to "0.5."black".
#' @param ann_dot_size Size of the dots in annotation plot. Defaults to 0.5.
#' @param ann_y_title Label for Y-axis in annotation plot. Defaults to "".
#' @param ann_rel_size Size of the annotation plot, relative to the main plot. Defaults to 0.4.
#' @param ann_one_scale Plot the annotation y-axis between 0 and 1. Boolean, defaults to FALSE.
#' @param curvature_intrachr_SVs Curvature for the arcs represeting intrachromosomal SVs. Defaults to -0.15.
#' @param curvature_interchr_SVs Curvature for the arcs represeting interchromosomal SVs. Defaults to -0.08.
#' @param max.cn Cap on the total copy number. Defaults to 8.
#' @param npc_now Controls spacing, needed to find SV partners on different chromosomes in the plot. Not recommended to change. Defaults to .00625 * 3.
#' @param scale_ticks Spacing of breaks in the x axis (in bp). Defaults to 20000000 (i.e., 20Mb).
#' @param size_interchr_SV_tip Size of the line indicating interchromosomal SVs involving chromosomes not shown (that is, not included in chr_selection). Defaults to .1
#' @param label_interchr_SV Flag to annotate the second chromosome of interchromosomal SVs involving chromosomes not shown. Defaults to FALSE.
#' @param size_sv_line Linewidth for SVs. Defaults to .1.
#' @param genome_version Reference genome used. Can be either hg19, hg38, T2T, mm10 or mm39. Defaults to "hg38".
#' @examples
#' @export

ReConPlot <- function(sv,
                      cnv,
                      title="",
                      genes=NULL,
                      chr_selection = NULL,
                      scaling_cn_SVs = 1/5,
                      scale_separation_SV_type_labels = 1/22,
                      pos_SVtype_description = 5000000,
                      window = 10000000 ,
                      xscale = 1e6,
                      percentage_increase_y_axis=0.10, #05,
                      legend_SV_types=TRUE,
                      size_chr_labels=7,
                      size_title=7,
                      size_text=5,
                      size_gene_label=1.5,
                      colour_band1="grey86", colour_band2= "antiquewhite1", 
                      colour_DEL = "orange", colour_h2hINV="forestgreen", 
				              colour_DUP="darkblue", colour_t2tINV="black",
				              colour_TRA="darkgray", 
				              colour_INS = "darkgoldenrod4", colour_SBE="bisque3",
				              color_minor_cn="#8491B4B2",
				              upper_limit_karyotype = -0.2,
				              karyotype_rel_size=.1,
					            custom_annotation=NULL,
					            ann_dot_col="black",
					            ann_dot_size=.5,
					            ann_y_title="",
					            ann_rel_size=.4,
					            ann_one_scale=FALSE,
                      curvature_intrachr_SVs=-0.15,
                      curvature_interchr_SVs=-0.08, 
					            max.cn=8,
                      npc_now=.00625 * 3,
                      scale_ticks=20000000,
					            size_interchr_SV_tip=.1,
					            label_interchr_SV=FALSE,
					            size_sv_line=0.1,
				              genome_version="hg38"
                        ){

  #----------------------------------------------------------------
  # check input data
  #----------------------------------------------------------------
	if ( sum(!( chr_selection$chr %in% paste0("chr",c(1:22,"X","Y"))))  >0){
		stop("Error: the input data contains non-supported chromosomes. Supported chromosomes include the autosomes (chr1-chr22) and sexual chromsomes (chrX and chrY)")
	}
	if ( sum( ! c(sv$chr1, sv$chr2) %in% c(paste0("chr",c(1:22,"X","Y")), ".")) >0){
		 stop("Error: the input SV data contains non-supported chromosomes. Supported chromosomes include the autosomes (chr1-chr22) and sexual chromsomes (chrX and chrY)")
	}
	if ( sum( !( cnv$chr %in% paste0("chr",c(1:22,"X","Y")) >0))){
		 stop("Error: the input copy number data contains non-supported chromosomes. Supported chromosomes include the autosomes (chr1-chr22) and sexual chromsomes (chrX and chrY)")
	}
  required_SV_columns = c("chr1"  , "pos1" ,  "chr2"  , "pos2" ,  "strands")
  if (sum (! (required_SV_columns %in% names(sv)) ) > 0 ){
	  stop("Error: the input SV dataframe is missing required columns. The required columns are: chr1, pos1, chr2, pos2 and strands")
  }
  required_cn_columns = c("chr", "start", "end", "copyNumber","minorAlleleCopyNumber" )
  if (sum (! (required_cn_columns %in% names(cnv)) ) > 0 ){
	   stop("Error: the input copy number dataframe is missing required columns. The required columns are: chr, start, end, copyNumber, and minorAlleleCopyNumber")
  }
  if(is.data.frame(cnv) == FALSE){stop("Error: the copy number data must be input in dataframe format")}
  if(is.data.frame(sv) == FALSE){stop("Error: the SV data must be input in dataframe format")}
  #Had problems with tibbles
  cnv <- as.data.frame(cnv)
  sv <- as.data.frame(sv)

  #Check custom annotation dataframe if exists#
  if (! is.null(custom_annotation)){
    required_custom_columns=c("chr", "pos", "y")
    if (sum (! (required_custom_columns %in% names(custom_annotation)) ) > 0 ){
      stop("Error: the input Custom annotation dataframe is missing required columns. The required columns are: chr, pos, y")
    }
  }
  
  idx = which(sv$pos1 == sv$pos2)
  if(length(idx)>0){sv$pos2[idx] = sv$pos2[idx]+1}

  #----------------------------------------------------------------
  # information for the karyotype
  #----------------------------------------------------------------
  #Check genome version
  if (genome_version == "hg38") {
    karyotype_data_now = karyotype_data[karyotype_data$chr %in% chr_selection$chr,]
    karyotype_data_now$y = rep(1,nrow(karyotype_data_now))
    if (nrow(karyotype_data_now)<7){
      karyotype_data_annot=karyotype_data_now}
    else{
      karyotype_data_annot=karyotype_data_now[seq(3,(nrow(karyotype_data_now)-3),3),]}
    karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=unique(chr_selection$chr))}
  else if (genome_version == "hg19") {
    karyotype_data_now = karyotype_data_hg19[karyotype_data_hg19$chr %in% chr_selection$chr,]
    karyotype_data_now$y = rep(1,nrow(karyotype_data_now))
    if (nrow(karyotype_data_now)<7){
      karyotype_data_annot=karyotype_data_now}
    else{
      karyotype_data_annot=karyotype_data_now[seq(3,(nrow(karyotype_data_now)-3),3),]}
    karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=unique(chr_selection$chr))}
  else if (genome_version == "T2T") {
    karyotype_data_now = karyotype_data_T2T[karyotype_data_T2T$chr %in% chr_selection$chr,]
    karyotype_data_now$y = rep(1,nrow(karyotype_data_now))
    if (nrow(karyotype_data_now)<7){
      karyotype_data_annot=karyotype_data_now}
    else{
      karyotype_data_annot=karyotype_data_now[seq(3,(nrow(karyotype_data_now)-3),3),]}
    karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=unique(chr_selection$chr))}
  else if (genome_version == "mm10") {
    karyotype_data_now = karyotype_data_mm10[karyotype_data_mm10$chr %in% chr_selection$chr,]
    karyotype_data_now$y = rep(1,nrow(karyotype_data_now))
    if (nrow(karyotype_data_now)<7){
      karyotype_data_annot=karyotype_data_now}
    else{
      karyotype_data_annot=karyotype_data_now[seq(3,(nrow(karyotype_data_now)-3),3),]}
    karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=unique(chr_selection$chr))}
  else if (genome_version == "mm39") {
    karyotype_data_now = karyotype_data_mm39[karyotype_data_mm39$chr %in% chr_selection$chr,]
    karyotype_data_now$y = rep(1,nrow(karyotype_data_now))
    if (nrow(karyotype_data_now)<7){
      karyotype_data_annot=karyotype_data_now}
    else{
      karyotype_data_annot=karyotype_data_now[seq(3,(nrow(karyotype_data_now)-3),3),]}
    karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=unique(chr_selection$chr))}
  else {stop("Genome version not recognized")}
  

  #----------------------------------------------------------------
  # process SV data
  #----------------------------------------------------------------

  if (sum(sv$strands == "INS") >= 1) {
  sv.ins <- sv %>% filter(strands == "INS")
    sv.ins$colour <- colour_INS
  }
  
  if (sum(sv$strands == "SBE") >= 1) {
  sv.sbe <- sv %>% filter(strands == "SBE")
    sv.sbe$colour <- colour_SBE
  }

  
  sv <- sv %>% filter(strands != "SBE" & strands != "INS")
  
  sv$pos1 <- as.integer(sv$pos1)
  sv$pos2 <- as.integer(sv$pos2)
  #sv$ori = paste0(sv$strand1,sv$strand2)
  sv$chr1=as.character(sv$chr1)
  sv$chr2=as.character(sv$chr2)

  if (nrow(sv) >= 1){
    sv$colour <- apply(sv, 1, function(x){
      if (x["strands"] == "+-" | x["strands"] == "DEL"){
        c = colour_DEL
      } else if (x["strands"] == "++" | x["strands"] == "h2hINV"){
        c = colour_h2hINV
      } else if (x["strands"] == "--" | x["strands"] == "t2tINV"){
        c = colour_t2tINV
      } else if (x["strands"] == "-+" | x["strands"] == "DUP"){
        c = colour_DUP
      } else if (x["strands"] == "TRA"){
        c = colour_TRA
      }
      return(c)
    })
    sv$curve <- apply(sv, 1, function(x){
      if (as.character(x["chr1"]) == as.character(x["chr2"])){
      #JEVI: Note this part was creating problems in some of the Ly plots.
      #      I commented it out and c = curvature_intrachr_SVs always
        
        if (abs(as.integer(x["pos2"])-as.integer(x["pos1"])) <= 10000000){
          c = 2*curvature_intrachr_SVs
        } else {c = curvature_intrachr_SVs}
        # c = curvature_intrachr_SVs
      }
      else {c = curvature_intrachr_SVs}
      c
    })
  }

  #----------------------------------------------------------------
  # get intrachr SVs
  #----------------------------------------------------------------
  intraSV=c()
  for(i in chr_selection$chr){
    now = subset(sv, chr1 == i & chr2 == i)
    if(nrow(now)>0){intraSV = rbind(intraSV, now)}
  }
  
  #----------------------------------------------------------------
  # get interchr SVs
  #----------------------------------------------------------------
  # check if there are interchrSVs
  inter=sv[sv$chr1 != sv$chr2,]
  if( (sum(inter$chr1 %in% chr_selection$chr) + sum(inter$chr2  %in% chr_selection$chr)) >0){
    interFlag = TRUE
  } else {
    interFlag = FALSE
  }

  #--------------
  # get the x axis limits for all chrs
  #--------------
  for(i in 1:nrow(chr_selection)){
    chr.now = chr_selection$chr[i]
    start.now = chr_selection$start[i]
    end.now = chr_selection$end[i]
    if(is.na(start.now) | is.na(end.now)){
      # plot regions with SVs + window
      # find start
      pos1 = sv$pos1[sv$chr1==chr.now]; if(length(pos1)>0){pos1=min(pos1)}
      pos2 = sv$pos2[sv$chr2==chr.now]; if(length(pos2)>0){pos2=min(pos2)}
      start.now = c(pos1,pos2)
      if (nrow(sv.ins >= 1)) {
        pos1.ins <- sv.ins[sv.ins$chr1==chr.now]
        start.now <- c(start.now, pos1.ins)
      }
      if (nrow(sv.sbe >= 1)) {
        pos1.sbe <- sv.sbe[sv.sbe$chr1==chr.now]
        start.now <- c(start.now, pos1.sbe)
      }
      # which ones are numbers?
      idx= unlist(sapply(start.now, is.numeric))
      if(length(idx)>0){start.now=min(start.now[idx])}else{start.now=0}
      chr_selection$start[i] = start.now - window
      if(chr_selection$start[i]<0){chr_selection$start[i]=0}

      # now find the end for each chr if not provided
      end1 = sv$pos1[sv$chr1==chr.now]; if(length(end1)>0){end1=max(end1)}
      end2 = sv$pos2[sv$chr2==chr.now]; if(length(end2)>0){end2=max(end2)}
      end.now = c(end1,end2);
      # which ones are numbers?
      idx= unlist(sapply(end.now, is.numeric))
      if(length(idx)>0){end.now=max(end.now[idx])}else{end.now=250000000} # arbitrarily high number so the entire chr is displayed
      chr_selection$end[i] = end.now + window
    }
    #Set to maximum of chromosome if above it
    if (chr_selection$end[i] > max(karyotype_data_now$end[karyotype_data_now$chr == chr_selection$chr])){
      chr_selection$end[i] = max(karyotype_data_now$end[karyotype_data_now$chr == chr_selection$chr])
    }
  }
  if (interFlag){
    idx1= which(sv$chr1 != sv$chr2 & sv$chr1 %in% chr_selection$chr)
    idx2= which(sv$chr1 != sv$chr2 & sv$chr2 %in% chr_selection$chr)
    if (length(idx1)>0 | length(idx2)>0){
      interSV <- sv[unique(c(idx1,idx2)),] #  which(sv$chr1 != sv$chr2)  ,]
      interSV$pos1 <- as.integer(interSV$pos1)
      interSV$pos2 <- as.integer(interSV$pos2)
    }
    
    # get also the interSVs with other chromosomes
    idx1= which(sv$chr1 != sv$chr2 & sv$chr1 %in% chr_selection$chr &  !(sv$chr2 %in% chr_selection$chr))
    idx2= which(sv$chr1 != sv$chr2 & sv$chr2 %in% chr_selection$chr &  !(sv$chr1 %in% chr_selection$chr))
    if (length(idx1)>0 | length(idx2)>0){
      interSV_other_chrs <- sv[unique(c(idx1,idx2)),] #  which(sv$chr1 != sv$chr2)  ,]
      interSV_other_chrs$pos1 <- as.integer(interSV_other_chrs$pos1)
      interSV_other_chrs$pos2 <- as.integer(interSV_other_chrs$pos2)
    }

    
    #--------------
    # now get copy number for all chrs
    #--------------
    
    cnv.plot=c()
    for(i in 1:nrow(chr_selection)){
    main.cnv <- subset(cnv, chr == chr_selection$chr[i] &
                         ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                            (start <= chr_selection$start[i] & end >= chr_selection$start[i]) |
                            (start >= chr_selection$start[i] & start <= chr_selection$end[i])))
    main.cnv[main.cnv$start < chr_selection$start[i], "start"] <- chr_selection$start[i]
    main.cnv[main.cnv$end > chr_selection$end[i], "end"] <- chr_selection$end[i]
    cnv.plot <- rbind(main.cnv, cnv.plot)
    }


    #-----------------------------------
    # karyotype info now
    #-----------------------------------
    karyo.filt=c()
    for(i in 1:nrow(chr_selection)){ # note, there might be negative start values because we consider a window above
      karyo_subset <- subset(karyotype_data_now,
                             chr==chr_selection$chr[i] &
                               ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                                  (start <= chr_selection$start[i] & end >= chr_selection$end[i]) |
                                  (start >= chr_selection$start[i] & start <= chr_selection$end[i]) |
                                  (end >= chr_selection$start[i] & end <= chr_selection$end[i]) ))
      karyo_subset[karyo_subset$start < chr_selection$start[i],"start"] <- chr_selection$start[i]
      karyo_subset[karyo_subset$end > chr_selection$end[i],"end"] <- chr_selection$end[i]
      karyo.filt=rbind(karyo.filt,karyo_subset)
    }
    karyotype_data_now = karyo.filt
  }
    # there are no interSVs
  else {
    cnv.plot=c()
    for(i in 1:nrow(chr_selection)){
      main.cnv <- subset(cnv, chr == chr_selection$chr[i] &
                           ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                              (start <= chr_selection$start[i] & end >= chr_selection$start[i]) |
                              (start >= chr_selection$start[i] & start <= chr_selection$end[i])))
      main.cnv[main.cnv$start < chr_selection$start[i], "start"] <- chr_selection$start[i]
      main.cnv[main.cnv$end > chr_selection$end[i], "end"] <- chr_selection$end[i]
      cnv.plot <- rbind(main.cnv, cnv.plot)
    }
    #-----------------------------------
    # karyotype info now
    #-----------------------------------
    karyo.filt=c()
    for(i in 1:nrow(chr_selection)){ 
      # note, there might be negative start values because we consider a window above
      karyo_subset <- subset(karyotype_data_now,
                             chr==chr_selection$chr[i] &
                               ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                                  (start <= chr_selection$start[i] & end >= chr_selection$end[i]) |
                                  (start >= chr_selection$start[i] & start <= chr_selection$end[i]) |
                                  (end >= chr_selection$start[i] & end <= chr_selection$end[i]) ))
      karyo_subset[karyo_subset$start < chr_selection$start[i],"start"] <- chr_selection$start[i]
      karyo_subset[karyo_subset$end > chr_selection$end[i],"end"] <- chr_selection$end[i]
      karyo.filt=rbind(karyo.filt,karyo_subset)
    }
    karyotype_data_now = karyo.filt
  }
  
  # cnv.plot$chr = cnv.plot$chr
  # those with very high copy number value, modify
  cnv.plot$copyNumber[which(cnv.plot$copyNumber >= max.cn)] = max.cn
  max_y = max.cn # max(c( cnv.plot$copyNumber)) # max copy number value
  max_y_rectagle = max_y
  max_y = max_y + (max_y* percentage_increase_y_axis) # add 10% to y axis
  max_y_svs_1 = max_y + (scaling_cn_SVs* max_y)
  max_y_svs_2 = max_y + 2*(scaling_cn_SVs* max_y)
  max_y_svs_3 = max_y + 2.5*(scaling_cn_SVs* max_y) # for the interchr with chrs not displayed
  max_y_svs_4 = max_y + 2.55*(scaling_cn_SVs* max_y) # for the genes
  
  # max_y_svs_1 = max_y - 0.5
  #------------------------------------------------------------------------------------
  # plotting
  #------------------------------------------------------------------------------------
  p = ggplot()

  # horizontal rectangles
  # add background
  seq1=rep(seq(0,max_y_rectagle-1,2), each=length(chr_selection$chr))
  seq2=rep(seq(1,max_y_rectagle,2),each=length(chr_selection$chr))
  rectangle = data.frame(xmin=0,xmax=NA,
                          ymin=seq1,
                          ymax=seq2,
                          chr=factor(rep(chr_selection$chr,length(seq1)/length(chr_selection$chr)),levels=unique(chr_selection$chr)))
  # ## plot background
  p = p + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, group=chr),
                    alpha = .5, fill = colour_band1,
                    data = rectangle,
                    inherit.aes = F)

  p = p + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ymin+1, ymax = ymax+1, group=chr),
                    alpha = .5, fill = colour_band2,
                    data = rectangle,
                    inherit.aes = F)
  
  ## plot karyotype
  #karyotype_data_now$chr = as.vector(karyotype_data_now$chr)
  #karyotype_data_now$chr = factor(karyotype_data_now$chr, levels=chr_selection$chr)
  
  karyotype_space=max_y_rectagle*karyotype_rel_size
  lower_limit_karyotype=upper_limit_karyotype-(max_y_rectagle*karyotype_rel_size)
  p = p + geom_rect(data=karyotype_data_now,
                    mapping = aes(xmin = start, xmax = end, ymin = lower_limit_karyotype, ymax = upper_limit_karyotype, group=chr), # can change the default values for karyotype
                    fill = karyotype_data_now$color, color="black",size=.1)

   chromosome_labeller <- function(chr, value){
     chrm_name=gsub("chr", "", chr)
     return(paste0("Chromosome ", chrm_name, " (Mb)"))
     #vec = paste0("Chromosome ", chrm_name, " (Mb)")
	 #names(vec) = chr
     #return(vec) #paste0("Chromosome ", chrm_name, " (Mb)"))
   }

   
 #named_v =  paste0("sdf",chr_selection$chr) 
 # names(named_v) = chr_selection$chr
  #chromosome_labeller = as_labeller(named_v)

	#cnv.plot$chr = as.vector(cnv.plot$chr)
	#cnv.plot$chr_f = factor(as.vector(cnv.plot$chr), levels=chr_selection$chr)


  #define the breaks on the y axis (cn) to show:
   if (max.cn < 8) {break.step = 1}
   else if (max.cn >= 8 & max.cn <= 20) {break.step = 2}
   else if (max.cn > 20 & max.cn <= 50) {break.step = 5}
   else {break.step = 10}
   breaks_y = seq(0,max.cn, break.step)
  # breaks_y = c(0,1,2)
  # breaks_y = c(breaks_y, unique(floor(quantile(2:max_y, probs = seq(.1, .9, by = .1)))))
  ## plot copy number
  p = p +
    geom_segment(aes(x = start, y = minorAlleleCopyNumber, xend = end,
                     yend = minorAlleleCopyNumber, group=chr), # ICC group
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr,], colour=color_minor_cn,size=1) +
    geom_segment(aes(x = start, y = copyNumber, xend = end,
                     yend = copyNumber,  group=chr), # ICC group
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr,], colour="black",size=1) +
    facet_grid( . ~ factor(chr, levels = unique(chr_selection$chr)), # not converting to factor to make the labeller work
                scales = "free_x", space = "free_x", switch = "x", #) +
                # #Add different labels?
                 labeller = as_labeller(chromosome_labeller)) + #(chr=chr_selection$chr) ) +
                 #labeller = labeller(chr = chromosome_labeller )) +
                 #labeller = labeller(chr = chromosome_labeller(chr=chr_selection$chr) )) +
    ggtitle(title) +
    theme(#aspect.ratio = 1,
      text=element_text(size=size_text, colour="black"), #, face ="bold"),
      axis.text.x=element_text(size=size_text, colour="black"),
      axis.title.x =element_text(size=size_text, colour="black"),
      axis.title.y =element_text(size=size_text, colour="black"),
      plot.title =element_text(size=size_title, colour="black", margin = unit(c(-.15,0,0,0),'cm')), # to reduce space between title and plot:
      axis.text.y=element_text(size=size_text, colour="black"),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.clip = "off",
      plot.margin = unit(c(.25, .1, -0.1, .1), "cm"),
      panel.grid.major.x = element_blank(),
      #panel.spacing = unit(4, "lines"), ICC
      panel.spacing.x = unit(npc_now, "npc"),
      strip.text.x = element_text(size=size_chr_labels, face="bold"), #element_blank(), ## to remove the title of the facet panels
      strip.switch.pad.grid = unit(0, "cm"),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line( size=.05, color="black"),
      panel.border = element_blank() #rect(colour = "black", fill=NA, size=.1)
    ) +
    labs(x = "", y = "Copy number")  +
    coord_cartesian(clip = "off", expand=0) + # ICC
    #ylim(c(0,max_y),expand) +
    scale_y_continuous(limits = c(lower_limit_karyotype - 0.05,max_y_svs_3+0.5),  breaks=breaks_y, expand = 1,
                       labels=c(breaks_y[1:length(breaks_y)-1], paste0(">=",breaks_y[length(breaks_y)]))) + #seq(0,max_y,2)) +
    scale_x_continuous(expand=c(0,0),
                       labels=number_format(scale=1/xscale),
                       breaks=scaler(scale_ticks))
  
  # add lines where the SVs go
  p = p + geom_abline(slope = 0, intercept=max_y_svs_1, colour="black", size=.25)
  p = p + geom_abline(slope = 0, intercept=max_y_svs_2, colour="black", size=.25)

  

  #----------------------------------------------------------------------------
  # add intrachr
  #----------------------------------------------------------------------------
  #if (nrow(intraSV) >= 1){
  if (!is.null(intraSV)){
    if (nrow(intraSV) >= 1){
      for (i in 1:nrow(intraSV)){
        # get the maximum coordinates for the chr at hand
        max_coord = max(cnv.plot[cnv.plot$chr==intraSV$chr1[i],"end"])
        min_coord = min(cnv.plot[cnv.plot$chr==intraSV$chr1[i],"start"])
        # add max
        intraSV$max_y_sv = max_y_svs_2
        intraSV$max_y_sv[which(intraSV$strands %in% c("DEL","DUP","+-", "-+"))] = max_y_svs_1
        # correct the arc;
        intraSV$curve[which( intraSV$strands %in% c("DEL","h2hINV", "+-", "--"))] = abs(intraSV$curve[which( intraSV$strands %in% c("DEL","h2hINV", "+-", "--"))])
        
        # if one breakpoint outside or range, just plot a vertical line, if not the entire arc
        if (intraSV$pos1[i] >= min_coord & intraSV$pos2[i] <= max_coord ){ # within range
          p = p +
            # arc
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos1[i], xend=intraSV$pos2[i],
                       y=intraSV$max_y_sv[i], yend=intraSV$max_y_sv[i], curvature=intraSV$curve[i],
                       size=size_sv_line,colour=intraSV$colour[i])
          p = p +
            # vertical
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos1[i], xend=intraSV$pos1[i],
                       y=0, yend=intraSV$max_y_sv[i], curvature=0,  size=size_sv_line, colour=intraSV$colour[i]) +
            # vertical
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos2[i], xend=intraSV$pos2[i],
                       y=0, yend=intraSV$max_y_sv[i], curvature=0,  size=size_sv_line, colour=intraSV$colour[i])
        }
      
        # cases now with left size in range
        if (intraSV$pos1[i] >= min_coord & intraSV$pos1[i] < max_coord  &  intraSV$pos2[i] > max_coord){ # within range
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos1[i], xend=intraSV$pos1[i],
                       y=0, yend=max_y_svs_3-0.15, curvature=0,  size=size_sv_line, colour=intraSV$colour[i])
          # add diagonal line on top
          x_range = (chr_selection$end[which(chr_selection$chr==intraSV$chr1[i])] - chr_selection$start[which(chr_selection$chr==intraSV$chr1[i])]) * 0.01
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos1[i], xend=intraSV$pos1[i]+x_range, #2000000,
                       y=max_y_svs_3-size_interchr_SV_tip, yend= max_y_svs_3, angle=45, curvature=0,  size=size_sv_line, colour=intraSV$colour[i])
          if (label_interchr_SV) {
            p = p + geom_text(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                              x=intraSV$pos1[i]-x_range, #2000000,
                              y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=intraSV$colour[i],
                              label=intraSV$chr2[i])
          }
        }

        # `    # cases now with right size in range
        if (intraSV$pos1[i] < min_coord &  intraSV$pos2[i] <= max_coord & intraSV$pos2[i] > min_coord){ # within range
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x=intraSV$pos2[i], xend=intraSV$pos2[i],
                       y=0, yend=max_y_svs_3-0.15, curvature=0,  size=size_sv_line, colour=intraSV$colour[i])
          # add diagonal line on top
          x_range = (chr_selection$end[which(chr_selection$chr==intraSV$chr1[i])] - chr_selection$start[which(chr_selection$chr==intraSV$chr1[i])]) * 0.01
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr2[i]),
                       x=intraSV$pos2[i], xend=intraSV$pos2[i]-x_range, #2000000,
                       y=max_y_svs_3-size_interchr_SV_tip, yend=max_y_svs_3 , angle=45, curvature=0,  size=size_sv_line, colour=intraSV$colour[i])
          if (label_interchr_SV) {
            p = p + geom_text(data = data.frame(cov = 1, chr = intraSV$chr2[i]),
                              x=intraSV$pos2[i]-x_range, #2000000,
                              y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=intraSV$colour[i],
                              label=intraSV$chr1[i])
          }
        }
      }
    }
  }
  
  #----------------------------------------------------------------
  # plot interchr SVs involving chrs in levels_chrs
  #----------------------------------------------------------------
  if (interFlag){
    # # first those involving the middle chromosome  ## missing those with other chrs
    idx= which(interSV$chr1 %in% chr_selection$chr & interSV$chr2 %in% chr_selection$chr)

    if(length(idx)>0){
      interSV=interSV[idx,]
      # add max
      interSV$max_y_sv = max_y_svs_2
      interSV$max_y_sv[which( interSV$strands %in% c("DEL","DUP","+-", "-+"))] = max_y_svs_1
      # change curvature
      interSV$curve[which( interSV$strands %in% c("DEL","h2hINV", "+-", "--"))] = abs(interSV$curve[which( interSV$strands %in% c("DEL","h2hINV", "+-", "--"))])
      info_chrs = data.frame(chrs=chr_selection$chr)
      # determine offset: we do so by: (1) compute the total size of the chrs displayed
      total_chr_size = 0; size_in_plot=c()
      for(ii in chr_selection$chr){
        size_now=max(cnv.plot[cnv.plot$chr==ii,"end"]) -min(cnv.plot[cnv.plot$chr==ii,"start"])
        size_in_plot= c(size_in_plot, size_now)
        total_chr_size = total_chr_size +  size_now
      }
      info_chrs$size_in_plot  = size_in_plot
      # now we compute the fraction of the plot that is "space", which corresponds to approx npc_now
      gap = (total_chr_size*npc_now)
      for (i in 1:nrow(interSV)){  ## check as well that the chrs interact with the selected one!
        # in whcih position of the chr list we find the main chr?
        # we assume as main chr the first in the list of chrs
        position_chr_1 = which(chr_selection$chr == interSV$chr1[i])
        position_chr_2 = which(chr_selection$chr == interSV$chr2[i])

        # now define which chr is the leftmost in the plot (that is, in chr_selection$chr), so we always draw SVs to the right
        #if (position_chr_1 < position_chr_2){} # the first chr in the SV object is the leftmost, no need to do anything
        if (position_chr_1 > position_chr_2){
          chr1_now=interSV$chr1[i]
          pos1_now=interSV$pos1[i]
          chr2_now=interSV$chr2[i]
          pos2_now=interSV$pos2[i]
          # change the order in object, so below we do not need to worry about which column to use depending on chr order in the plot
          interSV$chr1[i] = chr2_now; interSV$chr2[i] = chr1_now; interSV$pos1[i] = pos2_now; interSV$pos2[i] = pos1_now
          chrs_now = interSV[i,c("chr1","chr2")]
        }else{chrs_now = interSV[i,c("chr1","chr2")]}
        #----------------------------
        # we plot always to the right now
        #----------------------------
        # range in second chr
        min_range=  min(cnv.plot[cnv.plot$chr==chrs_now[1,2],"start"])
        size_leftmost_chr = max(cnv.plot[cnv.plot$chr==chrs_now[1,1],"end"])

        # compute offset
        offset = interSV$pos1[i] +  # the start in leftmost
          (size_leftmost_chr - interSV$pos1[i]) + # what we need to go of the leftmost chr till the end
          (interSV$pos2[i] - min_range) + # the position in the chr on the right + the value where it starts
          gap
        # are there other chrs in between the origin chr (main chr) and the target chr, which in this case is on the right??
        indexes = seq(1,length(chr_selection$chr))
        indexes = which(indexes>position_chr_1 & indexes<position_chr_2)
        if(length(indexes)>0){offset = offset + (gap * length(indexes)) + sum(info_chrs$size_in_plot[indexes])}

        # check that both breakpoints map to the copy number interval in the plot
        min_pos_chr1 =  min(cnv.plot[cnv.plot$chr==chrs_now[1,1],"start"])
        max_pos_chr1 =  max(cnv.plot[cnv.plot$chr==chrs_now[1,1],"end"])
        min_pos_chr2 =  min(cnv.plot[cnv.plot$chr==chrs_now[1,2],"start"])
        max_pos_chr2 =  max(cnv.plot[cnv.plot$chr==chrs_now[1,2],"end"])
        in_range_chr1 = (interSV$pos1[i] > min_pos_chr1 & interSV$pos1[i] < max_pos_chr1)
        in_range_chr2 = (interSV$pos2[i] > min_pos_chr2 & interSV$pos2[i] < max_pos_chr2)
        if(in_range_chr1 & in_range_chr2){
          # add vertical lines first
          p = p +
            # on the left
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x=interSV$pos1[i], xend=interSV$pos1[i],
                       y=0, yend=interSV$max_y_sv[i], curvature=0,  size=size_sv_line, colour=interSV$colour[i]) +
            # on the right
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x=interSV$pos2[i], xend=interSV$pos2[i],
                       y=0, yend=interSV$max_y_sv[i], curvature=0,  size=size_sv_line, colour=interSV$colour[i])
          # add arc
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = chrs_now[1,1]),
                       x=interSV$pos1[i], # the start in chr9
                       xend=offset,
                       y=interSV$max_y_sv[i], yend=interSV$max_y_sv[i], curvature=interSV$curve[i],
                       size=size_sv_line,colour=interSV$colour[i])
        }
        # the leftmost breakpoint is outside of the range
        if(!in_range_chr1){
          ## add vertical line
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x=interSV$pos2[i], xend=interSV$pos2[i],
                       y=0, yend=max_y_svs_3-size_interchr_SV_tip, curvature=0,  size=size_sv_line, colour=interSV$colour[i])
          # add diagonal line on top
          x_range = (chr_selection$end[which(chr_selection$chr==interSV$chr1[i])] - chr_selection$start[which(chr_selection$chr==interSV$chr1[i])]) * 0.01
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x=interSV$pos2[i], xend=interSV$pos2[i]-x_range, #2000000,
                       y=max_y_svs_3-size_interchr_SV_tip, yend=max_y_svs_3 , angle=45, curvature=0,  size=size_sv_line, colour=interSV$colour[i])
          if (label_interchr_SV) {
            p = p + geom_text(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                              x=interSV$pos2[i]-x_range, #2000000,
                              y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=interSV$colour[i],
                              label=interSV$chr1[i])
          }
        }

        # the rightmost breakpoint is outside of the range
        if(!in_range_chr2){
          ## add vertical line
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x=interSV$pos1[i], xend=interSV$pos1[i],
                       y=0, yend=max_y_svs_3-size_interchr_SV_tip, curvature=0,  size=size_sv_line, colour=interSV$colour[i])
          # add diagonal line on top
          x_range = (chr_selection$end[which(chr_selection$chr==interSV$chr1[i])] - chr_selection$start[which(chr_selection$chr==interSV$chr1[i])]) * 0.01
          p = p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x=interSV$pos1[i], xend=interSV$pos1[i]-x_range, #2000000,
                       y=max_y_svs_3-size_interchr_SV_tip, yend=max_y_svs_3 , angle=45, curvature=0,  size=size_sv_line, colour=interSV$colour[i])
          if (label_interchr_SV) {
            p = p + geom_text(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                              x=interSV$pos1[i]-x_range, #2000000,
                              y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=interSV$colour[i], 
                              label=interSV$chr2[i])
          }
        }
      }
    }
  }
  
  #----------------------------------------------------------------
  # plot lines with ticks for interchr SVs involving chrs not displayed
  #----------------------------------------------------------------
  
  if(interFlag){ #(interSV_other_chrs)){
    if (exists("interSV_other_chrs")){
      for(i in 1:nrow(interSV_other_chrs)){
  
        chr1_in = (interSV_other_chrs$chr1[i] %in% chr_selection$chr)
  
        # the leftmost breakpoint is outside of the range
        if((!chr1_in) ){
          # check that both breakpoints map to the copy number interval in the plot for chr2
          min_pos_chr2 =  min(cnv.plot[cnv.plot$chr==interSV_other_chrs$chr2[i],"start"])
          max_pos_chr2 =  max(cnv.plot[cnv.plot$chr==interSV_other_chrs$chr2[i],"end"])
          in_range_chr2 = (interSV_other_chrs$pos2[i] > min_pos_chr2 & interSV_other_chrs$pos2[i] < max_pos_chr2)
          if(in_range_chr2){
            ## add vertical line
            p = p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr2[i]),
                         x=interSV_other_chrs$pos2[i], xend=interSV_other_chrs$pos2[i],
                         y=0, yend=max_y_svs_3-size_interchr_SV_tip, curvature=0,  size=size_sv_line, colour=interSV_other_chrs$colour[i])
            # add diagonal line on top
            x_range = (chr_selection$end[which(chr_selection$chr==interSV_other_chrs$chr2[i])] - chr_selection$start[which(chr_selection$chr==interSV_other_chrs$chr2[i])]) * 0.01
            p = p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr2[i]),
                         x=interSV_other_chrs$pos2[i], xend=interSV_other_chrs$pos2[i]-x_range, #2000000,
                         y=max_y_svs_3-size_interchr_SV_tip, yend=max_y_svs_3 , angle=45, curvature=0,  size=size_sv_line, colour=interSV_other_chrs$colour[i])
          if (label_interchr_SV) {
              p = p + geom_text(data = data.frame(cov = 1, chr = interSV_other_chrs$chr2[i]),
                                x=interSV_other_chrs$pos2[i]-x_range, #2000000,
                                y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=interSV_other_chrs$colour[i],
                                label=interSV_other_chrs$chr1[i])
            }
          }
        }
  
        # the rightmost breakpoint is outside of the range
        chr2_in = (interSV_other_chrs$chr2[i] %in% chr_selection$chr)

        if( (!chr2_in)){
          # check that both breakpoints map to the copy number interval in the plot for chr1
          min_pos_chr1 =  min(cnv.plot[cnv.plot$chr==interSV_other_chrs$chr1[i],"start"])
          max_pos_chr1 =  max(cnv.plot[cnv.plot$chr==interSV_other_chrs$chr1[i],"end"])
          in_range_chr1 = (interSV_other_chrs$pos1[i] > min_pos_chr1 & interSV_other_chrs$pos1[i] < max_pos_chr1)
          if( in_range_chr1){
            # add vertical line
            p = p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr1[i]),
                         x=interSV_other_chrs$pos1[i], xend=interSV_other_chrs$pos1[i],
                         y=0, yend=max_y_svs_3-size_interchr_SV_tip, curvature=0,  size=size_sv_line, colour=interSV_other_chrs$colour[i])
            # add diagonal line on top
            x_range = (chr_selection$end[which(chr_selection$chr==interSV_other_chrs$chr1[i])] - chr_selection$start[which(chr_selection$chr==interSV_other_chrs$chr1[i])]) * 0.01

            p = p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr1[i]),
                         x=interSV_other_chrs$pos1[i], xend=interSV_other_chrs$pos1[i]-x_range, #2000000,
                         y=max_y_svs_3-size_interchr_SV_tip, yend=max_y_svs_3 , angle=45, curvature=0,  size=size_sv_line, colour=interSV_other_chrs$colour[i])
            if (label_interchr_SV) {
              p = p + geom_text(data = data.frame(cov = 1, chr = interSV_other_chrs$chr1[i]),
                                x=interSV_other_chrs$pos1[i]-x_range, #2000000,
                                y=max_y_svs_3+size_interchr_SV_tip, size=size_text/.pt, colour=interSV_other_chrs$colour[i], 
                                label=interSV_other_chrs$chr2[i])
            }
          }
        }
      }
    }
  }
  #----------------------------------------------------------------------------
  # plot SBE and insertions
  #----------------------------------------------------------------------------  
  if (nrow(sv.ins >= 1)){ 
    for (i in 1:nrow(sv.ins)) {
      min_coord = min(cnv.plot[cnv.plot$chr==sv.ins$chr1[i],"start"])
      max_coord = max(cnv.plot[cnv.plot$chr==sv.ins$chr1[i],"end"])
      # add max
      sv.ins$max_y_sv = max_y_svs_3-size_interchr_SV_tip
      # correct the arc;
      # if one breakpoint outside or range, just plot a vertical line, if not the entire arc
      if (sv.ins$pos1[i] >= min_coord & sv.ins$pos1[i] <= max_coord) { # within range
        p = p +
          # vertical
          geom_curve(data = data.frame(cov = 1, chr = sv.ins$chr1[i]),
                     x=sv.ins$pos1[i], xend=sv.ins$pos1[i],
                     y=0, yend=sv.ins$max_y_sv[i], 
                     curvature=0, size=size_sv_line, colour=sv.ins$colour[i]) +
           geom_point(data= data.frame(cov = 1, chr = sv.ins$chr1[i]), 
                      x= sv.ins$pos1[i], y = sv.ins$max_y_sv[i],
                      size=1, shape = 24, colour="black", stroke = .1, fill = sv.ins$colour[i])
    }}}
  if (nrow(sv.sbe >= 1)){ 
    for (i in 1:nrow(sv.sbe)) {
      min_coord = min(cnv.plot[cnv.plot$chr==sv.sbe$chr1[i],"start"])
      max_coord = max(cnv.plot[cnv.plot$chr==sv.sbe$chr1[i],"end"])
      # add max
      sv.sbe$max_y_sv = max_y_svs_3-size_interchr_SV_tip
      # correct the arc;
      # if one breakpoint outside or range, just plot a vertical line, if not the entire arc
      if (sv.sbe$pos1[i] >= min_coord & sv.sbe$pos1[i] <= max_coord) { # within range
        p = p +
          # vertical
          geom_curve(data = data.frame(cov = 1, chr = sv.sbe$chr1[i]),
                     x=sv.sbe$pos1[i], xend=sv.sbe$pos1[i],
                     y=0, yend=sv.sbe$max_y_sv[i], 
                     curvature=0, size=size_sv_line, colour=sv.sbe$colour[i]) +
          geom_point(data= data.frame(cov = 1, chr = sv.sbe$chr1[i]), 
                     x= sv.sbe$pos1[i], y = sv.sbe$max_y_sv[i],
                     size=1, shape = 22, colour="black", stroke = .1, fill = sv.sbe$colour[i])
      }}}
    
    
  
  
  
  #----------------------------------------------------------------------------
  # highlight genes
  #----------------------------------------------------------------------------
  if(!is.null(genes)){
    for(gene in genes){
      if (genome_version == "hg38") {
        gene_coord_now = gene_coord[gene_coord$gene==gene,]
      }
      else if (genome_version == "hg19") {
        gene_coord_now = gene_coord_hg19[gene_coord_hg19$gene==gene,]
      }
      else if (genome_version == "T2T") {
        gene_coord_now = gene_coord_T2T[gene_coord_T2T$gene==gene,]
      }
      else if (genome_version == "mm10") {
        gene_coord_now = gene_coord_mm10[gene_coord_mm10$gene==gene,]
      }
      else if (genome_version == "mm39") {
        gene_coord_now = gene_coord_mm39[gene_coord_mm39$gene==gene,]
      }
      else {stop("Genome version not recognized")}
      #Some genes are not in the annotation list -> skip
      if (nrow(gene_coord_now) == 0) {next}
      # are the genes to be plotted in the selected chromosomes?

      if (gene_coord_now$chr %in% chr_selection$chr){
        if(nrow(gene_coord_now)>0 & gene_coord_now$chr %in% chr_selection$chr &
           gene_coord_now$start > chr_selection$start[which(chr_selection$chr==gene_coord_now$chr)] & # if the gene is outside plotting range discard
           gene_coord_now$start < chr_selection$end[which(chr_selection$chr==gene_coord_now$chr)]
        ){
          
          # add vertical line
          p =  p +
            #Substitute line for geom rect,
            #Keep line for small genes in large plotting regions geom_rect sometimes fails
            geom_curve(data = data.frame(cov = 1, chr = gene_coord_now$chr),
                       x=gene_coord_now$start, xend=gene_coord_now$start,
                       y=0, #max_y_svs_2 + (max_y_svs_3-max_y_svs_2)/2,
                       yend=max_y_svs_3-size_interchr_SV_tip, curvature=0,  size=0.25, colour="green",alpha=.75) +
            geom_rect(data = data.frame(cov = 1, chr = gene_coord_now$chr),
                      xmin=gene_coord_now$start, xmax=gene_coord_now$end,
                      ymin=0, #max_y_svs_2 + (max_y_svs_3-max_y_svs_2)/2,
                      ymax=max_y_svs_3-size_interchr_SV_tip,
                      size=0.25, fill = "green", alpha=.6)
          
          
          # # add text annotation
          # dat_text <- data.frame(label = gene,
          #                        chr= factor(gene_coord_now$chr,levels = chr_selection$chr),
          #                        pos=gene_coord_now$start, y= max_y_svs_3+ 0.5
          # )
          dat_text <- data.frame(label = gene,
                                 chr= factor(gene_coord_now$chr,levels = chr_selection$chr),
                                 pos=(gene_coord_now$start+gene_coord_now$end)/2, y= max_y_svs_4
          )
          p = p + geom_point(data= dat_text, mapping = aes(x= pos, y = max_y_svs_3-size_interchr_SV_tip),
                             size=.5 ,colour="darkblue")
          p = p + geom_text(data = dat_text, mapping = aes(x= pos, y = y, label = label), 
                            size=size_gene_label ,fontface="italic") 
          # add point
         
        }
      }
    }
  }
  
  # add SV type description (add at the end so it is on top of lines)
  if(legend_SV_types){
  separation_labels= (scale_separation_SV_type_labels) * max_y
  dat_text <- data.frame(label = c("t2tINV (-/-)","h2hINV (+/+)","DUP (-/+)","DEL (+/-)"),
                         chr= factor(chr_selection$chr[1],levels = chr_selection$chr), # we plot on first chr
                         pos=rep(chr_selection$start[1] + pos_SVtype_description,4),
                         y= c(max_y_svs_2+separation_labels,
                              max_y_svs_2-separation_labels,max_y_svs_1+separation_labels,max_y_svs_1-separation_labels),
                         colour = c(colour_t2tINV, colour_h2hINV, colour_DUP, colour_DEL)
  )
  p = p + geom_text(data = dat_text, mapping = aes(x= pos, y = y, label = label), # col = colour),
                     size=1.5, fontface="bold", 
  colour = c(colour_t2tINV, colour_h2hINV, colour_DUP, colour_DEL)
  )


  # put copy number title further down:
  p = p + theme(axis.title.y = element_text(hjust=0.2))
  }

  # to remove vertical white lines:
  p = p + theme(panel.grid.minor = element_line(colour = NA))

  
  #SNV PLOTTING#
  
  if (! is.null(custom_annotation)){
    # get the maximum coordinates for the chr at hand
    chr.coords <- cnv.plot %>% 
      group_by(.data$chr) %>% 
      summarise(start = min(.data$start),
                end = max(.data$end)) %>% 
      ungroup()
    custom_annotation$selected = apply(custom_annotation, 1, function(row){
      sel = FALSE
      chr_i = row["chr"]
      pos_i = row["pos"]
      if (chr_i %in% cnv.plot$chr){
        max_coord = chr.coords %>% filter(.data$chr == chr_i) %>% pull(.data$end)
        min_coord = chr.coords %>% filter(.data$chr == chr_i) %>% pull(.data$start)
        if (dplyr::between(as.numeric(pos_i),min_coord, max_coord)){
          sel=TRUE
        }
      }
      return(sel)
    })
    custom_annotation_selected <- custom_annotation %>% filter(.data$selected)
    if (nrow(custom_annotation_selected) > 0){

      #Add anchor to have same plot area
      anchor.df <- data.frame()
      for (chr_i in unique(chr.coords$chr)){
        max_coord = chr.coords %>% filter(.data$chr == chr_i) %>% pull(.data$end)
        min_coord = chr.coords %>% filter(.data$chr == chr_i) %>% pull(.data$start)
        row.min = c("chr"=chr_i, "pos"=min_coord, "y"=0)
        row.max = c("chr"=chr_i, "pos"=max_coord, "y"=0)
        rows=bind_rows(row.min,row.max)
        anchor.df <- bind_rows(anchor.df, rows)
      }
      custom_annotation_selected <- custom_annotation_selected %>% 
        mutate(y = as.numeric(.data$y),
               pos = as.numeric(.data$pos))
      anchor.df <- anchor.df %>% 
        mutate(y = as.numeric(.data$y),
               pos = as.numeric(.data$pos))
    
      ann.plot <- ggplot(custom_annotation_selected, aes(x=pos, y=y)) +
        geom_point(col=ann_dot_col, size=ann_dot_size) +
        geom_point(data=anchor.df, col=NA, na.rm=T) +
        facet_grid(. ~ factor(chr, levels = unique(chr_selection$chr)), space="free_x", scales="free_x") +
        theme_classic() +
        scale_x_continuous(expand=c(0,0)) +
        ylab(ann_y_title) +
        # scale_fill_manual(limits=c("true","suspect", "false"), values=c("black", "gray47", "gray76")) +
        theme(
          text=element_text(size=size_text, colour="black"), #, face ="bold"),
          axis.text.x= element_blank(),
          axis.ticks.x= element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=size_text, colour="black"),
          axis.text.y=element_text(size=size_text, colour="black"),
          axis.line =  element_line(colour = "black", size=.25),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=.25),
          panel.grid.major.y=element_line( size=.01, color="gray"),
          strip.background = element_blank(),
          strip.placement = "outside",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.spacing.x = unit(.00625 * 3, "npc"),
          strip.text.x = element_blank(), ## to remove the title of the facet panels
          plot.title = element_blank(),
          legend.position = "none"
        )
      if (ann_one_scale) {
        ann.plot = ann.plot + scale_y_continuous(limits=c(0,1), breaks = seq(0,1,.25))
      }
    }
  }
  if (exists("ann.plot")){
    merged.p <- plot_grid(p, ann.plot, 
                          align = "hv", axis="trbl",
                          ncol=1, rel_heights = c(1, ann_rel_size))
  } else {
    merged.p <- p
  }
  return(merged.p)
}




