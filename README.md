# ReConPlot
R package to visualize complex genomic rearrangements by plotting copy number profiles and structural variants.

# Installation
R CMD INSTALL ReConPlot

# How to use and examples
Please see the detailed tutorial and documentation of the package for examples and best practices for using ReConPlot to generate publication-quality figures.

## Quick start
First load the package
```
library(ggplot2)
library(ReConPlot)
```
You will need three dataframes:
1. SV data (with columns chr1, pos1, chr2, pos2 and strands (+- notation)
2. CN data (with columns chr, start, end, copyNumber and minorAlleleCopyNumber)
3. Chromosome selection with genomic region(s) to plot (with columns chr, start, end)
Any extra column in the data frames will not be read. 

```
#SV data
print(head(sv_data))

sample chr1      pos1  chr2      pos2 filter homlen homseq inslen strands
Test chr1 203476803  chr1 204967317   PASS      .      .      .      -+
Test chr1 203476817  chr1 204967172   PASS      .      .      .      +-
Test chr1 203507815  chr9  37245409   PASS      .      .      .      +-
Test chr1 203509239 chr12  69072120   PASS      3    GTA      .      ++
Test chr1 203509728 chr12  69071885   PASS      .      .      .      -+
Test chr1 203509991  chr1 203522717   PASS      1      A      .      +-
```

```
#CN data
print(head(cn_data))

chr    start      end copyNumber minorAlleleCopyNumber
chr1        1  9631965          2                     1
chr1  9631966  9631966          3                     1
chr1  9631967 11239516          2                     1
chr1 11239517 11239533          3                     1
chr1 11239534 22578082          2                     1
chr1 22578083 27086500          2                     1
```
```
#Chromosome selection
chrs=c("chr1")
chr_selection = data.frame(
  chr=chrs,
  start=rep(0 ,length(chrs)),
  end=rep (250000000, length(chrs)) 
) 
print(chr_selection)

chr start     end
chr1     0 2.5e+08
```
You can then generate the plot with the ReconPlot function:
```
plot = ReConPlot(sv_data,
cn_data,
chr_selection=chr_selection,
legend_SV_types=T,
pos_SVtype_description=1000000,
scale_separation_SV_type_labels=1/23,
title="Example")

print(plot)
```
Please visit the [tutorial](Tutorial/tutorial.pdf) for advanced usage.



# Contact
If you have any comments or suggestions please raise an issue or contact us:\
Jose Espejo Valle-Inclan: jespejo@ebi.ac.uk\
Isidro Cortes-Ciriano: icortes@ebi.ac.uk
