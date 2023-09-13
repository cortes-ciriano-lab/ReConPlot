  <!-- badges: start -->
  [![R-CMD-check](https://github.com/cortes-ciriano-lab/ReConPlot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cortes-ciriano-lab/ReConPlot/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# ReConPlot
R package to visualize complex genomic rearrangements by plotting copy number profiles and structural variants.

If you use ReConPlot please cite our [pre-print](https://www.biorxiv.org/content/10.1101/2023.02.24.529890v2)
```
ReConPlot – an R package for the visualization and interpretation of genomic rearrangements
Jose Espejo Valle-Inclán, Isidro Cortés-Ciriano
bioRxiv 2023.02.24.529890; doi: https://doi.org/10.1101/2023.02.24.529890
```

# Installation
You can clone this repository and install using the following command in the command line:
```
R CMD INSTALL ReConPlot/
```
Or use devtools to install directly from GitHub within R:
```
> devtools::install_github("cortes-ciriano-lab/ReConPlot")
```

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
p = ReConPlot(sv_data,
cn_data,
chr_selection=chr_selection,
legend_SV_types=T,
pos_SVtype_description=1000000,
scale_separation_SV_type_labels=1/23,
title="Example")
print(p)
```
The ReCon plots are ggplot objects and can be modified after generation as such, and they can be easily saved to a PDF file.
In our experience, the following dimensions work well for publication-quality figures and written reports. If using an annotation plot in combination with the ReCon plot, the height might need to be increased.

```
ggsave(filename = "example_ReConPlot.pdf", plot = p, width = 19, height = 5, units = "cm")
```

Please visit the [tutorial](Tutorial/tutorial.pdf) for advanced usage.

# Contact
If you have any comments or suggestions please raise an issue or contact us:\
Jose Espejo Valle-Inclan: jespejo@ebi.ac.uk\
Isidro Cortes-Ciriano: icortes@ebi.ac.uk
