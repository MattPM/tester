Normalizing and denoising protein expression data from droplet-based
single-cell profiling
================

### Analysis code to reproduce manuscript results

Mulè MP, Martins AJ, Tsang JS. Normalizing and denoising protein
expression data from droplet-based single cell profiling. bioRxiv.
2020;2020.02.24.963603.  
repository author: Matt Mulè <mulemp@nih.gov>  
general contact: John Tsang <john.tsang@nih.gov>

### The dsb R package

Analysis in this paper supports the dsb R package for CITE-seq data
denoising and normalization [Download the dsb package from the NIAID
github](https://github.com/niaid/dsb)

### to install the dsb package

install dsb directly from github:

``` r
require(devtools); devtools::install_github(repo = 'niaid/dsb')
library(dsb)
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

### Introduction

This analysis code pipeline is for reproducing the analysis in the
paper: **Normalizing and denoising protein expression data from droplet
based single cell profiling** The paper introduces the experiments and
statistical modeling underlying the dsb normalization method for
CITE-seq protein data. The dsb R package (see above) is a lightweight R
package and is the first dedicated method for normalizing and denoising
antibody based protein expression data from droplet based single cell
experiments (CITE-seq, REAP-seq, Mission Bio Tapestri etc). The method
was developed in [John Tsang’s
Lab](https://www.niaid.nih.gov/research/john-tsang-phd) by Matt Mulè,
Andrew Martins and John Tsang.

This analysis workflow presented in the preprint can be run on a laptop
with 16GB RAM. Note in this reproducible workflow, the directory
containing the dsb\_normalization.rproj file acts like the root
directory by using the ‘here’ R package. To reproduce analysis results,
simply download the github repository, and add data described below into
a folder called data/V2\_Data from any folder from which you run the
analysis. Each R script is a self-contained analysis script reading the
data from the /data folder and writing figures or results files within a
generated\_data directory within each analyis subdirectory as outlined
below.

### Software package versions used in this analysis

Please see full session info at the end of this script

R version 3.5.3 mclust\_5.4.5  
reticulate\_1.12  
umap\_0.2.3.1  
magrittr\_1.5  
forcats\_0.4.0  
stringr\_1.4.0  
dplyr\_0.8.5  
purrr\_0.3.3  
readr\_1.3.1  
tidyr\_1.0.2  
tibble\_2.1.1  
tidyverse\_1.2.1  
Seurat\_2.3.4  
Matrix\_1.2-15  
cowplot\_0.9.4  
ggplot2\_3.1.1  
here\_0.1

### R packages used in this this analysis

``` r
# modeling  / analysis / helper functions  
library(Seurat)
library(tidyverse)
library(magrittr)
library(here)
library(mclust)
library(dsb) # see above for installation instructions.

# visualization
library(ggrepel)
library(ggridges)
library(pals)

# For umap only, must have python virtual env for this particular pipeline, for example run: 
# python must be installed 
virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "umap-learn")
use_virtualenv("r-reticulate")
library(umap)
library(reticulate)
library(umap)
```

The data directory: *data/V2\_Data/* contains the starting raw data
needed for analysis. Prior to running analysis confirm the data are in
the correct repository.

``` r
suppressMessages(library(here))
# all STARTING data 
data_ = c(
  "CITEseq_raw_PMID32094927_seurat2.4.rds", 
  "unstained_control_singlets.rds"
)
stopifnot(data_ %in% list.files(here("data/V2_Data/")))

# background drop data 
data_ = c(
"adt_neg_dmx_list.rds",
"adt_neg_dmx.rds",
"adt_neg_full_list.rds",
"adt_neg_full_qc.rds"
)
stopifnot(data_ %in% list.files(here("data/V2_Data/background_data/")))

# pipeline scripts scripts 
scripts_ = c(
  "protein_annotation_functions.r",
  "geneplot4.R"
)
stopifnot(scripts_ %in% list.files(here("V2/functions/")))
```

### dsb normalize PBMC data from 20 individuals.

V2/dsb\_normalize\_cluster\_pipeline/

``` r
source(here("V2/dsb_normalize_cluster_pipeline/1_dsb_normalize.r"))
# creates data drectory here("V2/dsb_normalize_cluster_pipeline/generated_data/") containing results 
# "dsb2_normlog10_list_bybatch.rds"
# "h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"
```

### run umap

reads data
from:  
V2/dsb\_normalize\_cluster\_pipeline/generated\_data/h1\_d0\_singlets\_ann\_Seurat2.4\_dsbnorm.rds

``` r
source(here("V2/dsb_normalize_cluster_pipeline/2_run_umap.r"))
```

### Figure generation

``` r
source(here("V2/dsb_normalize_cluster_pipeline/3_figure_generation.r"))
source(here("V2/dsb_normalize_cluster_pipeline/4_manual_gate_plots.r"))
```

### dsb process

This section illustrates the dsb process in steps and assesses
underlying modeling assumptions

``` r
source(here("V2/dsb_process_plots/6_mean_isotype_v_mean_control.R"))
source(here("V2/dsb_process_plots/6a_isotype_figure_generation.r"))
source(here("V2/dsb_process_plots/7_neg_control_plots.R"))
source(here("V2/dsb_process_plots/8_mixture_fits.r"))
```

### batch norm analysis

Assessment of single vs multi batch normalization and sensitivity of
each normalization scheme to defining background with hashing or library
size distributions.

``` r
source(here("V2/parameter_sensitivity/empty_drop_threshold_batch.r"))
```

### 10X data analysis

Read Cell Ranger raw output, select negative drops, run DSB
normalization, run each normalization modeling step separately, cluster
cells and plot distributions across clusters and on biaxial gates. Test
underlying modeling assumptions for each external dataset.

``` r

# 10K v3 data 
source(here("V2/10x_analysis/10x_pbmc_10k_V3.r"))
source(here("V2/10x_analysis/10x_pbmc_10k_V3_figure_generation.r"))
# 5k V3 data 
source(here("V2/10x_analysis/10x_pbmc_5k_V3.r"))
source(here("V2/10x_analysis/10x_pbmc_5k_V3_figure_generation.r"))
# Next Gem data 
source(here("V2/10x_analysis/10x_pbmc_NextGem.r"))
source(here("V2/10x_analysis/10x_pbmc_NextGem_figure_generation.r"))
# 5 prime data 
source(here("V2/10x_analysis/10x_pbmc_5prime_5k.r"))
source(here("V2/10x_analysis/10x_pbmc_5prime_5k_figure_generation.r"))
```

# Sessioninfo

R version 3.5.3

attached base packages: stats graphics grDevices utils datasets methods
base

other attached packages: mclust\_5.4.5 reticulate\_1.12 umap\_0.2.3.1
magrittr\_1.5 forcats\_0.4.0 stringr\_1.4.0 dplyr\_0.8.5 purrr\_0.3.3
readr\_1.3.1 tidyr\_1.0.2 tibble\_2.1.1 tidyverse\_1.2.1 Seurat\_2.3.4
Matrix\_1.2-15 cowplot\_0.9.4 ggplot2\_3.1.1 here\_0.1

loaded via a namespace (and not attached): readxl\_1.3.1 snow\_0.4-3
backports\_1.1.4 Hmisc\_4.2-0 plyr\_1.8.4 igraph\_1.2.4.1
lazyeval\_0.2.2 splines\_3.5.3 inline\_0.3.15 digest\_0.6.25
foreach\_1.4.4 htmltools\_0.3.6 lars\_1.2 rsconnect\_0.8.16
gdata\_2.18.0 checkmate\_1.9.3 cluster\_2.0.7-1 mixtools\_1.1.0
ROCR\_1.0-7 modelr\_0.1.4 matrixStats\_0.54.0 R.utils\_2.8.0
askpass\_1.1 prettyunits\_1.0.2 colorspace\_1.4-1 rvest\_0.3.4
haven\_2.1.0 xfun\_0.7 callr\_3.2.0 crayon\_1.3.4 jsonlite\_1.6
survival\_2.43-3 zoo\_1.8-6 iterators\_1.0.10 ape\_5.3 glue\_1.3.1
gtable\_0.3.0 pkgbuild\_1.0.3 kernlab\_0.9-27 rstan\_2.19.3
prabclus\_2.3-1 DEoptimR\_1.0-8 scales\_1.0.0 mvtnorm\_1.0-10
bibtex\_0.4.2 Rcpp\_1.0.1 metap\_1.1 dtw\_1.20-1 htmlTable\_1.13.1
foreign\_0.8-71 bit\_1.1-14 proxy\_0.4-23 SDMTools\_1.1-221.1
Formula\_1.2-3 stats4\_3.5.3 tsne\_0.1-3 StanHeaders\_2.21.0-1
htmlwidgets\_1.3 httr\_1.4.0 gplots\_3.0.1.1 RColorBrewer\_1.1-2
fpc\_2.2-1 acepack\_1.4.1 modeltools\_0.2-22 ica\_1.0-2 pkgconfig\_2.0.2
loo\_2.3.1 R.methodsS3\_1.7.1 flexmix\_2.3-15 nnet\_7.3-12
tidyselect\_0.2.5 rlang\_0.4.5 reshape2\_1.4.3 cellranger\_1.1.0
munsell\_0.5.0 tools\_3.5.3 cli\_1.1.0 generics\_0.0.2 broom\_0.5.2
ggridges\_0.5.1 evaluate\_0.14 yaml\_2.2.0 npsurv\_0.4-0 processx\_3.3.1
knitr\_1.23 bit64\_0.9-7 fitdistrplus\_1.0-14 robustbase\_0.93-5
caTools\_1.17.1.2 RANN\_2.6.1 packrat\_0.5.0 pbapply\_1.4-0
nlme\_3.1-137 R.oo\_1.22.0 xml2\_1.2.0 hdf5r\_1.2.0 compiler\_3.5.3
rstudioapi\_0.10 png\_0.1-7 lsei\_1.2-0 stringi\_1.4.3 ps\_1.3.0
lattice\_0.20-38 vctrs\_0.2.4 pillar\_1.4.1 lifecycle\_0.1.0
Rdpack\_0.11-0 lmtest\_0.9-37 data.table\_1.12.2 bitops\_1.0-6
irlba\_2.3.3 gbRd\_0.4-11 R6\_2.4.0 latticeExtra\_0.6-28
KernSmooth\_2.23-15 gridExtra\_2.3 codetools\_0.2-16 MASS\_7.3-51.1
gtools\_3.8.1 assertthat\_0.2.1 openssl\_1.4 rprojroot\_1.3-2
withr\_2.1.2 diptest\_0.75-7 parallel\_3.5.3 doSNOW\_1.0.16 hms\_0.4.2
grid\_3.5.3 rpart\_4.1-13 class\_7.3-15 rmarkdown\_1.13
segmented\_0.5-4.0 Rtsne\_0.15 lubridate\_1.7.4 base64enc\_0.1-3

### data initialization script.

The code below was run outside of this workflow for convenience, it is
shown here for reference. This step reformats the dataset hosted at the
[figshare
repository](https://nih.figshare.com/collections/Data_and_software_code_repository_for_Broad_immune_activation_underlies_shared_set_point_signatures_for_vaccine_responsiveness_in_healthy_individuals_and_disease_activity_in_patients_with_lupus_Kotliarov_Y_Sparks_R_et_al_Nat_Med_DOI_https_d/4753772)
that is associated with the
[manuscript](https://doi.org/10.1038/s41591-020-0769-8) where we
initially reported this dataset. The script below removes the protein
normalized data slot which used an earlier version of the DSB packge for
normaliation. In addition normaized RNA data, metadata, clustering
graphs tSNE results are removed to reduce object size and, cell type
annotations as reported in our (link above) are added to the object as
metadata.

``` r
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

datapath = here("data/V2_Data/")

# clear unneeded data slots tsne cell embeddings, unneeded metadata 
h1 = readRDS(file = "data/V2_Data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds")
h1@dr$protpca = NULL
h1@dr$tsne_p10 = NULL
h1@dr$tsne_p50 = NULL
h1@dr$tsne_p90 = NULL
h1@dr$tsne_p130 = NULL
h1@dr$tsne_p170 = NULL

h1@assay$HTO = NULL
h1@assay$CITE@data = NULL
h1@assay$CITE@scale.data = NULL
h1@data = NULL
# remove unused cells RNA data from raw slot
h1 = h1 %>% SubsetData(subset.raw = TRUE)

# remove unused additional metadata 
mds = read_delim(file = here("init_data_ignore/md_strip.txt"),delim = "\t")
vrm = mds$vars_
for (i in 1:length(vrm)) {
  vr = vrm[i]
  h1@meta.data[[vr]] = NULL
}

# add cluster labels from Fig 6 of https://www.nature.com/articles/s41591-020-0769-8 for reference in dsb paper.
cmd = read_delim(file = paste0(datapath, "clustree_node_labels_withCellTypeLabels.txt"), delim = '\t')
cmd1 = cmd %>% filter(clustering == "p3_dist_1")
cmd3 = cmd %>% filter(clustering == "p3_dist_3")

mdn = h1@meta.data %>% 
  mutate(celltype_label_1 = plyr::mapvalues(x = h1@meta.data$p3_dist_1, from = cmd1$cluster, to =cmd1$`Cell Type label`)) %>% 
  mutate(celltype_label_3 = plyr::mapvalues(x = h1@meta.data$p3_dist_3, from = cmd3$cluster, to =cmd3$`Cell Type label`)) %>% 
  select(barcode_check, celltype_label_3, celltype_label_1) %>% 
  column_to_rownames("barcode_check")

h1 = h1 %>% AddMetaData(metadata = mdn)

# save this starting object for DSB paper analysis  
saveRDS(h1, file = paste0(datapath, "CITEseq_raw_PMID32094927_seurat2.4.rds"))
```

### Project directory structure

├── README.Rmd  
├── README.md  
├── V2  
│ ├── 10x\_analysis  
│ │ ├── 10x\_pbmc\_10k\_V3.r  
│ │ ├── 10x\_pbmc\_10k\_V3\_figure\_generation.r  
│ │ ├── 10x\_pbmc\_5k\_V3.r  
│ │ ├── 10x\_pbmc\_5k\_V3\_figure\_generation.r  
│ │ ├── 10x\_pbmc\_5prime\_5k.r  
│ │ ├── 10x\_pbmc\_5prime\_5k\_figure\_generation.r  
│ │ ├── 10x\_pbmc\_NextGem.r  
│ │ ├── 10x\_pbmc\_NextGem\_figure\_generation.r  
│ ├── dsb\_normalize\_cluster\_pipeline  
│ │ ├── 1\_dsb\_normalize.r  
│ │ ├── 2\_run\_umap.r  
│ │ ├── 3\_figure\_generation.r  
│ │ ├── 4\_manual\_gate\_plots.r  
│ │ ├── generated\_data  
│ │ │ ├── MU1Denoised\_dsb\_Mtx.rds  
│ │ │ ├── NonDenoised\_dsb\_Mtx.rds  
│ │ │ ├── dsb2\_normlog10\_list\_bybatch.rds  
│ │ │ ├── h1\_d0\_singlets\_ann\_Seurat2.4\_dsbnorm.rds  
│ │ │ └── h1\_sng\_metadata\_umapdim\_prots\_dataframe.rds  
│ │ └── manual\_gates  
│ │ ├── clr\_gates.r  
│ │ └── dsb\_gates.r  
│ ├── dsb\_process\_plots  
│ │ ├── 6\_mean\_isotype\_v\_mean\_control.R  
│ │ ├── 6a\_isotype\_figure\_generation.r  
│ │ ├── 7\_neg\_control\_plots.R  
│ │ ├── 8\_mixture\_fits.r  
│ │ └── generated\_data  
│ │ ├── b1\_dsb\_norm\_adt\_mtx.rds  
│ │ ├── b2\_dsb\_norm\_adt\_mtx.rds  
│ │ ├── dsb\_norm\_adt\_mtx.rds  
│ │ ├── isotype\_values\_dsb.rds  
│ │ ├── logstats.rds  
│ │ ├── mixture\_model\_metadata\_merged.rds  
│ │ ├── multi\_component\_model\_comparison.rds  
│ │ └── normadt.rds  
│ ├── functions  
│ │ ├── analysis\_functions.R  
│ │ ├── geneplot4.R  
│ │ └── protein\_annotation\_functions.r  
│ └── parameter\_sensitivity  
│ ├── empty\_drop\_threshold\_batch.r  
│ └── generated\_data  
├── data  
│ ├── 10x\_data  
│ │ ├── 10x\_pbmc10k\_V3  
│ │ │ ├── pbmc\_10k\_protein\_v3 — Cell Ranger.htm  
│ │ │ └── raw\_feature\_bc\_matrix  
│ │ │ ├── barcodes.tsv.gz  
│ │ │ ├── features.tsv.gz  
│ │ │ └── matrix.mtx.gz  
│ │ ├── 10x\_pbmc5k\_NextGem  
│ │ │ ├── 5k\_pbmc\_protein\_v3\_nextgem — Cell Ranger.htm  
│ │ │ └── raw\_feature\_bc\_matrix  
│ │ │ ├── barcodes.tsv.gz  
│ │ │ ├── features.tsv.gz  
│ │ │ └── matrix.mtx.gz  
│ │ ├── 10x\_pbmc5k\_V3  
│ │ │ ├── 5k\_pbmc\_protein\_v3 - 5k Peripheral blood mononuclear cells
(PBMCs) from a healthy donor.html  
│ │ │ ├── 5k\_pbmc\_protein\_v3\_molecule\_info.h5  
│ │ │ └── raw\_feature\_bc\_matrix  
│ │ │ ├── barcodes.tsv.gz  
│ │ │ ├── features.tsv.gz  
│ │ │ └── matrix.mtx.gz  
│ │ └── 10x\_pbmc\_5prime\_5k  
│ │ ├── raw\_feature\_bc\_matrix  
│ │ │ ├── barcodes.tsv.gz  
│ │ │ ├── features.tsv.gz  
│ │ │ └── matrix.mtx.gz  
│ │ └── vdj\_v1\_hs\_pbmc2\_5gex\_protein — Cell Ranger.htm  
│ ├── V2\_Data  
│ │ ├── CITEseq\_raw\_PMID32094927\_seurat2.4.rds  
│ │ ├──
H1\_day0\_scranNorm\_adtbatchNorm\_dist\_clustered\_TSNE\_labels.rds  
│ │ ├── background\_data  
│ │ │ ├── adt\_neg\_dmx.rds  
│ │ │ ├── adt\_neg\_dmx\_list.rds  
│ │ │ ├── adt\_neg\_full\_list.rds  
│ │ │ └── adt\_neg\_full\_qc.rds  
├── dsb\_normalization.Rproj

### Project options

make rmd readme used: usethis::use\_readme\_rmd()

### notes

A review of this code has been conducted, no critical errors exist, and
to the best of the authors knowledge, there are no problematic file
paths, no local system configuration details, and no passwords or keys
included in this code. For questions about the dsb software package,
please open an issue at the [dsb github
repository](https://github.com/niaid/dsb).

Primary author(s): Matt Mulè  
Organizational contact information: General: john.tsang AT nih.gov,
code: mulemp AT nih.gov \[permanent address mattmule AT gmail\] Date of
release: Oct 7 2020  
Version: NA  
License details: NA  
Description: code to reproduce analysis of manuscript  
Usage instructions: Provided in this markdown  
Example(s) of usage: NA  
Proper attribution to others, when applicable: NA

### code check

``` r

# checked repo for PII and searched for strings with paths 

# code check 
library(lintr)
fcn = suppressMessages(list.files(here("functions"), ".r", full.names = TRUE))
tenx = suppressMessages(list.files(here("V2/10x_analysis/"),".r", full.names = TRUE))
pipel = suppressMessages(list.files(here("V2/dsb_normalize_cluster_pipeline/"),pattern = ".r", full.names = TRUE))
process = suppressMessages(list.files(here("V2/dsb_process_plots/"),pattern = ".r", full.names = TRUE))
param = suppressMessages(list.files(here("V2/parameter_sensitivity/"),pattern = ".r", full.names = TRUE))

# code check
scp = c(fcn, tenx, pipel, process, param) %>% as.list()
lt = suppressMessages((  
  lapply(scp, lintr::lint)
  )


# ignoring warnings on stylistic choices that do not impact code: 
# use of = not the <- operator for assignment 
# trailing whitespace esp. after "+" in ggplot commands 
# use of space after commas e.g. df[1, ] not df[1,]
# use of i %in% 1:length() in for loops instead of seq_along
```
