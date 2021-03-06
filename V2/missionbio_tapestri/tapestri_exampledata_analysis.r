suppressMessages(library(tidyverse))
suppressMessages(library(dsb))
suppressMessages(library(Seurat))
suppressMessages(library(here))

# proj title for fig paths 
project_title = "missionbio tapestri"

## change paths !! 
figpath = here("V2/missionbio_tapestri/figures/"); dir.create(figpath)
datapath = here("V2/missionbio_tapestri/generated_data/"); dir.create(datapath)
test = read.delim(file = here("git_ignore/AML-4-cell-line-multiomics-adt-counts.tsv"), sep = "\t",header = T)

# make dataframe
test = as.data.frame(test)
prot = test %>% spread(ab_description, raw) 
prot[is.na(prot)] <-  0 

# transpose into cells x prot matrix 
prot = prot %>% 
  column_to_rownames("cell_barcode") %>% 
  t %>% 
  as.data.frame()

# calculate library size of droplets to make rough thresholds for cell containing and ambient droplets 
prot_size = log10(Matrix::colSums(prot, na.rm = TRUE))
md = as.data.frame(prot_size)
md$bc = colnames(prot)
hist(md$prot_size, breaks = 100)

# define a vector of background / empty droplet barcodes based on protein library size
background_drops = md[md$prot_size < 2.5 & md$prot_size > 1.4, ]$bc
negative_mtx_rawprot = prot[ , background_drops] %>%  as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size 
positive_cells = md[md$prot_size > 2.7, ]$bc
cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()

# no isotype data available 
# normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot,
  empty_drop_matrix = negative_mtx_rawprot,
  denoise.counts = FALSE,
  use.isotype.control = FALSE
)




##########################
# visualization of drop distribution 
plot_layer = list(theme_bw() , 
                  ggsci::scale_fill_d3(), ggsci::scale_color_d3() ,
                  geom_histogram(aes(y=..count..), alpha=0.5, bins = 50,position="identity"),
                  geom_density(alpha = 0.5), 
                  ylab("Number of Drops"),  xlab("log10 protein library size"), 
                  theme(axis.title.x = element_text(size = 14)),
                  theme(plot.title = element_text(face = "bold",size = 14)),
                  theme(legend.position = c(0.8, 0.7), legend.margin = margin(0,0,0,0))
)
pv = md  %>% filter(bc %in% colnames(cells_mtx_rawprot)) %>% mutate(class = "cell_containing")
nv = md  %>% filter(bc %in% colnames(negative_mtx_rawprot)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = prot_size, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", nrow(test), "\n", 
    "estimated cell containing drops  = ", ncol(cells_mtx_rawprot), "\n",
    "negative droplets = ", ncol(negative_mtx_rawprot)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = prot_size, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)



##########################
# clustering analysis 

# calculate distance matrix for clustering 
p_dist = dist(t(dsb_norm_prot))
p_dist = as.matrix(p_dist)

# Graph based clustering 
s = CreateSeuratObject(raw.data =  cells_mtx_rawprot, min.cells = 0, min.genes = 0)
s = FindClusters(s, resolution = 0.5, distance.matrix = p_dist)

# add CLR norm 
s = SetAssayData(object = s, assay.type = "CLR", slot = "raw.data",new.data = cells_mtx_rawprot)
s = NormalizeData(s, assay.type = "CLR", normalization.method = "genesCLR")

# heatmap of protein values 
prots = rownames(s@raw.data)
adt_data = cbind(s@meta.data, as.data.frame(t(dsb_norm_prot)))
adt_plot = adt_data %>% 
  group_by(res.0.5) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("res.0.5") %>% 
  t %>% 
  as.data.frame

x = pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"),
                   filename = paste0(figpath,project_title, "_tapestri_heatmap.pdf"),
                   fontsize_row = 8, border_color = NA, width = 4, height = 4)


# umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(dsb_norm_prot), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(dsb_norm_prot)))
df_clr = cbind(s@meta.data, umap_res, as.data.frame(t(s@assay$CLR@data)))
saveRDS(df_dsb, file = paste0(datapath, project_title,"df_dsb.rds"))
saveRDS(df_clr, file = paste0(datapath, project_title,"df_clr.rds"))

#######################
### Visualizations 

# clusters by umap dims 
p = ggplot(df_dsb, aes(x = UMAP_1, y= UMAP_2, color = res.0.5 )) + theme_bw() + 
  geom_point(shape = 16) + 
  ggsci::scale_color_d3(alpha = 0.7)
ggsave(p, filename = paste0(figpath,"tapestri_umap_dsb.pdf"), width = 4, height = 4)

# long format for DSB heatmap 
index1 = colnames(df_dsb)[7]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# plots 
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2, color = DSB)) +
  geom_point(shape = 16) + 
  theme(legend.key.size = unit(0.8, units = "cm"), 
        legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 20, face = "bold")) + 
  scale_color_viridis_c(option = "B") +
  facet_wrap(~ prot,nrow = 2) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15, face = "bold")) + 
  theme(legend.position = "right")
ggsave(p, filename = paste0(figpath,  "clusters_DSBdist.png"), width = 12, height = 7.4)

# scatter plot flow 
mg = list(
  geom_point(shape = 16),
  geom_density_2d(), 
  geom_vline(xintercept = 0, color = "red3"),
  geom_hline(yintercept = 0, color = "red3") 
  )

p1 = ggplot(df_dsb, aes(x = CD19_0D, CD3_5A)) + mg
p2 = ggplot(df_dsb, aes(x = CD34_9C, CD30_6E)) + mg
p = plot_grid(p1, p2)
ggsave(p, filename = paste0(figpath, "scatter_prots_tapestri.pdf"),width = 6, height = 3 )


# vln plot 
dsb$prot = factor(dsb$prot, levels = x$tree_row$labels[x$tree_row$order])
dsb$res.0.5 = factor(dsb$res.0.5, levels = x$tree_col$labels[x$tree_col$order])
p = ggplot(dsb, aes(x = prot, y = DSB, fill = prot, color = prot)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 10)) + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) +
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, "vln_prots_tapestri.pdf"),width = 10, height = 4)


#### add CLR plot 
# save results dataframe 


# long format for CLR heatmap 
index1 = colnames(df_dsb)[7]; index2 = colnames(df_dsb)[ncol(df_dsb)]
clr = df_clr %>% gather(prot, CLR, index1:index2)
clr$prot = factor(dsb$prot, levels = x$tree_row$labels[x$tree_row$order])
clr$res.0.5 = factor(dsb$res.0.5, levels = x$tree_col$labels[x$tree_col$order])
p = ggplot(clr, aes(x = prot, y = CLR, fill = prot, color = prot)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 10)) + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) +
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  coord_flip() 
ggsave(p, filename = paste0(figpath, "CLR_vln_prots_tapestri.pdf"),width = 10, height = 4)
p = ggplot(clr, aes(x = prot, y = CLR, fill = prot)) + 
  theme_bw() +
  ggsci::scale_fill_d3() + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) + 
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1) +
  coord_flip()
p

ggsave(p, filename = paste0(figpath, "CLRvln_prots_tapestri.pdf"),width = 6, height = 3 )



# > sessionInfo()
# R version 3.5.3 Patched (2019-03-11 r77192)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] here_0.1        Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   dsb_0.1.0       forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5    
# [9] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    ggridges_0.5.1      class_7.3-15        modeltools_0.2-22   rprojroot_1.3-2    
# [7] mclust_5.4.5        htmlTable_1.13.1    base64enc_0.1-3     proxy_0.4-23        rstudioapi_0.10     npsurv_0.4-0       
# [13] bit64_0.9-7         flexmix_2.3-15      mvtnorm_1.0-10      lubridate_1.7.4     xml2_1.2.0          codetools_0.2-16   
# [19] splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3      
# [25] jsonlite_1.6        packrat_0.5.0       ica_1.0-2           broom_0.5.2         cluster_2.0.7-1     kernlab_0.9-27     
# [31] png_0.1-7           R.oo_1.22.0         compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1   
# [37] lazyeval_0.2.2      limma_3.38.3        cli_1.1.0           lars_1.2            acepack_1.4.1       htmltools_0.3.6    
# [43] tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          reshape2_1.4.3      RANN_2.6.1         
# [49] Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         ape_5.3             gdata_2.18.0        nlme_3.1-137       
# [55] iterators_1.0.10    fpc_2.2-1           lmtest_0.9-37       gbRd_0.4-11         xfun_0.7            rvest_0.3.4        
# [61] irlba_2.3.3         lifecycle_0.1.0     gtools_3.8.1        DEoptimR_1.0-8      zoo_1.8-6           MASS_7.3-51.1      
# [67] scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12    
# [73] pbapply_1.4-0       gridExtra_2.3       segmented_0.5-4.0   rpart_4.1-13        latticeExtra_0.6-28 stringi_1.4.3      
# [79] foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        dtw_1.20-1          Rdpack_0.11-0      
# [85] SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38    
# [91] ROCR_1.0-7          htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          magrittr_1.5       
# [97] R6_2.4.0            snow_0.4-3          gplots_3.0.1.1      generics_0.0.2      Hmisc_4.2-0         pillar_1.4.1       
# [103] haven_2.1.0         foreign_0.8-71      withr_2.1.2         mixtools_1.1.0      fitdistrplus_1.0-14 survival_2.43-3    
# [109] nnet_7.3-12         tsne_0.1-3          hdf5r_1.2.0         modelr_0.1.4        crayon_1.3.4        KernSmooth_2.23-15 
# [115] grid_3.5.3          readxl_1.3.1        data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7     
# [121] R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0  
# 
# 
# 
