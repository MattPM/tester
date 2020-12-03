# 10x analysis 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
library(dsb)
library(parallelDist)
library(here)

# set params 
project_title = "10X PBMC5k V3"
expected_cells = "~5,000"
res = 0.8
kparam = 40

# thresholds 
max_neg_logprotumi = 2.5
min_cell_logprotumi = 2.8
genemax = 3000
genemin = 80
mtthresh = 0.2

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/")
datapath = here("V2/10x_analysis/generated_data/")
dir.create(figpath); dir.create(datapath)

source("V2/functions/preprocessing_functions.R")
source("V2/functions/analysis_functions.R")

# read raw data
raw = Read10X_MPM("data/10x_data/10x_pbmc5k_V3/raw_feature_bc_matrix/")

# Separate RNA and protein matrix 
prot = raw[grep(rownames(raw), pattern = "Total"), ]
rna = raw[rownames(raw)[rownames(raw) %ni% rownames(prot)], ]

# calculate metadata 
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE)
pctmt = colSums(rna[mtgene, ])/colSums(rna)
umi  = colSums(rna)
log10umi = log10(umi)
umiprot = colSums(prot)
log10umiprot = log10(umiprot)
nGene = colSums(rna > 0)

# check to see if there are protein detected in drops with no RNA 
pdrop = names(umiprot)[names(umiprot) %ni% names(umi)]
stopifnot(length(pdrop) == 0)

# Confirm RNA assay reads in additional barcodes with no data 
rnadrop = names(umi)[names(umi) %ni% names(umiprot)]
qrna = quantile(umi[rnadrop])

# subset matrices 
int = intersect(names(umiprot), names(umi))
bcmax = length(int)
pctmt = pctmt[int]; umi = umi[int]; log10umi = log10umi[int]; nGene = nGene[int]

# combine into metadata 
md = as.data.frame(cbind(pctmt, umi, log10umi, nGene, umiprot, log10umiprot))

hist(md$log10umiprot[md$log10umiprot < 5], breaks = 100)

# define negative background subset protein matrix to empty drop subset
neg_drops2 = md %>% rownames_to_column("bc") %>% 
  filter(log10umiprot < max_neg_logprotumi & log10umiprot > 1.4)  %>% 
  filter(nGene < genemin) %$% bc

neg_prot2 = prot[ , neg_drops2] %>%  as.matrix()

# subset out outlier drops from positive protein matrix RNA based; add absolute ceiling for conservative cell estimate. 
positive_cells = md %>% rownames_to_column("bc") %>% 
  filter(log10umiprot > min_cell_logprotumi) %>% 
  filter(nGene < genemax & nGene > 200) %>% 
  filter(pctmt < mtthresh) %$%  
  bc

# subset protein matrix for cells 
pos_prot = prot[ , positive_cells] %>% as.matrix()


###########
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
pv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(pos_prot)) %>% mutate(class = "cell_containing")
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot2)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot2)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)

# save raw 
saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))

#############
# DSB normalize 
# define isotypes 
isotypes = rownames(pos_prot)[30:32]
# threshold 1 normalization
mtx2 = DSBNormalizeProtein(cell_protein_matrix = pos_prot,
                           empty_drop_matrix = neg_prot2,
                           denoise.counts = TRUE,
                           use.isotype.control = TRUE,
                           isotype.control.name.vec = isotypes)

# inspect distributions (plot in next script)
mgl = list(theme_bw(),  geom_point(size  = 0.3, alpha = 0.4) , geom_density_2d(color = 'red', size = 0.3))
dfplot = as.data.frame(t(mtx2))
ggplot(dfplot, aes(x = CD3_TotalSeqB, y = CD19_TotalSeqB  )) + mgl
############# 

# cluster 
rna_cells = rna[ ,positive_cells]
md_cells = md %>% rownames_to_column("bc") %>% filter(bc %in% positive_cells) %>% column_to_rownames("bc")
s = CreateSeuratObject(raw.data = rna_cells, min.cells = 40, min.genes = genemin, meta.data = md_cells)
s = SetAssayData(s, assay.type = "CITE", slot = "data",new.data = mtx2)

##### cluster 
prot = rownames(s@assay$CITE@data)
prot_subset = setdiff(prot, isotypes)

# Subset sd background normalized denoised protein 
s2_adt = GetAssayData(s, assay.type = "CITE", slot = "data")
s2_adt3 = s2_adt[prot_subset, ]
p3_dist = parDist(t(s2_adt3))
p3_dist = as.matrix(p3_dist)

# Protein clustering 
s = FindClusters(s, 
                 distance.matrix = p3_dist,
                 k.param = kparam,
                 print.output = F, 
                 resolution = res,
                 random.seed = 1,
                 algorithm = 3,
                 modularity.fxn = 1)
s = StashIdent(s, save.name = "clusters")

# run umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(s2_adt3), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(s@assay$CITE@data)))
saveRDS(df_dsb, file = paste0(datapath,project_title, "dsb_merged_result.RDS"))



# # split protein and rna data into separate matrix 
# prot = s@raw.data[grep(rownames(s@raw.data), pattern = "TotalSeq"), ]
# rna = s@raw.data[rownames(s@raw.data)[rownames(s@raw.data) %ni% rownames(prot)], ]
# dim(rna); dim(prot)
# stopifnot(dim(prot)[1] > 0 &  dim(rna)[1] > 0 )
# 
# # make new object with separate assays
# s = CreateSeuratObject(raw.data = rna)
# 
# # define neg cells 
# hist(log(s@meta.data$nUMI  + 1 ), breaks = 1000)
# hist(s@meta.data$nUMI[s@meta.data$nUMI < 100], breaks = 1000)
# ndrop = dim(s@meta.data)[1]
# 
# # Plot
# p1 = ggplot(s@meta.data, aes(x = log10(nUMI + 1 ) )) +
#   theme_bw()+ 
#   theme(text = element_text(size = 8)) + 
#   ggtitle(paste0(project_title, " raw_feature_bc_matrix: ", ndrop, " droplets")) + 
#   geom_density(fill = "#3e8ede") + 
#   geom_vline(xintercept = c(2.4, 1.5 ),   linetype = "dashed") + 
#   annotate("text", x = 0.5, y=2, label = " region 1: \n void of data ") + 
#   annotate("text", x = 2, y=2, label = " region 2: \n background drops \n define 'empty_drop_matrix' \n with these drops ") + 
#   annotate("text", x = 4, y=2, label = " region 3: \n cell containing droplets \n zomed in on next plot") 
# p1
# p2 = ggplot(s@meta.data %>% filter(log10(nUMI + 1) > 2.4), aes(x = log10(nUMI + 1 ) )) +
#   theme_bw() + 
#   theme(text = element_text(size = 8)) + 
#   geom_density(fill = "#3e8ede") + 
#   ggtitle(paste0(project_title, " drops containing cells "))  
# p3 = plot_grid( p1 , p2 ) 
# p3
# ggsave(p3, filename = paste0(figpath,project_title, "background.pdf"), width = 11, height = 4)
# 
# # add logumi to metadata for thresholding 
# md = s@meta.data %>%
#   rownames_to_column("bc") %>% 
#   mutate(log10umi = log10(nUMI + 1 )) %>%
#   select(bc, log10umi) %>% 
#   column_to_rownames("bc")
# s = AddMetaData(s,metadata = md)
# 
# # define negative drops based on threshold above 
# neg_drops2 = WhichCells(s, subset.name =  "log10umi",  accept.low = 1.5 , accept.high = 2.4)
# 
# # subset protein by negative drops   
# neg_prot2 = prot[ , neg_drops2 ]
# neg_prot2 = as.matrix(neg_prot2)
# 
# # subset out outlier drops from positive protein matrix 
# positive_cells = WhichCells(s, subset.name =  "log10umi", accept.low = 2.4, accept.high = 4.4)
# positive_prot = prot[ , positive_cells] 
# pos_prot = as.matrix(positive_prot)
# 
# # remove drops with no / minimal protein library based on protein 
# csp = log10(colSums(pos_prot))
# nop = csp[csp < 2.8]  %>% names( )
# pos_prot = pos_prot[ ,colnames(pos_prot) %ni% nop]
# 
# ### combined library plot 
# index1 = colnames(pos_prot)[1]; index2 = colnames(pos_prot)[length(colnames(pos_prot))]
# pos_tidy = pos_prot %>%
#   as.data.frame() %>% 
#   rownames_to_column("protein") %>% 
#   gather(cell, UMI, index1:index2 ) %>% 
#   group_by(cell) %>% 
#   summarize(log10_UMI = log10(sum(UMI))) %>% 
#   mutate(class = "cell_containing")
# index1 = colnames(neg_prot2)[1]; index2 = colnames(neg_prot2)[length(colnames(neg_prot2))]
# neg_tidy = neg_prot2 %>%
#   as.data.frame() %>% 
#   rownames_to_column("protein") %>% 
#   gather(cell, UMI, index1:index2 ) %>% 
#   group_by(cell) %>% 
#   summarize(log10_UMI = log10(sum(UMI))) %>% 
#   mutate(class = "background")
# ddf = rbind(pos_tidy, neg_tidy)
# # plot joint distributions 
# p = ggplot(ddf, aes(x = log10_UMI, fill = class, color = class )) +
#   theme_bw()+ 
#   ggsci::scale_fill_d3() + 
#   ggsci::scale_color_d3() + 
#   annotate("text", color = "darkorange", fontface = "bold", size = 3.7, x = 4, y=3, label = paste0("PRE-QC \n barcodes with cells: \n ", dim(pos_prot)[2], "\n expect loading: \n " ,expected_cells) ) + 
#   annotate("text", color = "dodgerblue3", fontface = "bold",  size = 3.7, x = 1, y=3, label = paste0("background drops: \n ", length(neg_drops2), " \n theretical max: \n", dim(raw)[2] ) ) + 
#   geom_histogram(aes(y=..density..), alpha=0.5, bins = 50,position="identity")+
#   geom_density(alpha = 0.5) + 
#   ylab("density") + xlab("log10 protein library size") + 
#   theme(axis.title.x = element_text(size = 14)) + 
#   ggtitle(project_title) + 
#   theme(plot.title = element_text(face = "bold",size = 20)) + 
#   theme(legend.position = c(0.2, 0.3), legend.margin = margin(0,0,0,0))
# p
# ggsave(p, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)
# p = p + xlab("") + ylab("")
# ggsave(p, filename = paste0(figpath,project_title, "protein_nolab_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)
# 
# # save raw 
# saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
# saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))
# saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
# saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))
# 
# 
# # DSB normalize 
# # define isotypes 
# isotypes = rownames(pos_prot)[30:32]
# # threshold 1 normalization
# mtx2 = DSBNormalizeProtein(cell_protein_matrix = pos_prot,
#                            empty_drop_matrix = neg_prot2,
#                            denoise.counts = TRUE,
#                            use.isotype.control = TRUE, 
#                            isotype.control.name.vec = isotypes)
# 
# # add protein data to Seurat object. 
# s2 = s %>% SubsetData(cells.use = colnames(mtx2), subset.raw = TRUE)
# pos_prot2 = pos_prot[ ,s2@cell.names]
# s2 = SetAssayData(s2, assay.type = "CITE", slot = "raw.data", new.data = pos_prot2)
# s2 = SetAssayData(s2, assay.type = "CITE", slot = "data",new.data = mtx2)
# 
# # ADD slot for CLR norm
# s2 = SetAssayData(s2, assay.type = "CITE_CLR", slot = "raw.data", new.data = pos_prot2)
# s2 = NormalizeData(s2,assay.type = "CITE_CLR", normalization.method = "genesCLR")
# 
# ##### cluster 
# 
# 
# # remove isotype controls
# prot = rownames(s2@assay$CITE@data)
# prot_subset = setdiff(prot, isotypes)
# 
# # Subset sd background normalized denoised protein 
# s2_adt = GetAssayData(s2, assay.type = "CITE", slot = "data")
# 
# # subset used markers 
# s2_adt3 = s2_adt[prot_subset, ]
# 
# #get distance matrix and cluster 
# p3_dist = parDist(t(s2_adt3))
# p3_dist = as.matrix(p3_dist)
# 
# # cluster based on surface protein
# s2 = FindClusters(s2, 
#                   distance.matrix = p3_dist,
#                   k.param = kparam,
#                   print.output = F, 
#                   resolution = res,
#                   random.seed = 1,
#                   algorithm = 3,
#                   modularity.fxn = 1)
# s2 = StashIdent(s2, save.name = "clusters")
# 
# 
# # run umap 
# library(reticulate); use_virtualenv("r-reticulate")
# library(umap)
# 
# # set umap config
# config = umap.defaults
# config$n_neighbors = 40
# config$min_dist = 0.4
# 
# # run umap
# ump = umap(t(s2_adt3), config = config)
# umap_res = ump$layout %>% as.data.frame() 
# colnames(umap_res) = c("UMAP_1", "UMAP_2")
# 
# # save results dataframe 
# df_dsb = cbind(s2@meta.data, umap_res, as.data.frame(t(s2@assay$CITE@data)))
# df_clr = cbind(s2@meta.data, umap_res, as.data.frame(t(s2@assay$CITE_CLR@data)))
# 
# # save_outputs 
# saveRDS(df_dsb, file = paste0(datapath,project_title, "dsb_merged_result.RDS"))
# saveRDS(df_clr, file = paste0(datapath,project_title, "clr_merged_result.RDS"))
# 
# 


