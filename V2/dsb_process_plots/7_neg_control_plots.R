# Version 2 with added droplets from preprocessing - uncounted drops below the 35,000 rank ADT threshold. 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here)) 


# savepaths 
figpath = here("V2/dsb_process_plots/figures/")
datapath = here("V2/dsb_process_plots/generated_data/")

# define a random subset of proteins to label 
plabel = c("CD56_PROT", "CD294_PROT", "CD40_PROT", "HLA-DR_PROT", "CD278_PROT", 
           "CD197_PROT", "CD25_PROT","CD62L_PROT", "CD28_PROT", "CD5_PROT", 
           "CD45RA_PROT", "CD8_PROT", "HLA-ABC_PROT", "CD18_PROT", "IgM_PROT", 
           "CD244_PROT", "CD16_PROT", "CD38_PROT", "TCRgd_PROT", "IgD_PROT")


#  load unstained control cell data 
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))

# summarize unstained control data by protein stats. 
uadt = 
  un@assay$CITE@raw.data %>% 
  as.matrix() %>% t %>%
  as.data.frame() %>% 
  rownames_to_column("barcode_check") %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_unstained_control = median(count),
            mean_unstained_control = mean(count), 
            var_unstained_control = var(count)
            ) 


# full ADT output without hashing based subsetting of negative drops. 
adt_neg_full_qc = readRDS(file = here("data/V2_Data/background_data/adt_neg_full_qc.rds"))

#  Empty / background droplets summarize data by protein stats. 
adt_neg_full_summary =
  adt_neg_full_qc %>% 
  rename(barcode_check = bc) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_negative_droplet = median(count),
            mean_negative_droplet = mean(count),
            var_negative_droplet = var(count)
            ) 

# cells summary stats 
raw = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")
cells_logadt = raw@assay$CITE@raw.data %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode_check") %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_cells = median(count),
            mean_cells = mean(count), 
            var_cells = var(count)
  ) 

# add gaussian mixture mu 1 across cells for supplemental figure on aternate background estimation
adtraw = raw@assay$CITE@raw.data
adtlog = log10(adtraw + 1)
library(mclust)
prot_model = apply(adtlog, 1, function(x) { mclust::Mclust(x, G=2, warn = TRUE , verbose = TRUE) } )

# extract mean and variance for each protein within each gaussian population 
# mean 1 
p_mean1 = lapply(prot_model, function(x){x$parameters$mean[1]}) %>% unlist()
p_var1 = lapply(prot_model, function(x){x$parameters$variance$sigmasq[1]}) %>% unlist() 
# mean 2 
p_mean2 = lapply(prot_model, function(x){x$parameters$mean[2]}) %>% unlist()
p_var2 = lapply(prot_model, function(x){x$parameters$variance$sigmasq[2]}) %>% unlist() 
# merge 
log_stats = as.data.frame(cbind(p_var1, p_var2, p_mean1, p_mean2)) %>% rownames_to_column("prot")
saveRDS(log_stats,file = paste0(datapath, "logstats.rds"))
log_stats = readRDS(here("V2/dsb_process_plots/generated_data/logstats.rds"))
# merge data by protein.   
merged = full_join(adt_neg_full_summary, uadt, by = "prot")
merged = full_join(merged, cells_logadt, by = "prot")
merged = full_join(merged, log_stats, by = "prot")

## plot background mean by unstained and empty 
corp = list(geom_point(shape = 21, size = 2.8, fill = "red3"), 
  geom_abline(slope = 1, linetype = "dashed"),  
  geom_smooth(method = "lm", color = "black"), 
  ggpubr::stat_cor(method = "pearson"), 
  theme_bw()
)

# mean A vs negative droplet mean 
p1 = ggplot(merged, aes(x = p_mean1, y = mean_negative_droplet)) + 
  corp + 
  xlab("mean of background: stained cells") + 
  ylab(" Empty droplet mean log10 + 1 protein") + 
  ggtitle("negative stained cell mean A vs empty droplets") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
  aes(label = prot), size = 2.5,force = TRUE,fontface = "bold",  segment.size = 0.2, box.padding = 0.5) 

# mean A vs unstained control 
p2 = ggplot(merged, aes(x = p_mean1, y = mean_unstained_control))  +
  corp + 
  xlab("mean of background: stained cells") + 
  ylab("unstained controls mean log10 + 1 protein") + 
  ggtitle("negative stained cell mean A vs unstained control") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
  aes(label = prot), size = 2.5,force = TRUE,fontface = "bold",  segment.size = 0.2, box.padding = 0.5) 
pg = plot_grid(p2,p1, ncol = 1)
ggsave(pg,filename = paste0(figpath, "mean1_cellestimate_vs_neg_unstained.pdf"), height = 8.5, width = 4.5)


# mean in unstained control cells vs mean in empty droplets. 
ndropsfull = ndrops = adt_neg_full_qc$bc %>% unique %>% length()

# mean unstained vs mean empty defined by protein library size distribution for main Fig 1B 
corp2 = list(geom_point(shape = 21, size = 2.8, fill = "red3"), 
             geom_smooth(method = "lm", color = "black"), 
             ggpubr::stat_cor(method = "pearson"), theme_bw()
)
pmain = ggplot(merged, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp2 + 
  ylab("unstained controls mean log10 + 1 protein") + 
  xlab(" Empty droplet mean log10 + 1 protein") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 2.5,force = TRUE,fontface = "bold", 
                           segment.size = 0.2, box.padding = 0.5)  
  # ggtitle(paste0(ndropsfull, " total background drops"))
ggsave(pmain, filename = paste0(figpath, "/mean_stainedvs_unstained.pdf"), height = 4.2, width = 4.2)

# mean unstained vs mean emty drops for supplemental figure comparing definitions of empty drops 
# Library size distribution defined background 
p1 = ggplot(merged, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp + 
  ylab("unstained controls mean log10 + 1 protein") + 
  xlab(" Empty droplet mean log10 + 1 protein") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 2.5,force = TRUE,fontface = "bold", 
                           segment.size = 0.2, box.padding = 0.5) +  
  ggtitle(paste0(ndropsfull, " total background drops"))

######## Background defined with demultiplexing subset 
adt_neg_dmx = readRDS(file = here("data/V2_Data/background_data/adt_neg_dmx.rds")) %>%
  t() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("bc")

#  Empty / background droplets summarize data by protein stats. 
adt_neg_sub_summary =
  adt_neg_dmx %>% 
  rename(barcode_check = bc) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_negative_droplet = median(count), 
            mean_negative_droplet = mean(count), 
            var_negative_droplet = var(count)) 

# new merged data with hashing definition of background 
merged2 = full_join(adt_neg_sub_summary, uadt, by = "prot")

# mean in unstained control cells vs mean in empty droplets. 
ndrops = adt_neg_dmx$bc %>% unique %>% length()
p2 = ggplot(merged2, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp + 
  ylab("unstained controls mean log10 + 1 protein") + 
  xlab(" Empty droplet mean log10 + 1 protein") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 2.5,force = TRUE,fontface = "bold", 
                           segment.size = 0.2, box.padding = 0.5) +  
  ggtitle(paste0(ndrops, " background drops defined by hashing"))

# combine plot 
p5 = plot_grid(p2, p1, ncol = 1)
ggsave(p5, filename = paste0(figpath, "/mean_stainedvs_unstained_BKGRND_COMPARE.pdf"), height = 8.5, width = 4.5)


######## 
# mixture model metadata 
md = readRDS(file = "V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds")
md = md %>% mutate(log_prot_lib = log10(nUMI_Prot))
av_lib_mu1 = md %>% select(mean.1, log_prot_lib)
av_lib_mu1 = apply(av_lib_mu1, 1, mean)
md$av_lib_mu1 = av_lib_mu1

# supplemental figs on correlated component of mu1 and mu2 
bin_plot =  list(geom_bin2d(bins = 200, show.legend = FALSE), 
                 theme_bw(), 
                 ggpubr::stat_cor(method = "pearson"), 
                 scale_fill_viridis_c(), 
                 geom_smooth(method = "lm")
                 )

p6 = ggplot(md, aes(x = av_lib_mu1, y = mean.2)) + bin_plot + xlab("average of µ1 and library size") + ylab("µ2")
p7 = ggplot(md, aes(x = mean.1, y = mean.2)) + bin_plot  + xlab("µ1") + ylab("µ2")
pg = plot_grid(p7,p6)
ggsave(pg,filename = paste0(figpath, "mean1_mean2cor.pdf"), width = 6, height = 3)


###########################################

# unused 

# # mean variance trend global 
# mv = list(theme_bw(),geom_point(shape = 21, fill = "red3"), xlim(c(0,2.7)), ylim(c(0,0.65)))
# p1 = ggplot(merged, aes(x = mean_negative_droplet, y = var_negative_droplet)) + mv
# p2 = ggplot(merged, aes(x = mean_unstained_control, y = var_unstained_control))  + mv
# p3 = ggplot(merged, aes(x = mean_cells, y = var_cells)) + mv + 
#   ylab("var: all stained cells") + 
#   xlab("mean: all stained cells")
# p4 = ggplot(merged, aes(x = p_mean1, y = p_var1)) + mv +
#   ylab("var: cells negative for prot.") + 
#   xlab("mean: cells negative for prot.")
# p = plot_grid(p1,p2,p4,p3)
# ggsave(p, filename = paste0(figpath, "mean_variance_controls.pdf"), width = 5, height = 5)
# 
# 
# 
# #### uncertain noize floor plots
# # calculate for each cell the percentage of proteins in that cell with a value > 1
# # repeat for raw and dsb normalized to illustrate the global zero centering
# raw_prot = readRDS(file = here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds"))
# raw_prot = raw_prot@assay$CITE@raw.data
# index1 = colnames(raw_prot)[1]; index2 = colnames(raw_prot)[length(colnames(raw_prot))]
# prots = rownames(raw_prot)
# ddf = raw_prot %>% 
#   as.matrix %>% as.data.frame %>% 
#   rownames_to_column("prot") %>% 
#   gather(cell, count, index1:index2) %>% 
#   group_by(cell) %>% 
#   summarize(raw_pg0 = FSA::perc(count, val = 1, dir = "gt")) %>% 
#   mutate(norm = "raw")
# 
# muraw = ddf %$% raw_pg0 %>% mean
# p = ggplot(ddf, aes(x = raw_pg0, fill = norm)) + 
#   geom_density(show.legend = FALSE) + 
#   theme_bw()+ 
#   ggsci::scale_fill_d3() + 
#   xlim(c(0,100)) + 
#   geom_vline(xintercept = muraw, linetype= "dashed", size = 2.5) + 
#   theme(axis.text.x = element_text(size = 18)) + 
#   xlab(" ") +
#   ggtitle(paste0( "RAW CITE-seq counts \n", "each cell: % proteins > 1")) + 
#   annotate("text", fontface = "bold", size = 6.7, x = 30, y=0.03, label = paste0("distribution: \n", ncol(norm_prot), "\n cells")) +
#   theme(plot.title = element_text(size = 18))
# ggsave(filename = paste0(figpath,"pct_g1_RAWnorm_distribution.pdf"),width = 3.5, height = 3.5)  
# 
# 
# # normalized distribution
# norm_prot = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
# norm_prot = norm_prot@assay$CITE@data
# index1 = colnames(norm_prot)[1]; index2 = colnames(norm_prot)[length(colnames(norm_prot))]
# prots = rownames(norm_prot)
# ddf2 = norm_prot %>% 
#   as.matrix %>% as.data.frame %>% 
#   rownames_to_column("prot") %>% 
#   gather(cell, count, index1:index2) %>% 
#   group_by(cell) %>% 
#   summarize(dsb_pg0 = FSA::perc(count, val = 1, dir = "gt")) %>% 
#   mutate(norm = "dsb")
# 
# mudsb = ddf2$dsb_pg0 %>% mean
# p = ggplot(ddf2, aes(x = dsb_pg0, fill = norm)) + 
#   geom_density(show.legend = FALSE) + 
#   theme_bw()+ 
#   ggsci::scale_fill_jama(alpha = 0.6) + 
#   xlim(c(0,100)) + 
#   geom_vline(xintercept = mudsb, linetype= "dashed", size = 2.5) + 
#   theme(axis.text.x = element_text(size = 18)) + 
#   xlab(" ") +
#   ggtitle(paste0( "DSB Normalized  \n", "each cell % proteins > 1")) + 
#   annotate("text", fontface = "bold", size = 6.7, x = 75, y=0.03, label = paste0("distribution: \n", ncol(norm_prot), "\n cells")) +
#   theme(plot.title = element_text(size = 18))
# ggsave(filename = paste0(figpath,"pct_g1_dsbnorm_distribution.pdf"),width = 3.5, height = 3.5)  



