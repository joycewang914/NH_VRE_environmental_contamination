# Investigate the effect of low and high-risk on gut microbiome using NMDS

library(lubripack)
libraries = c("vegan", "ggplot2", "grid", "ggpubr", "tidyverse", "ggvegan", "gridExtra", "ggrepel")
lubripack(libraries)
source("lib/mothur_taxonomy_munging.R")


final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")
abx_risk_ind = structure(c('none', 'low', 'high'), names = 0:2)
final_phase2$abx_30d_prior_risk = sapply(final_phase2$abx_30d_prior, FUN = function(x){
  abx_risk_ind[as.character(x)]
})

tax_filename = "input_data/eswab.final.0.03.cons.taxonomy"
taxonomy = mothur_tax_munging_func(tax_filename)

# Get count data
otu_ind = colnames(final_phase2)[grepl("Otu", colnames(final_phase2))]
raw_otu_ind = otu_ind[grep("RA", otu_ind, invert = T)]
otu_df = apply(final_phase2[,raw_otu_ind], 1:2, as.numeric)
otu_df = otu_df[rowSums(otu_df, na.rm = T) > 0, ] # remove entries without microbiome data

####
# RUN ANALYSIS
####

# add 1 count to all otu groups for CLR transformation
data_offset = otu_df + 1

result.filter = low.count.removal(data_offset, percent = 0.05)
data.filter <- result.filter$data.filter
top_otu_ind = names(sort(colSums(data.filter > 1), decreasing = T)[1:10]) # top 10 prevalence
# log transform
log_top_otu_df = log1p(data.filter)

# set.seed(123)
# otu_BC.nmds = metaMDS(log_top_otu_df, distance='bray', k=2, trymax=1000, autotransform = F)
# Run 329 stress 0.1911776 
# ... Procrustes: rmse 0.0002504323  max resid 0.003903709 
# ... Similar to previous best
# *** Solution reached

# saveRDS(object = otu_BC.nmds, file = "microbiome_otu_df_nmds.RDS")
# 
otu_BC.nmds = readRDS("input_data/microbiome_otu_df_nmds.RDS")
stressplot(otu_BC.nmds)

# Anosim to test microbiome clustering by antibiotics  grouping
all_compare_ano = anosim(log_top_otu_df, grouping = final_phase2[rownames(log_top_otu_df), "abx_30d_prior_risk"], distance = "bray", permutations = 9999)

# anosim to test microbiome clustering by colonization pattern
all_compare_ano_by_col = anosim(log_top_otu_df, grouping = final_phase2[rownames(log_top_otu_df), "vre_col"], distance = "bray", permutations = 9999)


# pairwise abx comparisons
abx_pairs = t(combn(abx_risk_ind, 2))

pairwise_anosim = apply(abx_pairs, 1,  FUN = function(x){
  abx_ind = intersect(rownames(log_top_otu_df), rownames(final_phase2)[final_phase2$abx_30d_prior_risk %in% x])
  pair_anosim = anosim(as.matrix(log_top_otu_df[abx_ind, ]), final_phase2[abx_ind, "abx_30d_prior_risk"], distance = 'bray', permutations = 9999)
  pair_anosim$signif
})

# pairwise vre_col comparisons
vre_col_pairs = t(combn(0:3, 2))

vre_col_pairwise_anosim = apply(vre_col_pairs, 1,  FUN = function(x){
  vre_col_ind = intersect(rownames(log_top_otu_df), rownames(final_phase2)[final_phase2$vre_col %in% x])
  pair_anosim = anosim(as.matrix(log_top_otu_df[vre_col_ind, ]), final_phase2[vre_col_ind, "vre_col"], distance = 'bray', permutations = 9999)
  pair_anosim$signif
})

names(vre_col_pairwise_anosim) = apply(vre_col_pairs, 1, FUN = function(x){
  paste0(x[1], "_", x[2])
})
# Plot microbiome on NMDS
# first panel showing distribution of microbiome samples
abx_otu_BC.nmds = data.frame(NMDS1 = scores(otu_BC.nmds)[,1], 
                             NMDS2 = scores(otu_BC.nmds)[,2],
                             abx = as.factor(final_phase2[rownames(scores(otu_BC.nmds)), "abx_30d_prior_risk"]))

abx_otu_BC.nmds$Exposure = sapply(abx_otu_BC.nmds$abx, FUN = function(x){
  if (x == 'high'){full_x = "High-risk antibiotics"}
  if (x == 'low'){full_x = "Low-risk antibiotics"}
  if (x == 'none'){full_x = "No antibiotics"}
  full_x
})

full_abx_risk_ind = sapply(abx_risk_ind, FUN = function(x){
  if (x == 'high'){full_x = "High-risk antibiotics"}
  if (x == 'low'){full_x = "Low-risk antibiotics"}
  if (x == 'none'){full_x = "No antibiotics"}
  full_x
})

p1 = ggplot(abx_otu_BC.nmds, aes(x = NMDS1, y=NMDS2, color=Exposure, fill=Exposure)) +
  stat_ellipse(level=0.75, geom='polygon', alpha=0.2, show.legend=F) +
  geom_point() +
  coord_fixed()+
  labs(x = 'NMDS1', y = 'NMDS2') +
  scale_color_manual(name=NULL,
                     breaks=full_abx_risk_ind,
                     values = c("grey", "blue", "red"),
                     labels=full_abx_risk_ind) +
  scale_fill_manual(name = NULL,
                    breaks = full_abx_risk_ind,
                    values = c("lightgrey", "dodgerblue", "pink"),
                    labels = full_abx_risk_ind
  ) +
  theme_classic() +
  theme(legend.text=element_text(size=20), 
        legend.title=element_text(size=20), 
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.position = "top",
        legend.direction = "vertical")+
  guides(colour = guide_legend(override.aes = list(size=10)))
  # theme(
  #   legend.key.size=unit(0.5, 'cm'),
  #   legend.position=c(0.95, 0.95),
  #   legend.background = element_rect(fill='NA',color='black'),
  #   legend.margin = margin(t=-2, r=3, b=3, l=3)) +
  # labs(subtitle = paste0('stress = ', round(otu_BC.nmds$stress, 3),
  #                        "\nANOSIM R = ", round(all_compare_ano$statistic, 3),
  #                        "\nANOSIM significance ", ifelse(all_compare_ano$signif < 0.001, "< 0.001", paste0("= ", round(all_compare_ano$signif, 3)))))

ggsave(path = "output/", filename = "abx_nmds_plot.pdf", width = 20, height = 10, device='pdf', dpi=1400)

# second panel showing species driving community separation
fort = fortify(otu_BC.nmds)
fort$Label_genus = sapply(fort$Label, FUN = function(x){
  x = as.character(x)
  new_label = ifelse(grepl("Otu", x), taxonomy[taxonomy$OTU %in% x, "genus"], x)
  new_label
})
p2 = ggplot() +
  coord_fixed() +
  geom_abline(intercept=0, slope=0, linetype = 'dashed', size=.8, colour='grey')+
  geom_vline(aes(xintercept=0), linetype='dashed', size=0.8, colour = 'grey')+
  geom_point(data=subset(fort, Score=='sites'),
             mapping = aes(x = NMDS1, y = NMDS2),
             colour = 'lightgrey',
             alpha = 0.5) +
  geom_segment(data = subset(fort[fort$Label %in% top_otu_ind,], Score=='species'),
               mapping = aes(x=0, y =0, xend = NMDS1, yend = NMDS2),
               arrow=arrow(length=unit(0.015, 'npc'),
                           type = 'closed'),
               colour='#35B5AC',
               size=0.8) +
  geom_text_repel(aes(x = fort[fort$Label %in% top_otu_ind,"NMDS1"], 
                      y = fort[fort$Label %in% top_otu_ind,"NMDS2"], 
                      label = fort[fort$Label %in% top_otu_ind,"Label_genus"]),
                  size = 6, 
                  force_pull = 0.2, 
                  min.segment.length = 0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour='black'),
        legend.margin = margin(t=-2, r=3, b=3, l=3),
        legend.text=element_text(size=45), 
        legend.title=element_text(size=45), 
        axis.title.y = element_text(size=45),
        axis.title.x = element_text(size=45), 
        axis.text = element_text(size=45))

# grid.arrange(p1, p2, ncol = 1)
# g = arrangeGrob(p1, p2, ncol = 1)

ggsave(path = "output/", filename = "abx_nmds_vectors.pdf", width = 10, height = 10, device='pdf', dpi=1400)
