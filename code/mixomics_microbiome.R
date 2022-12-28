# Use MixOmics to identify differences between microbiomes exposed and non-exposed to any antibiotics
library(mixOmics)
source("lib/mothur_taxonomy_munging.R")

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

tax_filename = "input_data/eswab.final.0.03.cons.taxonomy"
taxonomy = mothur_tax_munging_func(tax_filename)

# need count data
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
# check the number of variables kept after filtering
# in this particular case we had already filtered the data so no was change made, but from now on we will work with 'data.filter'

abx_groups = as.factor(sapply(rownames(otu_df), FUN = function(x){
  ifelse(final_phase2[x,"abx_30d_prior"] > 0, "Yes", "No")
  final_phase2[x,"abx_30d_prior"]
}))


splsda.result <- splsda(data.filter, abx_groups, keepX = c(50,30)) # run the method


# extract the variables used to construct the first latent component
otu_ind = selectVar(splsda.result, comp = 1)$name 
otu_genus = sapply(names(result.filter$keep.otu), FUN = function(x){
  taxonomy[taxonomy$OTU %in% x, "genus"]
})

# depict weight assigned to each of these variables
file ="output/mixomics_plotloadings.png"
png(file, 640, 480)
par(mar = c(5, 18, 2, 1))

loading_cols = structure(c("#FFCB05", "#00274C"), names = c("Yes", "No"))

top_driver_ind = names(sort(abs(splsda.result$loadings$X[,1]), decreasing = T)[1:20])
top_driver_val = splsda.result$loadings$X[top_driver_ind,1]
names(top_driver_val) = otu_genus[names(top_driver_val)]
ordered_top_drivers = sort(top_driver_val, decreasing = T)
ordered_top_drivers_cols = sapply(ordered_top_drivers, FUN = function(x){
  ifelse(x < 0, loading_cols["No"], loading_cols["Yes"])
})

bp = barplot(ordered_top_drivers, horiz = T, las = 2, xlab = "Component 1", border = NA, col = ordered_top_drivers_cols, xlim = c(round(min(ordered_top_drivers), 1), round(max(ordered_top_drivers), 1)))

segments(x0 = round(min(ordered_top_drivers), 1), y0 = min(bp), x1 = round(min(ordered_top_drivers), 1), y1 = max(bp), col = "black", lwd = 2)

for (b in bp){
  print(b)
  segments(x0 = round(min(ordered_top_drivers), 1) - 0.01,
           y0 = b,
           x1 = round(min(ordered_top_drivers), 1),
           y1 = b,
           col = "black",
           xpd = TRUE)
}

legend("topright",
       legend=names(loading_cols), 
       xpd = TRUE, 
       bty = "n", 
       border = NA, 
       col = loading_cols, 
       fill = loading_cols,
       title = "                       Prior antibiotic exposure", 
       cex = 1.25, 
       adj = 0, 
       xjust = 0,  
       y.intersp = 1.25)


dev.off()
