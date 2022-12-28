### Generate correlation plot ###
library(corrplot)
library(fmsb)

####
# RUN ANALYSIS
####
site_code = read.csv("input_data/Micro Coding Scheme.csv", header = TRUE)

sites = site_code[site_code$X.1 %in% c("N","O", "H", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "Eq1", "Eq2", "Eq3", "Eq4", "Eq5", "Eq6", "Eq7"), c("X.1", "X.2")]

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

# Add perirectal/groin-specific VRE colonization
final_phase2$VREsiteR_G = apply(final_phase2, 1, FUN = function(x){
  ifelse(sum(as.numeric(as.character(x["VREsiteR"])), as.numeric(as.character(x["VREsiteG"])), na.rm = T) > 0, 1, 0)})

# Subset to environmental sites sampled at least 400 visits (~50% of all visits) and at least 10% prevalence
keep_sites = c(gsub("sampled", "", names(which(rev(sort(colSums(apply(final_phase2[,grep(paste(paste0("sampled",as.character(sites$X.1)), collapse = "|"), colnames(final_phase2))], 1:2, as.numeric), na.rm = T))) >= 400))), "R_G")


vre_sites = colSums(apply(final_phase2[,paste0("VREsite", keep_sites)], 1:2, as.numeric), na.rm = T) >=40

keep_sites = names(which(vre_sites))

vre_kappa_cor_mat = matrix(0, nrow = length(keep_sites), ncol = length(keep_sites),dimnames = list(c(keep_sites), c(keep_sites)))

vre_kappa_lowerCI_mat = matrix(0, nrow = length(keep_sites), ncol = length(keep_sites),dimnames = list(c(keep_sites), c(keep_sites)))

vre_kappa_upperCI_mat = matrix(0, nrow = length(keep_sites), ncol = length(keep_sites),dimnames = list(c(keep_sites), c(keep_sites)))

for (k1 in rownames(vre_kappa_cor_mat)){
  for (k2 in colnames(vre_kappa_cor_mat)){
    if (k1 == k2){next}
    
    kappa_cor = Kappa.test(final_phase2[,k1], final_phase2[,k2])
    vre_kappa_cor_mat[k1,k2] = kappa_cor$Result$estimate
    vre_kappa_lowerCI_mat[k1,k2] = kappa_cor$Result$conf.int[1]
    vre_kappa_upperCI_mat[k1,k2] = kappa_cor$Result$conf.int[2]
  }
}

rownames(vre_kappa_cor_mat) = colnames(vre_kappa_cor_mat) = rownames(vre_kappa_lowerCI_mat) = colnames(vre_kappa_lowerCI_mat) = rownames(vre_kappa_upperCI_mat) = colnames(vre_kappa_upperCI_mat) = sapply(keep_sites,FUN = function(x){
  ifelse(x %in% "VREsiteR_G", "Groin/Perirectal", as.character(sites[sites$X.1 %in% gsub("VREsite", "", x), "X.2"]))
})

# Make correlation plot
file = "output/site_contamination_corrplot.pdf"
pdf(file)

cr = corrplot(vre_kappa_cor_mat, order = "hclust", 
              type="upper",  method = "color", tl.cex = 1, tl.col="black", number.cex = .75, tl.srt = 45, mar = c(0,0,1,0), outline = T, hclust.method = "median")

conf_pos_mat = matrix(FALSE, nrow = nrow(cr$corr), ncol = ncol(cr$corr), dimnames = list(rownames(cr$corr), colnames(cr$corr)))
conf_pos_mat[upper.tri(conf_pos_mat)] = T

conf_pos = apply(which(upper.tri(conf_pos_mat), arr.ind = T), 1, FUN = function(x){
  row_x = rownames(conf_pos_mat)[x[1]]
  col_x = colnames(conf_pos_mat)[x[2]]
  paste0(row_x, "_", col_x)
})

for (p in conf_pos){
  xi = strsplit(p, "_")[[1]][1]
  yi = strsplit(p, "_")[[1]][2]
  
  corr_clust = vre_kappa_cor_mat[xi, yi]
  upperCI = vre_kappa_upperCI_mat[xi, yi]
  lowerCI = vre_kappa_lowerCI_mat[xi, yi]
  CI = paste0(round(lowerCI, 2), "-", round(upperCI, 2))
  
  xi_pos = which(colnames(conf_pos_mat) %in% yi)
  yi_pos = which(rev(rownames(conf_pos_mat)) %in% xi)
  
  text(xi_pos, yi_pos + .4, round(corr_clust,2), cex = 1, pos = 1)   
  text(xi_pos, yi_pos + 0, CI, cex = 0.7, pos = 1)
}

dev.off()