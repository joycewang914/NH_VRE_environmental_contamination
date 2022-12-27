# Use consort diagram to clarify data usage
# Use R package "consort" to generate plot
# https://cran.r-project.org/web/packages/consort/vignettes/consrot_diagram.html
library(consort)
library(grid)
library(tidyverse)

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")
site_code = read.csv("input_data/Micro Coding Scheme.csv", header = TRUE)
sites = site_code[site_code$X.1 %in% c("O", "H", "G", "R",  "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "Eq1", "Eq2", "Eq3", "Eq4", "Eq5", "Eq6", "Eq7"), c("X.1", "X.2")]

# Organize relevant data to be entered in consort diagram
# patient ID
num_pts = unique(final_phase2$studyID)

# patient visit
num_visits = nrow(final_phase2)

# patient visit breakdown
num_visit_breakdown = table(sapply(unique(final_phase2$studyID), FUN = function(x){sum(final_phase2$studyID %in% x)}))
num_visit_breakdown_list = paste0(paste0(num_visits, " visits total\n"), paste(sapply(1:length(num_visit_breakdown), FUN = function(x){paste0(num_visit_breakdown[x], " patients: ", names(num_visit_breakdown)[x], " visits")}), collapse = "\n"))

# environmental samples
num_env_samples = sum(final_phase2[,grep("sampledE", colnames(final_phase2))] > 0, na.rm = T)

# patient samples
num_pt_samples = sum(final_phase2[,grep(paste("sampled", c("G", "H", "O", "R"), collapse = "|", sep = ""), colnames(final_phase2))] > 0, na.rm = T)

# sample breakdown
sample_breakdown = colSums(final_phase2[,grepl(paste(sapply(sites$X.1, FUN = function(x){paste0("sampled", x)}), collapse = "|"),colnames(final_phase2))], na.rm = T)
names(sample_breakdown) = sapply(names(sample_breakdown), FUN = function(x){
  as.character(sites[sites$X.1 %in% gsub("sampled", "", x),"X.2"])
})

sample_breakdown_gt0 = rev(sort(sample_breakdown[sample_breakdown  > 0]))
sample_breakdown_gt100 = c(sample_breakdown_gt0[sample_breakdown_gt0 > 100], Other = sum(sample_breakdown_gt0[sample_breakdown_gt0 <100]))

sample_breakdown_gt100_list = paste(sapply(names(sample_breakdown_gt100), FUN = function(x){paste0(x,": ", sample_breakdown_gt100[x], " samples")}), collapse = "\n")

# microbiome samples
num_microbiome_samples = sum(!is.na(final_phase2$invsimp))
num_microbiome_sample_breakdown = table(sapply(num_pts, FUN = function(x){sum(!is.na(final_phase2[final_phase2$studyID %in% x, "invsimp"]))}))
num_microbiome_breakdown_list = paste0(paste(sapply(1:length(num_microbiome_sample_breakdown), FUN = function(x){paste0(num_microbiome_sample_breakdown[x], " patients: ", names(num_microbiome_sample_breakdown)[x], " samples")}), collapse = "\n"))

# Plot
options(txt_gp = gpar(cex = 0.8))

txt1 <- paste0("Enrolled patients (n=", length(num_pts), ")")
txt1_side <- num_visit_breakdown_list
txt2 = paste0("Patient samples (n=", num_pt_samples, ")\nEnvironmental samples (n=", num_env_samples, ")")
txt2_side = sample_breakdown_gt100_list
txt3 = paste0("Perirectal E-swabs\n(n=", num_microbiome_samples, ")")
txt3_side = num_microbiome_breakdown_list

g <- add_box(txt = txt1) %>%
  add_side_box(txt = txt1_side) %>%
  add_box(txt = "Sample collection") %>%
  add_split(txt = c("VRE culture data", "Gut microbial community data")) %>%
  add_box(c(txt2, txt3)) %>%
  add_box(c(txt2_side, txt3_side)) %>%
  add_box(c("Longitudinal analysis", "Cross-sectional analysis"))


png("output/consort_diagram.png", width = 29, 
    height = 30, res = 300, units = "cm", type = "cairo") 
plot(g)
dev.off() 
