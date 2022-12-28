# General stats on patient samples

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")
phase2_bl = final_phase2[final_phase2$visit %in% 0, ]

# number of patient samples
sum(final_phase2[,grep(paste("sampled", c("G", "H",  "O", "R"), collapse = "|", sep = ""), colnames(final_phase2))] > 0, na.rm = T)

# number of environmental samples
sum(final_phase2[,grep("sampledE", colnames(final_phase2))] > 0, na.rm = T)

# reason for antibiotic use
path_abx = read.csv("/Users/joycewang/Desktop/Mody_NH_collections/Pathways/Antibiotic Use.csv")


abx_reasons = sapply(unique(final_phase2$studyID), FUN = function(x){
  reason = path_abx[path_abx$studyID %in% x & path_abx$visit %in% 0, "abxindic"]
  if (length(reason) == 0){
    reason = "no abx"
  }
  paste(as.character(unique(tolower(reason))), collapse = "; ")
})

# make a table comparing reasons for antibiotic use in colonization groups 2 & 3
# raw_abx_indic = abx_reasons[phase2_bl[phase2_bl$vre_col %in% 2:3, "studyID"]]
raw_abx_indic = abx_reasons

# manually set patients 2505 and 2514 to no abx since no information is provided
raw_abx_indic[setdiff(phase2_bl[phase2_bl$antibiotic %in% 1, "studyID"], phase2_bl[phase2_bl$abx_30d_prior > 0, "studyID"] )] = 'no abx'

curated_abx_indic = sapply(raw_abx_indic, FUN = function(x){

  curated_reasons = c()
  if (grepl('uti|uri|cystitis|pcp', tolower(x))){curated_reasons = c(curated_reasons, 'UTI-related')}
  if (grepl('pneu|pleu|copd|resp', tolower(x))){curated_reasons = c(curated_reasons,'Respiratory-related')}
  if (grepl('osteo|bone', tolower(x))){curated_reasons = c(curated_reasons,'Bone-related')}
  if (grepl('s*g', tolower(x))){curated_reasons = c(curated_reasons, 'Surgery-related')}
  if (grepl('sepsis|gastritis|c-diff|cellulitis|abscess', tolower(x))){curated_reasons = c(curated_reasons, 'Sepsis/enteric/soft-tissue infections')}
  if (grepl('no abx', tolower(x))){curated_reasons = c('No antibiotics')}
  if (length(curated_reasons) == 0){
    curated_reasons = "Unknown"
  }
  if (length(curated_reasons) > 1){
    curated_reasons = "More than one"
  }
  curated_reasons
})

# function summarizing p value
p_val_denote = function(int){
  if (int >= 0.1){p_symbol = ""}
  if (int < 0.1){p_symbol = "*"}
  if (int < 0.05){p_symbol = "**"}
  if (int < 0.01){p_symbol = "***"}
  p_symbol
}

# make summary table with p values
curated_reason_mat = matrix(0, nrow = 8, ncol = 4,
                           dimnames = list(
                             c("UTI-related", "Respiratory-related", "Bone-related", "Surgery-related", "Sepsis/enteric/soft-tissue infections", "More than one", "Unknown", "No antibiotics"),
                             c(paste0("Uncolonized\n(N=", length(phase2_bl[phase2_bl$vre_col %in% 0, "studyID"]), ")"),
                               paste0("Environment-only\n(N=", length(phase2_bl[phase2_bl$vre_col %in% 1, "studyID"]), ")"),
                               paste0("Patient-only\n(N=", length(phase2_bl[phase2_bl$vre_col %in% 2, "studyID"]), ")"), 
                               paste0("Both\n(N=", length(phase2_bl[phase2_bl$vre_col %in% 3, "studyID"]), ")"))))


for (r in rownames(curated_reason_mat)){
  
  for (o in 0:3){
    
    
    group_types = table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% o, "studyID"]])
    
    if (!r %in% names(group_types)){group_num = 0}else{group_num = group_types[r]}
    
    group_pct = round(group_num / length(phase2_bl[phase2_bl$vre_col %in% o, "studyID"]) * 100, 1)
    
  if (o == 0){
    curated_reason_mat[r, o+1] =paste0(group_num, " (", group_pct, ")")
  }else{
    stat_table = matrix(c(table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% 0, "studyID"]])[r],
                          sum(table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% 0, "studyID"]])[!names(table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% o, "studyID"]])) %in% r]),
                          table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% o, "studyID"]])[r],
                          sum(table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% o, "studyID"]])[!names(table(curated_abx_indic[phase2_bl[phase2_bl$vre_col %in% o, "studyID"]])) %in% r])), nrow=2, ncol=2, byrow = T)
    
    stat_table[is.na(stat_table)] = 0
    stat_result = fisher.test(stat_table)
    p_val = p_val_denote(stat_result$p.value)
    
    curated_reason_mat[r,o+1] = paste0(group_num, " (", group_pct, ")",p_val)
  }
  }
}
rownames(curated_reason_mat) = paste0(rownames(curated_reason_mat), " (%)")

write.csv(file = "output/antibiotic_use_reason.csv", x = curated_reason_mat)
