### Generate table 1 including descriptive statistics of patients at enrollment###

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

all_factors = c("age_bl", "sex_bl", "charlson", "psmscore", "hospstay", "urincath", "abx_30d_prior")

binary = c("sex_bl", "urincath") #sex: M1F0

phase2_bl = final_phase2[final_phase2$visit %in% 0, ]

p_val_denote = function(int){
  if (int >= 0.1){p_symbol = ""}
  if (int < 0.1){p_symbol = "*"}
  if (int < 0.05){p_symbol = "**"}
  if (int < 0.01){p_symbol = "***"}
  p_symbol
}
  
var_list = list()
for (x in all_factors){
  
  for (v in 0:3){
    
    # all statistical comparisons with vre_col = 0
    
    if (!x %in% binary & !x %in% "abx_30d_prior"){
      v_mean = round(mean(phase2_bl[phase2_bl$vre_col %in% v, x], na.rm = T), 1)
      v_sd = round(sd(phase2_bl[phase2_bl$vre_col %in% v, x], na.rm = T), 1)
      v_pval = t.test(phase2_bl[phase2_bl$vre_col %in% 0, x], phase2_bl[phase2_bl$vre_col %in% v, x])$p.value
      v_pval_symbol = p_val_denote(v_pval)
      var_list[[x]][[as.character(v)]] = paste0(v_mean, " (", v_sd, ")", v_pval_symbol)
      }
  
    if (x %in% "abx_30d_prior"){
      
      for (a in 1:2){
        abx_ind = ifelse(a == 1, "low_risk_abx", "high_risk_abx")
        
        a_num = sum(phase2_bl$vre_col %in% v & phase2_bl[,x] %in% a)
        a_prop = a_num / sum(phase2_bl$vre_col %in% v)
        a_table = table(phase2_bl[,x], phase2_bl$vre_col)
        a_table_subset = a_table[rownames(a_table) %in% c(0, a), colnames(a_table) %in% c(0, v)]
        a_pval = ifelse(v == 0, 1, fisher.test(a_table_subset)$p.value)
        a_pval_symbol = p_val_denote(a_pval)
        
        var_list[[abx_ind]][[as.character(v)]] = paste0(a_num, " (", round(a_prop*100, 1), ")", a_pval_symbol)
      }
    }
    
    if (x %in% binary){
      v_num = sum(phase2_bl$vre_col %in% v & phase2_bl[,x] %in% 1)
      v_prop = round(v_num / sum(phase2_bl$vre_col %in% v)*100, 1)
      v_table = table(phase2_bl[,x], phase2_bl$vre_col)
      v_pval = ifelse(v == 0, 1, fisher.test(v_table[, colnames(v_table) %in% c(0, v)])$p.value)
      v_pval_symbol = p_val_denote(v_pval)
      var_list[[x]][[as.character(v)]] = paste0(v_num, " (", v_prop, ")", v_pval_symbol)
      }
    }
  }
  

clinical_var_mat = do.call(rbind, var_list)

clinical_vars = rownames(clinical_var_mat)
rownames(clinical_var_mat) = sapply(clinical_vars, FUN = function(x){
  temp = x
  if (x %in% "age_bl"){temp = "Age (Mean, SD)"}
  if (x %in% "sex_bl"){temp = "Male sex (%)"}
  if (x %in% "urincath"){temp = "Urinary catheter use within 30 days (%)"}
  if (x %in% "hospstay"){temp = "Length of prior hospital stay (Mean, SD)"}
  if (x %in% "charlson"){temp = "Charlson score (Mean, SD)"}
  if (x %in% "psmscore"){temp = "Physical self-maintenance score (Mean, SD)"}
  if (x %in% "low_risk_abx"){temp = "Exposure to low-risk antibiotic (%)"}
  if (x %in% "high_risk_abx"){temp = "Exposure to high-risk antibiotic (%)"}
  
  temp
  
})

colnames(clinical_var_mat) = sapply(0:3, FUN = function(x){
  if (x == 0){col_type = "Uncolonized"}
  if (x == 1){col_type = "Environment-only"}
  if (x == 2){col_type = "Patient-only"}
  if (x == 3){col_type = "Both"}
  paste0(col_type, "\n(N = ", table(phase2_bl$vre_col)[as.character(x)], ")" )
})

write.csv(file = "output/table1.csv", x = clinical_var_mat)
  