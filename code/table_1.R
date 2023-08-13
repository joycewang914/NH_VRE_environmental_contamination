### Generate table 1 including descriptive statistics of patients at enrollment###
setwd("/Users/joycewang/Desktop/Snitkin_lab/Manuscripts/Pathways_mbiome_environment/NH_env_microbiome/analysis/")

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

all_factors = c("age_bl", "sex_bl", "charlson", "psmscore", "hospstay", "urincath", "abx_30d_prior", "race_bl", "intervention", "novisits")
#sex: M1F0; race: white1black2; intervention: control0intervention1

binary = c("urincath",  "intervention") 

phase2_bl = final_phase2[final_phase2$visit %in% 0, ]
phase2_bl$race_bl[phase2_bl$race_bl %in% "99"] = 1

p_val_denote = function(int){
  if (int >= 0.1){p_symbol = ""}
  if (int < 0.1){p_symbol = "*"}
  if (int < 0.05){p_symbol = "**"}
  if (int < 0.01){p_symbol = "***"}
  p_symbol
}
  
var_list = list()
for (x in all_factors){
  print(x)
  
  for (v in 0:3){
    
    # all statistical comparisons with vre_col = 0
    
    if (!x %in% binary & !x %in% c("abx_30d_prior", "race_bl", "sex_bl")){
      v_mean = round(mean(phase2_bl[phase2_bl$vre_col %in% v, x], na.rm = T), 1)
      v_sd = round(sd(phase2_bl[phase2_bl$vre_col %in% v, x], na.rm = T), 1)
      v_pval = t.test(phase2_bl[phase2_bl$vre_col %in% 0, x], phase2_bl[phase2_bl$vre_col %in% v, x])$p.value
      v_pval_symbol = p_val_denote(v_pval)
      var_list[[x]][[as.character(v)]] = paste0(v_mean, " (", v_sd, ")", v_pval_symbol)
      }
  
    if (x %in% c("abx_30d_prior", "race_bl", "sex_bl")){
      
      for (a in names(table(phase2_bl[[x]]))){
        
        if (x == "abx_30d_prior"){
          if (a == "0"){next}
          a_ind = ifelse(a == 1, "low_risk_abx", "high_risk_abx")
          a_num = sum(phase2_bl$vre_col %in% v & phase2_bl[,x] %in% a)
          a_prop = a_num / sum(phase2_bl$vre_col %in% v)
          a_table = table(phase2_bl[,x], phase2_bl$vre_col)
          a_table_subset = a_table[rownames(a_table) %in% c(0, a), colnames(a_table) %in% c(0, v)]
          a_pval = ifelse(v == 0, 1, fisher.test(a_table_subset)$p.value)
        }
        # only perform one fisher test for race and sex since there are only two groups (African American and Women are the control group, respectively)
        if (x == "race_bl"){
          a_ind = ifelse(a == 1, "Non-hispanic White", "African American")
          a_num = sum(phase2_bl$vre_col %in% v & phase2_bl[,x] %in% a)
          a_prop = a_num / sum(phase2_bl$vre_col %in% v)
          a_table = table(phase2_bl[,x], phase2_bl$vre_col)
          a_pval = ifelse(v == 0, 1, fisher.test(a_table[, colnames(v_table) %in% c(0, v)])$p.value)
        }
        
        if (x == "sex_bl"){
          a_ind = ifelse(a == 1, "Men", "Women")
          a_num = sum(phase2_bl$vre_col %in% v & phase2_bl[,x] %in% a)
          a_prop = a_num / sum(phase2_bl$vre_col %in% v)
          a_table = table(phase2_bl[,x], phase2_bl$vre_col)
          a_pval = ifelse(v == 0, 1, fisher.test(a_table[, colnames(v_table) %in% c(0, v)])$p.value)}

        a_pval_symbol = p_val_denote(a_pval)
        
        var_list[[a_ind]][[as.character(v)]] = paste0(a_num, " (", round(a_prop*100, 1), ")", a_pval_symbol)
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
clinical_var_mat = rbind(clinical_var_mat, Race = rep(NA ,4), Sex = rep(NA, 4))

clinical_var_mat = clinical_var_mat[c("age_bl", "Sex", "Women", "Men", "Race", "Non-hispanic White", "African American", "charlson", "psmscore", "hospstay", "urincath","low_risk_abx", "high_risk_abx", "intervention", "novisits"),]

clinical_vars = rownames(clinical_var_mat)
rownames(clinical_var_mat) = sapply(clinical_vars, FUN = function(x){
  temp = x
  if (x %in% "age_bl"){temp = "Age (Mean, SD)"}
  if (x %in% "Sex"){temp = "Sex (%)"}
  if (x %in% "Race"){temp = "Race (%)"}
  if (x %in% "urincath"){temp = "Urinary catheter use within 30 days (%)"}
  if (x %in% "hospstay"){temp = "Length of prior hospital stay (Mean, SD)"}
  if (x %in% "charlson"){temp = "Charlson score (Mean, SD)"}
  if (x %in% "psmscore"){temp = "Physical self-maintenance score (Mean, SD)"}
  if (x %in% "low_risk_abx"){temp = "Exposure to low-risk antibiotic (%)"}
  if (x %in% "high_risk_abx"){temp = "Exposure to high-risk antibiotic (%)"}
  if (x %in% "Women"){temp = " Women"}
  if (x %in% "Men"){temp = " Men"}
  if (x %in% "African American"){temp = " African American"}
  if (x %in% "Non-hispanic White"){temp = " Non-hispanic white"}
  if (x %in% "intervention"){temp = "Intervention (%)"}
  if (x %in% "novisits"){temp = "Number of follow-up visits (Mean, SD)"}
  
  temp
  
})

colnames(clinical_var_mat) = sapply(0:3, FUN = function(x){
  if (x == 0){col_type = "Uncolonized"}
  if (x == 1){col_type = "Environment-only"}
  if (x == 2){col_type = "Patient-only"}
  if (x == 3){col_type = "Both"}
  paste0(col_type, "\n(N = ", table(phase2_bl$vre_col)[as.character(x)], ")" )
})


write.csv(file = "output/table1.csv", x = clinical_var_mat, na = "")
  