setwd("/Users/joycewang/Desktop/Snitkin_lab/Manuscripts/Pathways_mbiome_environment/NH_env_microbiome/Analysis/")

# Construct GEE models to identify drivers for environmental contamination

source("lib/QIC_binom_geeglm.R")
library(geepack)
library(pstools)

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")
final_phase2$wave = as.numeric(factor(final_phase2$visit))

# difference between patient-only and both colonization
outcomes = c(2, 3)
mods = list()

# subset to rows with patient-only or both colonization, antibiotic, PSMS, sex, and intervention information
all_factors = c("age_bl", "sex_bl", "charlson", "psmscore", "hospstay", "urincath", "abx_30d_prior", "vre_col")

not_na_ind = sapply(all_factors, FUN = function(x){
  which(!is.na(final_phase2[,x]))
})

keep_ind = Reduce(intersect, not_na_ind)

pt_both_df = final_phase2[intersect(which(final_phase2$vre_col %in% outcomes), keep_ind),]

# model including PSMS score

# choosing data correlation structure 
# follow this lecture: https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm

corr_structure = c("independence", "exchangeable", "ar1")
corstr_mods = list()
for (corstr in corr_structure){
  print(corstr)
  geeglmfit = geeglm(vre_col == 3 ~ age_bl + charlson + hospstay + factor(urincath) + factor(abx_30d_prior > 0) + psmscore + factor(sex_bl) + factor(intervention), id=studyID, family=binomial, corstr=corstr, data=pt_both_df, scale.fix = T, waves=wave)
  corstr_mods[[corstr]] = geeglmfit
}

model_compare = lapply(corstr_mods, FUN = function(m){
  QIC.binom.geeglm(m, corstr_mods[["independence"]])
})

# Both QIC and CIC rate the AR1 model as best
all_mod_fit = gee_stepper(corstr_mods$ar1, formula(corstr_mods$ar1))

# selected variables
names(coef(all_mod_fit))[!names(coef(all_mod_fit)) %in% "(Intercept)"]

all_step_mod = geeglm(vre_col == 3 ~ factor(abx_30d_prior > 0) + psmscore + factor(sex_bl), id=studyID, family=binomial, corstr='ar1', data=pt_both_df, scale.fix = T, waves=wave)

mod_df = tidy(all_step_mod, conf.int = T, conf.level = 0.95)

mod_mat = matrix(NA, nrow = length(mod_df$term[!mod_df$term %in% "(Intercept)"]), ncol = 3,
                 dimnames = list(mod_df$term[!mod_df$term %in% "(Intercept)"], 
                                 c("aOR", "95% CI", "P value")))
for (m in rownames(mod_mat)){
  mod_mat[m,"aOR"] = as.character(round(exp(mod_df[mod_df$term %in% m,"estimate"]), 2))
  mod_mat[m, "95% CI"] = paste0(round(exp(mod_df[mod_df$term %in% m,"conf.low"]), 2), " - ", round(exp(mod_df[mod_df$term %in% m,"conf.high"]), 2))
  mod_mat[m, "P value"] = ifelse(mod_df[mod_df$term %in% m,"p.value"] < 0.01, "< 0.01", as.character(round(mod_df[mod_df$term %in% m,"p.value"], 2)))
}

rownames(mod_mat) = sapply(rownames(mod_mat), FUN = function(x){
  if (x == "factor(abx_30d_prior > 0)TRUE"){col_x = "Recent antibiotic exposure"}
  if (x == "psmscore"){col_x = "Physical self-maintenance scale"}
  if (x == "factor(sex_bl)1"){col_x = "Male sex"}
  col_x
})

write.csv(file = "output/multivariate_overall_psms.csv", x = mod_mat)



# determine which PSMS activity drives the inverse relationship between PSMS score and environmental contamination

psms_cats = c("psmstoil", "psmsfeed", "psmsdress", "psmsgroom", "psmsamb", "psmsbath")

psms_mods = list()

psms_mod_mat = matrix(NA, nrow = 3, ncol = length(psms_cats), dimnames = list(c('Prior antibiotic exposure', "PSMS activity score", "Male sex"), psms_cats))

for (p in psms_cats){
  print(p)
  p_formula = as.formula(paste0("vre_col == 3 ~ factor(abx_30d_prior > 0) +", p, "+ factor(sex_bl)"))
  mod.ar1 = geeglm(p_formula, id=studyID, family=binomial, corstr='ar1', data=pt_both_df, scale.fix = T, waves=wave)
  psms_mod_df = tidy(mod.ar1, conf.int = T)
  
  psms_mods[[p]] = psms_mod_df
  psms_mod_mat[,p] = sapply(psms_mod_df$term[!psms_mod_df$term %in% "(Intercept)"], FUN = function(x){
    paste0(round(psms_mod_df[psms_mod_df$term %in% x, "estimate"], 1), " (", round(psms_mod_df[psms_mod_df$term %in% x, "conf.low"], 1), " - ", round(psms_mod_df[psms_mod_df$term %in% x, "conf.high"], 1), ")\nP = ", round(psms_mod_df[psms_mod_df$term %in% x, "p.value"], 3))
  })
}
colnames(psms_mod_mat) = sapply(colnames(psms_mod_mat), FUN = function(s){
  if (s == "psmstoil"){activity = "Toileting"}
  if (s == "psmsfeed"){activity = "Feeding"}
  if (s == "psmsdress"){activity = "Dressing"}
  if (s == "psmsgroom"){activity = "Grooming"}
  if (s == "psmsamb"){activity = "Ambulation"}
  if (s == "psmsbath"){activity = "Bathing"}
  activity
})

write.csv(file = "output/individual_psms.csv", x = psms_mod_mat)

