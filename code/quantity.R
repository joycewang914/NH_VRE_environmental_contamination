# Compare VRE colonization burden between patient-only and both samples using VRE plate quantity as proxy
library(viridis)
library(reshape2)

micro = read.csv("/Users/joycewang/Desktop/Mody_NH_collections/Pathways/CID_Pathways_Micro_Data.csv", header = T)
final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

vre_quantR = t(apply(final_phase2, 1, FUN = function(x){
  pt = x["studyID"]
  visit = gsub(" ", "", x["visit"])
  perirectal = gsub(" ", "", x["VREsiteR"])
  vre_col = gsub(" ", "", x['vre_col'])

  perirectal_VRE_quant = NA
  
  if (!is.na(perirectal) & perirectal > 0){
    perirectal_VRE_quant = micro[micro$Resident.Study.ID.Number %in% pt & micro$Visit.Number %in% visit & micro$Specimen.Type %in% "R", "VRE.1.Quantity"]}
  
  cbind(pt, vre_col,perirectal_VRE_quant)
}))

colnames(vre_quantR) = c("Patient_ID", "VRE_col", "Quant")

vre_quantR_df = as.data.frame(vre_quantR)

# Quantification by colonizatoin pattern
vre_col_quant = table(vre_quantR_df$VRE_col, vre_quantR_df$Quant)[3:4,]
vre_col_quant_prop = t(apply(vre_col_quant, 1, FUN = function(x){x/sum(x)}))
rownames(vre_col_quant_prop) = c("Patient-only", "Both")

file ="output/perirectal_VRE_quant_vs_colonization.png"
png(file)
par(mar = c(3, 8, 5, 3))
prop_bp_cols = viridis(4)
bp = barplot(t(vre_col_quant_prop*100), beside = T, ylab="Proportion (%)",  cex.names=1.5, las=1, ylim=c(0,60), xpd = T, col = prop_bp_cols, border = F, cex.axis = 1.5, cex.lab = 1.5)

legend("topright", legend = rownames(t(vre_col_quant_prop)), border = F, bty = "n", fill = prop_bp_cols, title = "Quandrant", horiz = F, adj = 0, cex = 1.5, inset = -.1, xpd = T)

text(x = melt(bp)[,3], y= c(t(vre_col_quant_prop)*100)+2.5, c(t(vre_col_quant)), xpd = T, cex = 1.5)

dev.off()

KW_test = kruskal.test(as.numeric(vre_quantR_df[vre_quantR_df$VRE_col %in% 2:3, "Quant"]) ~ vre_quantR_df[vre_quantR_df$VRE_col %in% 2:3, "VRE_col"])
