# Compare VRE colonization quantity between patient-only and both samples

micro = read.csv("/Users/joycewang/Desktop/Mody_NH_collections/Pathways/CID_Pathways_Micro_Data.csv", header = T)

vre_quantG_R = apply(final_phase2, 1, FUN = function(x){
  pt = x["studyID"]
  visit = gsub(" ", "", x["visit"])
  groin = gsub(" ", "", x["VREsiteG"])
  perirectal = gsub(" ", "", x["VREsiteR"])
  print(c(pt, visit, groin, perirectal))
  
  groin_VRE_quant = NA
  perirectal_VRE_quant = NA
  
  if (!is.na(groin) & groin > 0){
    groin_VRE_quant = micro[micro$Resident.Study.ID.Number %in% pt & micro$Visit.Number %in% visit & micro$Specimen.Type %in% "G", "VRE.1.Quantity"]}
  
  if (!is.na(perirectal) &perirectal > 0){
    perirectal_VRE_quant = micro[micro$Resident.Study.ID.Number %in% pt & micro$Visit.Number %in% visit & micro$Specimen.Type %in% "R", "VRE.1.Quantity"]}
  
  temp_dat = c(groin_VRE_quant, perirectal_VRE_quant)
  print(temp_dat)
  temp_dat
  
})

vre_quant_culture = cbind(final_phase2[,c("VREsiteG", "VREsiteR",keep_sites, "vre_env", "vre_col")], t(vre_quantG_R))

colnames(vre_quant_culture)[colnames(vre_quant_culture) %in% c("1", "2")] = c("G_quant", "R_quant")


# Quantification by colonizatoin pattern
vre_col_quant = table(vre_quant_culture[,"R_quant"], vre_quant_culture[,"vre_col"])[,3:4]
vre_col_quant_prop = apply(vre_col_quant, 2, FUN = function(x){x/sum(x)})
colnames(vre_col_quant_prop) = c("Patient-only", "Both")

# file ="/Users/joycewang/Desktop/Snitkin_lab/Manuscripts/Pathways_mbiome_environment/Draft/v2/figures/perirectal_VRE_quant_vs_colonization.pdf"
# pdf(file)
prop_bp_cols = viridis(4)[3:4]
bp = barplot(t(vre_col_quant_prop*100), beside = T, ylab="Proportion (%)", xlab = "VRE in quadrant (Perirectal samples)", cex.names=1.5, las=1, ylim=c(0,60), xpd = T, col = prop_bp_cols, border = F, cex.axis = 1.5, cex.lab = 1.5)

legend("topright", legend = rownames(t(vre_col_quant_prop)), border = F, bty = "n", fill = prop_bp_cols, title = "VRE colonization pattern", horiz = T, adj = 0, cex = 1.5)

text(x = melt(bp)[,3], y= c(t(vre_col_quant_prop)*100)+2.5, c(t(vre_col_quant)), xpd = T, cex = 1.5)
# dev.off()

prop.test(table(vre_quant_culture$R_quant, vre_quant_culture$vre_col)[,3:4], correct = F)

KW_test = kruskal.test(vre_quant_culture[vre_quant_culture$vre_col %in% 2:3, "vre_col"] ~ vre_quant_culture[vre_quant_culture$vre_col %in% 2:3, "R_quant"])