# Temporal trend of VRE colonization
library(ggalluvial)

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

# Alluvial plot of colonization pattern
col_mat = tapply(1:nrow(final_phase2), INDEX = list(final_phase2$studyID, final_phase2$visit), FUN = function(x){
  final_phase2[x,"vre_col"]
} )

alluvial_df = melt(col_mat)
colnames(alluvial_df) = c("Patient", "Visit", "VRE_colonization")
alluvial_df$freq = 1
alluvial_df$Patient = as.factor(alluvial_df$Patient)
alluvial_df$VRE_colonization_category = sapply(alluvial_df$VRE_colonization, FUN = function(x){
  col_pat = NA
  if (x %in% 0:3){
    if (x == 0){col_pat = "Uncolonized"}
    if (x == 1){col_pat = "Environment-only"}
    if (x == 2){col_pat = "Patient-only"}
    if (x == 3){col_pat = "Both environment and patient"}
  }
  col_pat
})
alluvial_df$Visit = as.factor(alluvial_df$Visit)
alluvial_df$VRE_colonization_category = factor(alluvial_df$VRE_colonization_category, levels = c("Uncolonized", "Environment-only", "Patient-only", "Both environment and patient"))


ggplot(alluvial_df,
       aes(x = Visit, stratum = VRE_colonization_category, alluvium = Patient,
           y = freq, 
           fill = VRE_colonization_category, label = VRE_colonization)) +
  ylab("Number of Patients") + 
  labs(fill = "") + 
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5, lwd = .5) +
  theme(legend.position = "right",
        text = element_text(size=30), 
        legend.text=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"))



ggsave("output/alluvial_trend.png", width=20)
