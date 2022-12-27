# Descriptive stats on the relationsihp between microbial diveristy and antibioitc exposure

library(vioplot)

final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

abx_div = t.test(final_phase2$invsimp ~ as.factor(final_phase2$abx_30d_prior > 0))

table(final_phase2$abx_30d_prior > 0, final_phase2$vre_col)

fisher.test(table(final_phase2$abx_30d_prior > 0, final_phase2$vre_col)[,3:4])
