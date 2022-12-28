### Generate barplot indicating VRE colonization on hand in different colonization patterns ###
library(viridis)

####
# RUN ANALYSIS
####
final_phase2 = readRDS("input_data/2020-12-01_phase2_db.RDS")

#vre_col = colonization on O, R, G, and N

file ="output/VRE_hand_colonization_pattern.pdf"
pdf(file, onefile = T)

par(mar = c(10,10, 10, 15))
hand_table = table(final_phase2$vre_col, final_phase2$VREsiteH)
hand_table_freq = hand_table[,"1"]/rowSums(hand_table)
vre_patterns = c("Uncolonized", "Environment-only", "Patient-only", "Both")
vre_pattern_labels = apply(hand_table, 1, FUN = function(x){paste0("(", x[2], " / ", sum(x), ")")})

bp = barplot(t(hand_table_freq*100),  
             names.arg = vre_patterns,
             xlab = "VRE colonization on hand (%)", 
             cex.lab = 1.5, 
             cex.names = 1.5, 
             xpd = T, 
             horiz = T, 
             col = viridis(4)[1:4], 
             border = F,  
             plot = T, 
             las = 1, 
             cex.axis = 1.5, 
             beside = T)
text(y = bp, x = hand_table_freq *100 - 2, labels = vre_pattern_labels, xpd = T, cex = 1.5, xpd = T, pos = 4)
text(y = bp, x = -52, vre_patterns, xpd = T, cex =1.5, adj = 0)

segments(x0 = t(hand_table_freq*100)[4] + 28.5,y0 = bp[4], x1 = t(hand_table_freq*100)[4] + 28.5, y1 = bp[3], xpd = T, lwd = 3)

p_val = 0.001
both_pt_p = fisher.test(hand_table[3:4,])
text(y = (bp[4] + bp[3])/2, x = t(hand_table_freq*100)[4] + 45,
     labels = ifelse(both_pt_p$p.value < p_val, paste0("P < ", p_val), ""), xpd = T, cex = 1.5)

dev.off()