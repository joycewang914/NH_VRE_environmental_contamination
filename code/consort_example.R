# Betsy suggested making a consort diagram to clarify data usage
# Use R package "consort" to generate plot
# https://cran.r-project.org/web/packages/consort/vignettes/consrot_diagram.html
library(consort)

# Example 
# set up mock dataframe
set.seed(1001)
N <- 300

trialno <- sample(c(1000:2000), N)
exc1 <- rep(NA, N)
exc1[sample(1:N, 15)] <- sample(c("Sample not collected", "MRI not collected",
                                  "Other"), 15, replace = T, prob = c(0.4, 0.4, 0.2))

induc <- rep(NA, N)
induc[is.na(exc1)] <- trialno[is.na(exc1)]

exc2 <- rep(NA, N)
exc2[sample(1:N, 20)] <- sample(c("Sample not collected", "Dead",
                                  "Other"), 20, replace = T, prob = c(0.4, 0.4, 0.2))
exc2[is.na(induc)] <- NA

exc <- ifelse(is.na(exc2), exc1, exc2)

arm <- rep(NA, N)
arm[is.na(exc)] <- sample(c("Conc", "Seq"), sum(is.na(exc)), replace = T)
arm3 <- sample(c("Trt A", "Trt B", "Trt C"), N, replace = T)
arm3[is.na(arm)] <- NA

fow1 <- rep(NA, N)
fow1[!is.na(arm)] <- sample(c("Withdraw", "Discontinued", "Death", "Other", NA),
                            sum(!is.na(arm)), replace = T, 
                            prob = c(0.05, 0.05, 0.05, 0.05, 0.8))
fow2 <- rep(NA, N)
fow2[!is.na(arm) & is.na(fow1)] <- sample(c("Protocol deviation", "Outcome missing", NA),
                                          sum(!is.na(arm) & is.na(fow1)), replace = T, 
                                          prob = c(0.05, 0.05, 0.9))


df <- data.frame(trialno, exc1, induc, exc2, exc, arm, arm3, fow1, fow2)
rm(trialno, exc1, induc, exc2, exc, arm, arm3, fow1, fow2, N)

# Plot single arm
out <- consort_plot(data = df,
                    orders = c(trialno = "Population",
                               exc1    = "Excluded",
                               arm     = "Allocated",
                               fow1    = "Lost of Follow-up",
                               trialno = "Finished Followup",
                               fow2    = "Not evaluable for the final analysis",
                               trialno = "Final Analysis"),
                    side_box = c("exc1", "fow1", "fow2"),
                    cex = 0.9)
plot(out)

png("../figures/consort_diagram.png", width = 29, 
    height = 21, res = 300, units = "cm", type = "cairo") 
plot(out)
dev.off() 