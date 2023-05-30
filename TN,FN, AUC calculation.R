suppressWarnings ( suppressPackageStartupMessages ( library ( data.table ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( dplyr ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( fst ))) 
suppressWarnings ( suppressPackageStartupMessages ( library ( ggplot2 )))
library(dndscv)


library(readxl)
all_mutated_genes <- read_excel("all mutated genes.xlsx")

########################### get baseline negatives ##########################
########### exclude genes that not appear in all_mutated_genes from all_known_drivers
result <- summary_All_known_drivers_1042_genes$gene_name %in% all_mutated_genes$gene_name
summary_All_known_drivers_TF <-  cbind(summary_All_known_drivers_1042_genes[, 1],result)

summary_All_known_drivers_Compatible <- summary_All_known_drivers_TF[grepl("T",summary_All_known_drivers_TF$"result5"),]
summary_All_known_drivers_Compatible <- summary_All_known_drivers_Compatible[,1]

################ all mutated genes subtract known drivers as true negatives
total_non_driver <- setdiff(all_mutated_genes, summary_All_known_drivers_Compatible)
total_non_driver <- total_non_driver[!grepl("CDKN2A.p14arf|CDKN2A.p16INK4a",total_non_driver$gene_name),]


################## number of TN,FN ################################
colnames(major_SBS)[1] <- "gene_name"
majorSBS_non_driver <- setdiff(all_mutated_genes, major_SBS)
colnames(histology)[1] <- "gene_name"
histology_non_driver <- setdiff(all_mutated_genes, histology)

result <- majorSBS_non_driver$gene_name %in% total_non_driver$gene_name
table(result)
result1 <- histology_non_driver$gene_name %in% total_non_driver$gene_name
table(result1)


################## calculate AUC ####################################
#对于A方法：
#TPR_A = TP_A / (TP_A + FN_A) = 0.76 / (0.76 + 0.049885378) ≈ 0.9385
#FPR_A = FP_A / (FP_A + TN_A) = 0.24 / (0.24 + 0.950114622) ≈ 0.2015

#对于B方法：
#TPR_B = TP_B / (TP_B + FN_B) = 0.8 / (0.8 + 0.049835543) ≈ 0.9414
#FPR_B = FP_B / (FP_B + TN_B) = 0.2 / (0.2 + 0.950164457) ≈ 0.1738

TPR_A <- 0.9385
FPR_A <- 0.2015
TPR_B <- 0.9414
FPR_B <- 0.1738

roc_data <- data.frame(FPR = c(0, FPR_A, 1), TPR = c(0, TPR_A, 1))
AUC_A <- auc(roc_data$FPR, roc_data$TPR)
print(AUC_A)

roc_data <- data.frame(FPR = c(0, FPR_B, 1), TPR = c(0, TPR_B, 1))
AUC_B <- auc(roc_data$FPR, roc_data$TPR)
print(AUC_B)


################# McNemar's test  ##############
A_TP <-0.76
A_FN <- 0.049885378
B_TP <- 0.8
B_FN <-0.049835543
data <- matrix(c(A_TP, B_TP, A_FN, B_FN), nrow = 2, byrow = TRUE)
colnames(data) <- c("A_Positive", "B_Positive")
rownames(data) <- c("A_Negative", "B_Negative")

result <- mcnemar.test(data)

print(result)
