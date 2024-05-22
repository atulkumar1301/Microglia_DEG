#! /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/bin/Rscript
library(ROCR)
library(pROC)
library(bdpv)
library(data.table)
library(dplyr)
library(caret)
library(Rcpp)
library(RNOmni)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=10, nrow=7)) 
names(TABLE)<-c("Model", "Effect", "SE", "R2", "L95", "U95", "AIC", "BIC", "Number_of_Variants", "P")
Clinical_data <- fread (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/5_Data_Full_Imputed_Analysis_Taupet.txt"))
Full_File <- fread ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/CSF_pTau217/PCA_FILE.eigenvec") %>% merge(Clinical_data, by = 'IID')
PRS_1 <- fread ("p_value_0.05.sscore") %>% merge(Full_File, by = 'IID')
PRS_2 <- fread ("p_value_0.005.sscore") %>% merge(Full_File, by = 'IID')
PRS_3 <- fread ("p_value_0.0005.sscore") %>% merge(Full_File, by = 'IID')
PRS_4 <- fread ("p_value_5e-05.sscore") %>% merge(Full_File, by = 'IID')
PRS_5 <- fread ("p_value_5e-06.sscore") %>% merge(Full_File, by = 'IID')
PRS_6 <- fread ("p_value_5e-07.sscore") %>% merge(Full_File, by = 'IID')
PRS_7 <- fread ("p_value_5e-08.sscore") %>% merge(Full_File, by = 'IID')
df <- Clinical_data [, -1:-35]
for (p in colnames (df)) {
  i = 1
  data2 <- PRS_1
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 1"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS1 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS1)
  df1 <- df.residual (model_PRS1)
  
  TABLE[i, 2] <- summary (model_PRS1)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS1)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS1) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS1) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS1)
  TABLE[i, 8] <- BIC (model_PRS1)
  
  TABLE[i, 9] <- nrow (fread ("p_value_0.05.txt"))
  TABLE[i, 10] <- summary (model_PRS1)$coefficients[2, "Pr(>|t|)"]
  
  #=========================================================================================================
  
  i = 2
  
  data2 <- PRS_2
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 2"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS2 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS2)
  df1 <- df.residual (model_PRS2)
  
  TABLE[i, 2] <- summary (model_PRS2)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS2)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS2) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS2) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS2)
  TABLE[i, 8] <- BIC (model_PRS2)
  
  TABLE[i, 9] <- nrow (fread ("p_value_0.005.txt"))
  TABLE[i, 10] <- summary (model_PRS2)$coefficients[2, "Pr(>|t|)"]
  
  #================================================================================================
  
  i = 3
  
  data2 <- PRS_3
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 3"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS3 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS3)
  df1 <- df.residual (model_PRS3)
  
  TABLE[i, 2] <- summary (model_PRS3)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS3)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS3) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS3) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS3)
  TABLE[i, 8] <- BIC (model_PRS3)
  
  TABLE[i, 9] <- nrow (fread ("p_value_0.0005.txt"))
  TABLE[i, 10] <- summary (model_PRS3)$coefficients[2, "Pr(>|t|)"]
  
  #=============================================================================================
  
  i = 4
  
  data2 <- PRS_4
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 4"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS4 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS4)
  df1 <- df.residual (model_PRS4)
  
  TABLE[i, 2] <- summary (model_PRS4)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS4)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS4) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS4) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS4)
  TABLE[i, 8] <- BIC (model_PRS4)
  
  TABLE[i, 9] <- nrow (fread ("p_value_5e-05.txt"))
  TABLE[i, 10] <- summary (model_PRS4)$coefficients[2, "Pr(>|t|)"]
  
  #====================================================================================================
  
  i = 5
  
  data2 <- PRS_5
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 5"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS5 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS5)
  df1 <- df.residual (model_PRS5)
  
  TABLE[i, 2] <- summary (model_PRS5)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS5)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS5) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS5) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS5)
  TABLE[i, 8] <- BIC (model_PRS5)
  
  TABLE[i, 9] <- nrow (fread ("p_value_5e-06.txt"))
  TABLE[i, 10] <- summary (model_PRS5)$coefficients[2, "Pr(>|t|)"]
  
  #===================================================================================================
  
  i = 6
  
  data2 <- PRS_6
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 6"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS6 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS6)
  df1 <- df.residual (model_PRS6)
  
  TABLE[i, 2] <- summary (model_PRS6)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS6)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS6) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS6) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS6)
  TABLE[i, 8] <- BIC (model_PRS6)
  
  TABLE[i, 9] <- nrow (fread ("p_value_5e-07.txt"))
  TABLE[i, 10] <- summary (model_PRS6)$coefficients[2, "Pr(>|t|)"]
  
  #=========================================================================================================
  
  i = 7
  
  data2 <- PRS_7
  plot (density(data2[[p]]))
  data2$Normal <- RankNorm(data2[[p]])
  plot (density(data2$Normal))
  data2$PHENO<- data2$Normal
  TABLE [i, 1] <- "PRS 7"
  
  m1 <- mean (data2$SCORE1_AVG)
  sd1 <- sd (data2$SCORE1_AVG)
  data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1
  
  modeldata20 <- glm (data2$PHENO ~ 1, family = gaussian, data = data2)
  model_PRS6 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$E2 + data2$E4 + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10 + data2$Abnormal_CSF_Ab42_Ab40_Ratio, family = gaussian, data = data2)
  
  
  l0 <- deviance (modeldata20)
  df0 <- df.residual (modeldata20) #degree of freedom
  l1 <- deviance (model_PRS6)
  df1 <- df.residual (model_PRS6)
  
  TABLE[i, 2] <- summary (model_PRS6)$coefficients[2, "Estimate"]
  TABLE[i, 3] <- summary(model_PRS6)$coefficients[2, "Std. Error"]
  
  TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
  TABLE[i, 5] <- confint (model_PRS6) [2, 1] 
  TABLE[i, 6] <- confint (model_PRS6) [2, 2]
  TABLE[i, 7] <- AIC (model_PRS6)
  TABLE[i, 8] <- BIC (model_PRS6)
  
  
  TABLE[i, 9] <- nrow (fread ("p_value_5e-08.txt"))
  TABLE[i, 10] <- summary (model_PRS6)$coefficients[2, "Pr(>|t|)"]
  
  #=======================================================================================================
  
  TABLE$Bonferroni <- p.adjust(TABLE$P, method = "bonferroni", n = length(TABLE$P))
  TABLE$FDR <- p.adjust(TABLE$P, method = "fdr", n = length(TABLE$P))
  write.table (TABLE, file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/Microglia_DEG/BF_2/Tau_PET/Original_PRS/Non_APOE/Adjusted/PRS_Result_",p,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
