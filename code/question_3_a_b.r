setwd("~/+++ Faculdade +++/5ano_1sem/Bioinformática/Labs/Lab 6 (omics)/")
rm(list=ls()) #clean variables in the local environment

#Required libraries
library(tidyverse)
library(edgeR)
library(limma)
library(pROC)
library(mlr)
library(survival)
library(survminer)
library(dplyr)
library(glmnet)
library(naivebayes)

#### Group 3 (a) ####
# Read the the counts as matrix
read_counts <- as.matrix(read.table(file = 'TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE, row.names="Gene"))
samples_info <- as.matrix(read.table(file = 'TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE,row.names="Patient.ID"))

# Normalization of counts (for each samples) + dealing with bias
# Requires edgeR and limma
# Obtain factor for normalization with calcNormFactors
dge <- DGEList(counts=read_counts)
dge <- calcNormFactors(dge)
# Apply normalization factors and convert to log2 counts per million reads (CPM)
v <- voom(dge)
read_counts_norm <- v$E
 
# Select genes to analyze
genes_name <- c("ESR1", "ESR2", "PGR", "ERBB2")
read_counts_selection <- read_counts_norm[genes_name,]

# Classifications of the immunohistochemistry-based tests
#Check number of samples with information
n_samples <- length(samples_info[,1])
#Convert classification to NA- NA, unkonwn or indeterminate
protein_expression <- c("Estrogen.Receptor","Progesterone.Receptor", "HER2")
classification <- samples_info[,protein_expression]

i <- 1
while (i<=n_samples)
{
  
  if (!is.na(classification[i,1])){
    #Check Estrogen Receptor
    if (classification[i,1]=="positive"){
      classification[i,1] <- 1} else {
        if (classification[i,1]=="negative"){
          classification[i,1]<-0} else{
            classification[i,1] <- NA
          }
      }
  }
  
  if (!is.na(classification[i,2])){
    #Check Progesterone Receptor
    if (classification[i,2]=="positive"){
      classification[i,2] <- 1} else {
        if (classification[i,2]=="negative"){
          classification[i,2]<-0} else{
            classification[i,2] <- NA
          }
      }
  }
  
  if (!is.na(classification[i,3])){
    #Check HER2 protein
    if (classification[i,3]=="positive"){
      classification[i,3] <- 1} else {
        if (classification[i,3]=="negative"){
          classification[i,3]<-0} else{
            classification[i,3] <- NA
          }
      }
  }


  i=i+1
  
}


# Select tumour samples from read counts
# Requires tidyverse library
read_counts_patients_tb <- as_tibble(read_counts_selection,rownames = NA)
read_counts_patients <- read_counts_patients_tb %>% select(ends_with("01"))
n_samples_reads <- length(read_counts_patients[1,])
patients_id <- colnames(read_counts_patients)
read_counts_patients <- t(read_counts_patients)

# Select only samples present in both counts and sample information
patients_id <- gsub("[.]", "-", patients_id)
patients_id <- gsub("-01", "", patients_id)
classification <- classification[patients_id,]


# Divide in training and testing set
# 75% of the sample size for training
training_size <- floor(0.75 * n_samples_reads)
testing_size <- n_samples_reads - training_size

#Setting the seed makes the partitions reproducible
set.seed(123)
train_ind <- sample(seq_len(n_samples_reads), size = training_size)

train_x <- read_counts_patients[train_ind,]
train_y <- classification[train_ind,]

test_x <- read_counts_patients[-train_ind,]
test_y <- classification[-train_ind,]

#Train the classifier
#Estrogen 1
to_select = !is.na(train_y[,1])
train_x_es_1 = train_x[to_select,1]
train_y_es_1 = train_y[to_select,1]
#Plot distribution
dev.new()
plot(x=train_x_es_1, y=train_y_es_1,xlab="ESR1 Normalized counts", ylab="Classification")
#dev.copy(png,'Estrogen_1.png')
dev.off()
#Plot ROC curve
dev.new()
Est1 <- roc(train_y_es_1,train_x_es_1,
    plot=TRUE,
    legacy.axes=TRUE,
    xlab="False Positive Rate",
    ylab="True Positive Rate",
    lwd=4)
auc_Est1 <- round(Est1$auc, digits=2)
text(x=0.2,y=0.2,labels=paste("AUC=",auc_Est1))
par(pty="s")
distEst1 <- sqrt((1-Est1$specificities)^2 + (Est1$sensitivities-1)^2)
idxEst1 <- which.min(distEst1)
spec_minEst1 <- round(Est1$specificities[idxEst1], digits=2)
sens_minEst1 <- round(Est1$sensitivities[idxEst1], digits=2)
points(x=spec_minEst1,y=sens_minEst1, pch=19,cex=2,col="red")
text(spec_minEst1-0.15, sens_minEst1-0.05, paste("(",spec_minEst1,",",sens_minEst1,")"))
dev.copy(png,'Estrogen_1_ROC.png')


dev.off()
#Check thresholds
roc.info <- roc(train_y_es_1,train_x_es_1,
                legacy.axes=TRUE)
roc.df <- data.frame(
  tpr = roc.info$sensitivities,
  fpr = (1-roc.info$specificities),
  thresholds = roc.info$thresholds)

thresholds_es_1 = roc.df[roc.df$tpr>0.90,]
threshold_es_1 = thresholds_es_1[131,3] #by inspection

#Check performance
to_select = !is.na(test_y[,1])
test_x_es_1 = test_x[to_select,1]
test_y_es_1 = test_y[to_select,1]

TP_FN = sum(as.numeric(test_y_es_1)) #TP and FN
TP = sum(as.numeric(test_y_es_1[test_x_es_1>threshold_es_1])) #TP
true_positive_rate_es_1 <- TP/TP_FN #TPR (sensitivity)

FP = sum(as.numeric(test_y_es_1[test_x_es_1<=threshold_es_1])) #FP
FP_TN = length(test_y_es_1) #FP and TN
false_positive_rate_es_1 <- FP/FP_TN


#Estrogen 2
to_select = !is.na(train_y[,1])
train_x_es_2 = train_x[to_select,2]
train_y_es_2 = train_y[to_select,1]
#Plot distribution
dev.new()
plot(x=train_x_es_2, y=train_y_es_2,xlab="ESR2 Normalized counts", ylab="Classification")
#dev.copy(png,'Estrogen_2.png')
dev.off()
#Plot ROC curve
dev.new()
Est2 <- roc(train_y_es_2,train_x_es_2,
    plot=TRUE,
    legacy.axes=TRUE,
    xlab="False Positive Rate",
    ylab="True Positive Rate",
    lwd=4)
auc_Est2 <- round(Est2$auc, digits=2)
text(x=0.2,y=0.2,labels=paste("AUC=",auc_Est2))
par(pty="s")
distEst2 <- sqrt((1-Est2$specificities)^2 + (Est2$sensitivities-1)^2)
idxEst2 <- which.min(distEst1)
spec_minEst2 <- round(Est2$specificities[idxEst2], digits=2)
sens_minEst2 <- round(Est2$sensitivities[idxEst2], digits=2)
points(x=spec_minEst2,y=sens_minEst2, pch=19,cex=2,col="red")
text(spec_minEst2+0.05, sens_minEst2+0.05, paste("(",spec_minEst2,",",sens_minEst2,")"))
dev.copy(png,'Estrogen_2_ROC.png')
dev.off()


#Progesterone
to_select = !is.na(train_y[,2])
train_x_pr = train_x[to_select,3]
train_y_pr = train_y[to_select,2]
#Plot distribution
dev.new()
plot(x=train_x_es_1, y=train_y_es_1,xlab="PGR Normalized counts", ylab="Classification")
#dev.copy(png,'Progesterone.png')
dev.off()
#Plot ROC curve
dev.new()
Prog <- roc(train_y_pr,train_x_pr,
    plot=TRUE,
    legacy.axes=TRUE,
    xlab="False Positive Rate",
    ylab="True Positive Rate",
    lwd=4)
par(pty="s")
auc_Prog <- round(Prog$auc, digits=2)
text(x=0.2,y=0.2,labels=paste("AUC=",auc_Prog))
distProg <- sqrt((1-Prog$specificities)^2 + (Prog$sensitivities-1)^2)
idxProg <- which.min(distProg)
spec_minProg <- round(Prog$specificities[idxProg], digits=2)
sens_minProg <- round(Prog$sensitivities[idxProg], digits=2)
points(x=spec_minProg,y=sens_minProg, pch=19,cex=2,col="red")
text(spec_minProg-0.15, sens_minProg-0.05, paste("(",spec_minProg,",",sens_minProg,")"))
dev.copy(png,'Progesterone_ROC.png')
dev.off()

#Check thresholds
roc.info <- roc(train_y_pr,train_x_pr,
                legacy.axes=TRUE)
roc.df <- data.frame(
  tpr = roc.info$sensitivities,
  fpr = (1-roc.info$specificities),
  thresholds = roc.info$thresholds)

thresholds_pr = roc.df[roc.df$tpr>0.90,]
threshold_pr = thresholds_pr[154,3] #by inspection

#Check performance
to_select = !is.na(test_y[,2])
test_x_pr = test_x[to_select,3]
test_y_pr = test_y[to_select,2]

TP_FN = sum(as.numeric(test_y_pr)) #TP and FN
TP = sum(as.numeric(test_y_pr[test_x_pr>threshold_pr])) #TP
true_positive_rate_pr <- TP/TP_FN #TPR (sensitivity)

FP = sum(as.numeric(test_y_pr[test_x_pr<=threshold_pr])) #FP
FP_TN = length(test_y_pr) #FP and TN
false_positive_rate_pr <- FP/FP_TN


#HER2
to_select = !is.na(train_y[,3])
train_x_her = train_x[to_select,4]
train_y_her = train_y[to_select,3]
#Plot distribution
dev.new()
plot(x=train_x_her, y=train_y_her,xlab="ERBB2 Normalized counts", ylab="Classification")
#dev.copy(png,'HER2.png')
dev.off()
#Plot ROC curve
dev.new()
HER2 <- roc(train_y_her,train_x_her,
    plot=TRUE,
    legacy.axes=TRUE,
    xlab="False Positive Rate",
    ylab="True Positive Rate",
    lwd=4)
auc_HER2 <- round(HER2$auc, digits=2)
text(x=0.2,y=0.2,labels=paste("AUC=",auc_HER2))
par(pty="s")
distHER2 <- sqrt((1-HER2$specificities)^2 + (HER2$sensitivities-1)^2)
idxHER2 <- which.min(distHER2)
spec_minHER2 <- round(HER2$specificities[idxHER2], digits=2)
sens_minHER2 <- round(HER2$sensitivities[idxHER2], digits=2)
points(x=spec_minHER2,y=sens_minHER2, pch=19,cex=2,col="red")
text(spec_minHER2-0.15, sens_minHER2-0.05, paste("(",spec_minHER2,",",sens_minHER2,")"))
dev.copy(png,'HER2_ROC.png')
dev.off()



#Check thresholds
roc.info <- roc(train_y_her,train_x_her,
                legacy.axes=TRUE)
roc.df <- data.frame(
  tpr = roc.info$sensitivities,
  fpr = (1-roc.info$specificities),
  thresholds = roc.info$thresholds)

thresholds_her = roc.df[roc.df$tpr>0.90,]
threshold_her = thresholds_her[87,3] #by inspection

#Check performance
to_select = !is.na(test_y[,3])
test_x_her = test_x[to_select,4]
test_y_her = test_y[to_select,3]

TP_FN = sum(as.numeric(test_y_her)) #TP and FN
TP = sum(as.numeric(test_y_her[test_x_her>threshold_her])) #TP
true_positive_rate_her <- TP/TP_FN #TPR (sensitivity)

FP = sum(as.numeric(test_y_her[test_x_her<=threshold_her])) #FP
FP_TN = length(test_y_her) #FP and TN
false_positive_rate_her <- FP/FP_TN

#### Group 3 (b) ####
# Survival analysis
#Select PAM50 and survival columns
PAM50_Survival <- c("PAM50","Vital.status","Days.to.last.follow.up","Days.to.death")
survival_analysis_samples <- samples_info[,PAM50_Survival]
to_select = !is.na(survival_analysis_samples[,1])
survival_analysis_samples <- as.data.frame(survival_analysis_samples[to_select,],
                                           stringsAsFactors = FALSE)

#Generate life table
#Sort times & generate vector of events and censored data
survival_analysis_samples <- survival_analysis_samples[order(survival_analysis_samples$Days.to.last.follow.up,
                                        survival_analysis_samples$Days.to.death),]

n_alive = length(survival_analysis_samples[!is.na(survival_analysis_samples[,3]),1])
n_dead = length(survival_analysis_samples[,1]) - n_alive

survival_analysis_samples$futime <- c(survival_analysis_samples[1:n_alive,3],
                                      survival_analysis_samples[(n_alive+1):(n_dead+n_alive),4])

survival_analysis_samples$fustat <- c(rep(0, n_alive), 
                                      rep(1, n_dead))

survival_analysis_samples <- survival_analysis_samples[order(survival_analysis_samples$futime),]


#Generate survival object
surv_object <- Surv(time = as.numeric(survival_analysis_samples$futime), 
                    event = as.numeric(survival_analysis_samples$fustat))

#Fit Kaplan-Meier curves
fit_PAM50 <- survfit(surv_object ~ survival_analysis_samples$PAM50, data = survival_analysis_samples)

#Plot
dev.new()
ggsurvplot(fit_PAM50, data = survival_analysis_samples, 
           pval = TRUE,
           risk.table = TRUE,
           main = "Kaplan-Meier Survival Curve",
           submain = "Kaplan-Meier Survival Curves",
           xlab = "Time (days)", 
           ylab = "Survival probabilities")
dev.copy(png,'Survival_Analysis.png')
dev.off()

#Show median survival times
fit_PAM50

#### Group 3 (c) ####
#Select PAM50_genes
PAM50_genes <- c("UBE2T", "BIRC5", "NUF2","CDC6", "CCNB1", "TYMS", "MYBL2", "CEP55", "MELK", "NDC80",
                 "RRM2", "UBE2C", "CENPF", "PTTG1", "EXO1", "ORC6L", "ANLN", "CCNE1", "CDC20", "MKI67",
                 "KIF2C", "ACTR3B", "MYC", "EGFR", "KRT5", "PHGDH", "CDH3", "MIA", "KRT17", "FOXC1", 
                 "SFRP1", "KRT14", "ESR1", "SLC39A6", "BAG1", "MAPT", "PGR", "CXXC5", "MLPH", "BCL2",
                 "MDM2", "NAT1", "FOXA1", "BLVRA", "MMP11", "GPR160", "FGFR4", "GRB7", "TMEM45B","ERBB2")

# Select tumour samples from read counts
# Requires tidyverse library
read_counts_patients_tb <- as_tibble(read_counts_norm,rownames = NA)
read_counts_patients <- read_counts_patients_tb %>% select(ends_with("01"))
n_samples_reads <- length(read_counts_patients[1,])
patients_id <- colnames(read_counts_patients)
read_counts_patients <- t(read_counts_patients)

# Select only samples present in both counts and sample information
patients_id <- gsub("[.]", "-", patients_id)
patients_id <- gsub("-01", "", patients_id)
PAM_50_subtypes <- samples_info[patients_id,]
PAM_50_subtypes <- PAM_50_subtypes[,"PAM50"]

to_select = !is.na(PAM_50_subtypes)
read_counts_patients <- read_counts_patients[to_select,]
PAM_50_subtypes <- PAM_50_subtypes[to_select]

n_samples_reads_new = length(PAM_50_subtypes)

PAM_50_subtypes <- as.factor(PAM_50_subtypes)

#Select best set of genes using Lasso
#Select lambda by 10-folds cross validation
cv.regularizationmodel <- cv.glmnet(x = read_counts_patients,
                                    y = PAM_50_subtypes,
                                    family="multinomial",
                                    #type.multinomial="grouped",
                                    standardize = TRUE,
                                    alpha = 1.0,
                                    nfolds=10,
                                    parallel=TRUE)

idealLambda <- cv.regularizationmodel$lambda.1se

co <- coef(cv.regularizationmodel, s=idealLambda)
genes_LuminalA_i <- co[["Luminal A"]]@i[-1]
genes_LuminalB_i <- co[["Luminal B"]]@i[-1]
genes_BasalLike_i <- co[["Basal-like"]]@i[-1]
genes_HER2Enriched_i <- co[["HER2-enriched"]]@i[-1]
genes_NormalLike_i <- co[["Normal-like"]]@i[-1]

genes_id <- colnames(read_counts_patients)

genes_LuminalA <- genes_id[genes_LuminalA_i]
genes_LuminalB <- genes_id[genes_LuminalB_i]
genes_BasalLike <- genes_id[genes_BasalLike_i]
genes_HER2Enriched <- genes_id[genes_HER2Enriched_i]
genes_NormalLike <- genes_id[genes_NormalLike_i]

signature <- c(genes_LuminalA,
               genes_LuminalB,
               genes_BasalLike,
               genes_HER2Enriched,
               genes_NormalLike)
signature <- unique(unlist(strsplit(signature, " ")))

counter=0
j=1

while(j<=50){
  i=1
  while (i<=length(signature)){
    if(PAM50_genes[j]==signature[i])
    {counter = counter +1
    }
    i=i+1
  }
  j=j+1
}


#Compare performance using a Naive Bayes Classifier
#Divide in trainig and testing sets
# 75% of the sample size for training
training_size <- floor(0.75 * n_samples_reads_new)
testing_size <- n_samples_reads_new - training_size

#Setting the seed makes the partitions reproducible
set.seed(123)
train_ind <- sample(seq_len(n_samples_reads_new), size = training_size)

#PAM50
train_x <- read_counts_patients[train_ind,PAM50_genes]
train_y <- PAM_50_subtypes[train_ind]

test_x <- read_counts_patients[-train_ind,PAM50_genes]
test_y <- PAM_50_subtypes[-train_ind]

NB_PAM <- naive_bayes(x = train_x,
                      y = train_y, 
                      prior = NULL, 
                      laplace = 0.5)

prediction <- predict(NB_PAM, newdata = test_x)
CM_PAM <- as.matrix(table(Actual = test_y, Predicted = prediction))

#Signature
train_x <- read_counts_patients[train_ind,signature]
train_y <- PAM_50_subtypes[train_ind]

test_x <- read_counts_patients[-train_ind,signature]
test_y <- PAM_50_subtypes[-train_ind]

NB_SIG <- naive_bayes(x = train_x,
                      y = train_y, 
                      prior = NULL, 
                      laplace = 0.5)

prediction <- predict(NB_SIG, newdata = test_x)
CM_SIG <- as.matrix(table(Actual = test_y, Predicted = prediction))


#### Group 3 (d) ####
#Apply Naive Bayes for signature to NA samples
read_counts_patients_tb <- as_tibble(read_counts_norm,rownames = NA)
read_counts_patients <- read_counts_patients_tb %>% select(ends_with("01"))
n_samples_reads <- length(read_counts_patients[1,])

patients_id <- colnames(read_counts_patients)
patients_id <- gsub("[.]", "-", patients_id)
patients_id <- gsub("-01", "", patients_id)

read_counts_patients <- t(read_counts_patients)

NA_samples <- read_counts_patients[!to_select,signature]

prediction_NA <- predict(NB_SIG, newdata = NA_samples)

#Check suvival analysis and compare
#Survival analysis
#Select PAM50 and survival columns
PAM50_Survival <- c("PAM50","Vital.status","Days.to.last.follow.up","Days.to.death")
survival_analysis_samples <- samples_info[patients_id,PAM50_Survival]
survival_analysis_samples <- survival_analysis_samples[!to_select,PAM50_Survival]
survival_analysis_samples[,1] = as.character(prediction_NA)

survival_analysis_samples <- as.data.frame(survival_analysis_samples,
                                           stringsAsFactors = FALSE)

#Generate life table
#Sort times & generate vector of events and censored data
survival_analysis_samples <- survival_analysis_samples[order(survival_analysis_samples$Days.to.last.follow.up,
                                                             survival_analysis_samples$Days.to.death),]

n_alive = length(survival_analysis_samples[!is.na(survival_analysis_samples[,3]),1])
n_dead = length(survival_analysis_samples[,1]) - n_alive

survival_analysis_samples$futime <- c(survival_analysis_samples[1:n_alive,3],
                                      survival_analysis_samples[(n_alive+1):(n_dead+n_alive),4])

survival_analysis_samples$fustat <- c(rep(0, n_alive), 
                                      rep(1, n_dead))

survival_analysis_samples <- survival_analysis_samples[order(survival_analysis_samples$futime),]


#Generate survival object
surv_object <- Surv(time = as.numeric(survival_analysis_samples$futime), 
                    event = as.numeric(survival_analysis_samples$fustat))

#Fit Kaplan-Meier curves
fit_PAM50 <- survfit(surv_object ~ survival_analysis_samples$PAM50, data = survival_analysis_samples)

#Plot
dev.new()
ggsurvplot(fit_PAM50, data = survival_analysis_samples, 
           pval = TRUE,
           risk.table = TRUE,
           main = "Kaplan-Meier Survival Curve",
           submain = "Kaplan-Meier Survival Curves",
           xlab = "Time (days)", 
           ylab = "Survival probabilities")
dev.copy(png,'Survival_Analysis_d.png')
dev.off()