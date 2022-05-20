library(tibble)
library(dplyr)
library(survival)
library(survminer)

# Read the the counts as matrix
read_counts <- as.matrix(read.table(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE, row.names="Gene"))
samples_info <- as.matrix(read.table(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE,row.names="Patient.ID"))

# Normalization of counts (for each samples) + dealing with bias
# Requires edgeR and limma
# Obtain factor for normalization with calcNormFactors
dge <- DGEList(counts=read_counts)
dge <- calcNormFactors(dge)
# Apply normalization factors and convert to log2 counts per million reads (CPM)
v <- voom(dge)
read_counts_norm <- v$E



#### Group 3 (d) ####
#Apply Naive Bayes for signature to NA samples
read_counts_patients_tb <- as_tibble(read_counts_norm,rownames = NA)
read_counts_patients <- read_counts_patients_tb %>% select(ends_with("01"))
n_samples_reads <- length(read_counts_patients[1,])

patients_id <- colnames(read_counts_patients)
patients_id <- gsub("[.]", "-", patients_id)
patients_id <- gsub("-01", "", patients_id)

read_counts_patients <- t(read_counts_patients)

PAM_50_subtypes <- samples_info[patients_id,]
PAM_50_subtypes <- PAM_50_subtypes[,"PAM50"]

to_select = !is.na(PAM_50_subtypes)


NA_samples <- read_counts_patients[!is.na(PAM_50_subtypes),signature]

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
