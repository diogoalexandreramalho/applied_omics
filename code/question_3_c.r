library(glmnet)
library(naivebayes)
library(caret)
library(ggplot2)
library(lattice)


clinical_annotation <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE)
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)


names <- names(original_gene_counts[-1])
patient_ids <- gsub("\\.", '-', substr(names,1,nchar(names)-3))
id_samples <- data.frame("Patient.id" = patient_ids, "Sample" = names(original_gene_counts[-1]))

# adds column with samples names
x <- setNames(data.frame(t(original_gene_counts[,-1])), original_gene_counts[,1])
x <- data.frame(id_samples$Sample, x)
names(x)[1] <- "Sample"



substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

check_tumor_state <- function(names){ 
  get_type_str = substrRight(names,2)
  get_type_int = as.numeric(get_type_str)
  
  if (get_type_int <= 5 ){
    state = 1
  } else { 
    state = 0
  }
  return(state) 
}

# adds pam50 results column
pam50 <- clinical_annotation[match(id_samples$Patient.id, clinical_annotation$Patient.ID), "PAM50"]
x$pam50 <- pam50

# removes normal samples and keeps tumour ones
names_vec <- data.frame("Sample" = list(id_samples$Sample))
Tumor_binary <- apply(names_vec,1,check_tumor_state)
x <- data.frame(Tumor_binary, x)
x <- x[x$Tumor_binary==1,]

x$Tumor_binary <- NULL
x$Sample <- NULL

# removes samples that don't have a PAM50 result
x <- x[!is.na(x$pam50),]

y <- matrix(x$pam50)

x$pam50 <- NULL
x <- as.matrix(x)

names(y)[1] <- "PAM50"
val <- factor(y,labels = c(1,2,3,4,5), levels = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like"))
val <- data.frame(val)



#Select best set of genes using Lasso
#Select lambda by 10-folds cross validation
cv.regularizationmodel <- cv.glmnet(x,
                                    y,
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

genes_id <- colnames(x)

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


# How many genes from our set are in PAM50
pam50_genes <- read.csv(file = '~/Downloads/Lab-6-Omics/pam50_genes.txt', sep = '\t')

counter = 0

for(gene in signature) {
  if(gene %in% pam50_genes) {
    counter = counter + 1
  }
}

print(counter)


#Compare performance using a Naive Bayes Classifier
#Divide in trainig and testing sets
# 75% of the sample size for training
nr_observations <- dim(x)[1]
training_size <- floor(0.75 * nr_observations)
testing_size <- nr_observations - training_size

#Setting the seed makes the partitions reproducible
set.seed(123)
train_ind <- sample(seq_len(nr_observations), size = training_size)

#PAM50
train_x <- x[train_ind,pam50_genes$Genes]
train_y <- y[train_ind]

test_x <- x[-train_ind,pam50_genes$Genes]
test_y <- y[-train_ind]

NB_PAM <- naive_bayes(x = train_x,
                      y = train_y, 
                      prior = NULL, 
                      laplace = 0.5)

prediction <- predict(NB_PAM, newdata = test_x)
CM_PAM <- as.matrix(table(Actual = test_y, Predicted = prediction))

#Signature
train_x <- x[train_ind,signature]
train_y <- y[train_ind]

test_x <- x[-train_ind,signature]
test_y <- y[-train_ind]

NB_SIG <- naive_bayes(x = train_x,
                      y = train_y, 
                      prior = NULL, 
                      laplace = 0.5)

prediction <- predict(NB_SIG, newdata = test_x)
CM_SIG <- as.matrix(table(Actual = test_y, Predicted = prediction))










