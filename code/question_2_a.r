library(dplyr)

####################################
# Read coverage of 
####################################

# Load the dataset
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)

# Remove gene column
original_gene_counts$Gene <- NULL

# Sum the counts for each sample
counts_per_sample <- c(colSums(original_gene_counts))

# Histogram of read coverage
hist(counts_per_sample, col="darkgreen", xlab="Read counts", main="Distribution of read coverage")


####################################
# Read coverage of all patients 
# (Cumulative distribution for each patient through genes)
####################################

library(ggplot2)

# load dataset
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)

# remove gene name column
original_gene_counts$Gene <- NULL

# Adds mean counts per gene to dataset
original_gene_counts$GenExpression_Mean=rowMeans(original_gene_counts)

# sort all patients
patient_counts_sorted <- apply(original_gene_counts,2,sort,decreasing=T)

# relative cumsum calculus
check<-function(x){ return((cumsum(x)/sum(x))*100) }

# apply relative cumsum 
cumsum_patients <- apply(patient_counts_sorted,2,check)
cumsum_patients <- data.frame(list(cumsum_patients))

# Adds mean counts per gene to dataset
cumsum_patients$Cumulative_Mean <- rowMeans(cumsum_patients)

plot_all <- ggplot(cumsum_patients, aes(x=as.numeric(row.names(cumsum_patients))))

for (i in 1:878)
{
  plot_all <- plot_all + geom_line(aes_string(y = cumsum_patients[,i]),color = "red")
  print(i)
}

# Plot the genetic expression mean - Black
# plot_all <- plot_all + geom_line(aes_string(y = cumsum_patients[,879]),color = "black")
# Plot the cumulative mean - White
plot_all <- plot_all + geom_line(aes_string(y = cumsum_patients[,880]),color = "white")

plot_all <- plot_all + labs(x = "Number of genes", y= "Relative cumsum (%)")

print(plot_all)

##########################################################
# Distinct Genes (Reads Sequenced) 
#########################################################

library(ggplot2)

# load dataset
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)
# remove gene name column
original_gene_counts$Gene <- NULL

# Load the dataset
clinical_annotation<- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE)

# sort all patients
patient_counts_sorted <- apply(original_gene_counts,2,sort,decreasing=T)
# relative cumsum calculus
check<-function(x){ return((cumsum(x)/sum(x))*100) }
# apply relative cumsum 
cumsum_patients <- apply(patient_counts_sorted,2,check)
cumsum_patients <- data.frame(list(cumsum_patients))

# Finding the max outlier
num_of_line <- 125
line_find <- cumsum_patients[num_of_line,] 
index_max_outlier <- which(line_find == max(line_find))
max_outlier = patient_counts_sorted[,index_max_outlier]
total_counts_max_outlier = colSums(data.frame(max_outlier))

# Using Mean individual
mean_outlier <- rowMeans(cumsum_patients) 
prob_density <- seq(1,20502)
i_set = seq(1, 20502)
for (i in i_set){
  if (i == 1 ){
    prob_density[i] <- mean_outlier[i]/100
  } else {
    prob_density[i] <- (mean_outlier[i+1] - mean_outlier[i])/100
  }
}
prob_density <- prob_density
mean_outlier = prob_density*total_counts_max_outlier
mean_outlier = round(mean_outlier)
                     
# Building repeated vector
num_genes_considered = 20502
times_mean_outlier = mean_outlier[1:num_genes_considered]
rep_mean_outlier = rep(1:num_genes_considered, times_mean_outlier)
times_max_outlier = max_outlier[1:num_genes_considered]
rep_max_outlier = rep(1:num_genes_considered, times_max_outlier)

# Sampling the repeated vector
sup_num_sampled_reads = 100000
num_diff_genes_mean_vec = list()
num_diff_genes_max_vec = list()
for (i in 1:sup_num_sampled_reads)
{
  if (i %% 1000 == 0){
    print(i)
  }
  sample_mean_outlier = sample(rep_mean_outlier,i)
  unique_mean_outlier = unique(sample_mean_outlier)
  num_diff_genes_mean = length(unique_mean_outlier)
  num_diff_genes_mean_vec = append(num_diff_genes_mean_vec, num_diff_genes_mean)
  
  sample_max_outlier = sample(rep_max_outlier,i)
  unique_max_outlier = unique(sample_max_outlier)
  num_diff_genes_max = length(unique_max_outlier)
  num_diff_genes_max_vec = append(num_diff_genes_max_vec, num_diff_genes_max)
}

num_diff_genes <- data.frame("mean" = t(data.frame(num_diff_genes_mean_vec)), "max" = t(data.frame(num_diff_genes_max_vec)))

# Plot the number of different genes per number of sampled reads
plot_samples <- ggplot(num_diff_genes,aes(x = 1:sup_num_sampled_reads))
# Max line -> Red
plot_samples <- plot_samples + geom_line(aes_string(y = num_diff_genes[,2]),color = "red")
# Mean line -> Blue
plot_samples <- plot_samples + geom_line(aes_string(y = num_diff_genes[,1]),color = "blue")
plot_samples <- plot_samples + labs(x = "Number of Sampled Reads", y= "Number of Distinct Genes")
print(plot_samples)
