library(limma)
library(edgeR)
library(ggplot2)


# loads dataset
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)
original_gene_counts$Gene <- NULL


#####################
# NORMALIZE with voom
#####################
dge <- DGEList(counts=data.matrix(original_gene_counts))
dge <- calcNormFactors(dge)
v <- voom(dge, plot=TRUE)




###############
# BOXPLOTS
###############

# plot non-normalized dataset
ten_patient_counts <- original_gene_counts[,1:9]
ten_patient_counts$TCGA.GI.A2C8.11 <- original_gene_counts[,870]
ten_patient_counts[ten_patient_counts == 0] <- NA
boxplot(ten_patient_counts, log='y', pch=20, xlab = 'Samples', ylab = 'log(Read Counts)', main="Before normalization")


# plot normalized dataset
norm_ten_patient_counts <- data.frame(v$E[,1:9])
norm_ten_patient_counts$TCGA.GI.A2C8.11 <- v$E[,870]
boxplot(norm_ten_patient_counts, pch=20, xlab = 'Samples', ylab = 'log(Read Counts)', main="After normalization")



####################################################################
#               Variance per PC
####################################################################


# pca normalized
pca <- prcomp(t(v$E))

# plot the variance for first 10 PC components in descending order
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 2)
ylim <- c(0, 1.1*max(pca.var.per[1:10]))
xx <- barplot(pca.var.per[1:10], main="Variance in first 10 PC", ylim = ylim, xlab = "Principal Component", ylab = "Variance (%)")
freqs <- as.numeric(pca.var.per[1:10])
text(x = xx, y = freqs, label = freqs, pos = 3, cex = 0.8, col = "black")
axis(1, at=xx, labels=names(pca_components)[1:10], tick=FALSE, las=2, line=-0.5, cex.axis=0.8)


# plot the explained variance with PC increase
var_cumsum <- cumsum(pca.var.per)
plot(var_cumsum, type = "l", xlab = "Number of components", ylab = "Variance (%)", pch=20, ylim = c(0,100), main = "Cumulative variance with number of components")

#####################################################################################


####################################################################
#               Gene variance
####################################################################

# Genes associated to higer variance in PC1
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_scores_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_scores_ranked[1:10])
pca$rotation[top_10_genes,1]
genes_1 <- original_gene_counts[as.numeric(top_10_genes),1]
main_genes_PC1 <- data.frame("Gene" = genes_1, "Rot" = pca$rotation[top_10_genes,1])

# Genes associated to higer variance in PC2
loading_scores <- pca$rotation[,2]
gene_scores <- abs(loading_scores)
gene_scores_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_scores_ranked[1:10])
pca$rotation[top_10_genes,2]
genes_2 <- original_gene_counts[as.numeric(top_10_genes),1]
main_genes_PC2 <- data.frame("Gene" = genes_2, "Rot" = pca$rotation[top_10_genes,2])

main_genes <- data.frame("Genes PC1" = genes_1, "Rot" = pca$rotation[top_10_genes,1], "Genes PC2" = genes_2, "Rot" = pca$rotation[top_10_genes,2])
#####################################################################################



#####################################################################################
#                             Clinical Annotation
#####################################################################################

# loads dataset
clinical_annotation <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE)


# Group by tumor category
names_vec <- data.frame(names)
Tumor_binary <- apply(names_vec,1,check_tumor_state)
names(Tumor_binary)[1] <- "Tumor"
ggplot(pca_components) + geom_point(aes(x=pca_components[,1], y=pca_components[,2], color=Tumor_binary)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  ggtitle("Tumor")

# Group by read coverage
total_counts <- data.frame(counts_per_sample)
total_counts$label <- apply(data.frame(total_counts$counts_per_sample),1,check_read_coverage)
ggplot(pca_components) + geom_point(aes(x=pca_components[,1], y=pca_components[,2], color=total_counts$label)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  ggtitle("Read coverage")

check_read_coverage <- function(value){ 
  if (value <= 60000000 ){
    label = "black"
  } else if (value > 60000000 && value <= 90000000){ 
    label = "blue"
  } else if (value > 60000000 && value <= 90000000){ 
    label = "red"
  }else if (value > 90000000 && value <= 120000000){ 
    label = "yellow"
  } else {
    label = "green"
  }
  return(label) 
}

# Create tumor column function 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

check_tumor_state <- function(names){ 
  get_type_str = substrRight(names,2)
  get_type_int = as.numeric(get_type_str)
  
  if (get_type_int <= 10 ){
    state = 1
  }else{ 
    state = 0
  }
  return(state) 
}
  


#####################################################################################
#                             Dispersions
#####################################################################################

# pca non-normalized
pca_without_norm <- prcomp(t(original_gene_counts[-1]))
pca_without_norm_components <- data.frame(pca_without_norm$x)
ggplot(pca_without_norm_components) + geom_point(aes(x=pca_without_norm_components[,1], y=pca_without_norm_components[,2])) +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("Dispersion through PC1 and PC2 - not normalized")


for (col in clinical_names){
  Categories <- clinical_annotation[match(id_samples$Patient.id, clinical_annotation$Patient.ID), col]
  plot <- ggplot(pca_without_norm_components) + geom_point(aes(x=pca_without_norm_components[,1], y=pca_without_norm_components[,2], color=as.factor(Categories))) +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) + 
    xlim(NA, 3000000) +
    ylim(-1e+06, 1e+06)
  col <- gsub("\\.", ' ', col)
  plot <- plot + ggtitle(paste(col, "- Not normalized"))
  print(plot)
}




# pca normalized
pca_components <- data.frame(pca$x)
ggplot(pca_components) + geom_point(aes(x=pca_components[,1], y=pca_components[,2], color=a)) +
   xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
   ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
   ggtitle("Dispersion through PC1 and PC2 - normalized")

# Group by column in clinical annotation
names <- names(original_gene_counts[-1])
patient_ids <- gsub("\\.", '-', substr(names,1,nchar(names)-3))
id_samples <- data.frame("Patient.id" = patient_ids, "Sample" = names(original_gene_counts[-1]))
clinical_names <- names(clinical_annotation)

for (col in clinical_names){
  Categories <- clinical_annotation[match(id_samples$Patient.id, clinical_annotation$Patient.ID), col]
  plot <- ggplot(pca_components) + geom_point(aes(x=pca_components[,1], y=pca_components[,2], color=as.factor(Categories))) +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep = ""))
  col <- gsub("\\.", ' ', col)
  plot <- plot + ggtitle(col)
  print(plot)
}


######## Find the max outlier
n <- 878
Samples <- rep("black", n)
Samples[870] <- "red"
p <- ggplot(pca_components) + geom_point(aes(x=pca_components[,1], y=pca_components[,2], color=Samples)) +
  xlab("PC1") +
  ylab("PC2") +
  labs(color = "Samples") +
  ggtitle("Dispersion through PC1 and PC2 - normalized with max outlier as red point")

p <- p + scale_color_manual(values=c("#000000", "#FF0000")) 
p + scale_fill_discrete(name = "Samples", labels = c("Other samples", "Max outlier"))
print(p)


