library(limma)
library(edgeR)
library(ggplot2)
library(tibble)
library(ggrepel)

clinical_annotation<- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_ClinicalAnnotation.txt', sep = '\t', header = TRUE)
original_gene_counts <- read.csv(file = '~/Downloads/Lab-6-Omics/TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)
rownames(original_gene_counts) <- original_gene_counts$Gene 
original_gene_counts$Gene <- NULL

dge <- DGEList(counts=data.matrix(original_gene_counts))
dge <- calcNormFactors(dge)
v <- voom(dge, plot=TRUE)
normalized_counts = data.frame(v$E)

names_vec = names(original_gene_counts)
patient_ids <- gsub("\\.", '-', substr(names_vec,1,nchar(names_vec)-3))
id_samples <- data.frame("Patient.id" = patient_ids, "Sample" = names(original_gene_counts))

age <- clinical_annotation[match(id_samples$Patient.id, clinical_annotation$Patient.ID), "Age.at.diagnosis..years."]
  
# Create tumor column function 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
check_tumor_state <- function(names){ 
  get_type_str = substrRight(names,2)
  get_type_int = as.numeric(get_type_str)
  
  if (get_type_int <2 ){
    state = 1
  }else if (get_type_int >10 ){ 
    state = 0
  } else{
    state = NA
  }
  return(state) 
}

age <- data.frame(age)
names_vec <- data.frame(names_vec)
Tumor_binary <- data.frame(apply(names_vec,1,check_tumor_state))

# Delete Metastic Rows
indexes_NA = which(is.na(Tumor_binary$apply.names_vec..1..check_tumor_state.))
Tumor_binary <- na.omit(Tumor_binary)
normalized_counts <- normalized_counts[,-c(indexes_NA)]
names_vec <- data.frame(names_vec[-c(indexes_NA),])
age <- data.frame(age[-c(indexes_NA),])

diff_expr_data <- data.frame("Sample" = names_vec, "Age"= age, "Tumor" = Tumor_binary)
colnames(diff_expr_data) <- c("Sample", "Age", "Tumor")

design.matrix <- model.matrix(~Tumor*Age, data=diff_expr_data)
linearfit = lmFit(normalized_counts, design.matrix)
EBFit = eBayes(linearfit)

# Primary Tumor Vs Non-primary Tumor
mypvalue = 0.05
myFC = 3
NormalvsTumour <- topTable(EBFit,coef=2,sort.by="t",number=Inf)
NormalvsTumour_FC <- topTable(EBFit,coef=2,sort.by="t",number=Inf, p.value = mypvalue, lfc = myFC)

significance_ht <- NormalvsTumour$adj.P.Val
i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < mypvalue){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, NormalvsTumour$logFC, NormalvsTumour$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(NormalvsTumour), plot_ht) 
NvsT_FC <- cbind(gene=rownames(NormalvsTumour_FC), NormalvsTumour_FC) 
interm_FC <- (match(plot_ht$gene,NvsT_FC$gene, nomatch = NA))
gene_labels <- data.frame(NvsT_FC$gene[interm_FC])
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Normal vs Tumor(Primary)") +
  geom_text_repel(data=head(plot_ht, 11), aes(label = NvsT_FC.gene.interm_FC.)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 









# Primary Tumor Vs Non-primary Tumor * Age
mypvalue = 0.05
myFC = 0.05
NormalvsTumourvsAge <-topTable(EBFit,coef=4,sort.by="t",number=Inf)
NormalvsTumourvsAge_FC <- topTable(EBFit,coef=4,sort.by="t",number=Inf, p.value = mypvalue, lfc = myFC)

significance_ht <- NormalvsTumourvsAge$adj.P.Val
i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < mypvalue){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, NormalvsTumourvsAge$logFC, NormalvsTumourvsAge$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(NormalvsTumourvsAge), plot_ht) 
NvsTvsA_FC <- cbind(gene=rownames(NormalvsTumourvsAge_FC), NormalvsTumourvsAge_FC) 
interm_FC <- (match(plot_ht$gene,NvsTvsA_FC$gene, nomatch = NA))
gene_labels <- data.frame(NvsTvsA_FC$gene[interm_FC])
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color = guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Normal vs Tumor(Primary) * Age") +
  geom_text_repel(data = head(plot_ht, 38), aes(label = NvsTvsA_FC.gene.interm_FC.)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

# Produce file for GSEA NormalvsTumour
GoAnalysis <- data.frame (NormalvsTumourvsAge$t)
colnames(GoAnalysis)<- c("t")
GoAnalysis <- cbind(gene = rownames(NormalvsTumourvsAge), GoAnalysis)
GoAnalysis_sorted <- GoAnalysis[order(-GoAnalysis$t) ,]
write.table(GoAnalysis_sorted, '~/Downloads/Lab-6-Omics/GSEAtable_AGE.txt', sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE) 



