# Applied Omics

*Project developed in the Bioinformatics course, at IST (Lisboa).*

The objective of this project was to analyze a clinically annotated breast cancer transcriptomics dataset from The Cancer Genome Atlas (TCGA), namely to estimate gene expression in a cohort of tumors and some matched normal samples and use that information to find molecular features potentially relevant for understanding the biology of the disease and for its clinical management. 

Some of the files required as input or output of some of code scripts can not be found in this repository as those could be found online or were too large. Those files were:
+ TCGA-C8-A138-01_1.fastq
+ TCGA-C8-A138-01_2.fastq
+ TCGA_BRCA_ClinicalAnnotation.txt
+ TCGA_BRCA_Gene_ReadCounts.txt
+ GRCh38_latest_rna.fna
+ mart_export.txt
+ transcripts.idx 

The structure of this project was guided by several questions.

Therefore, inside the directory **code**, you can find the code scripts (in ***R***) used to answer each question. 

In the **report.pdf** file, you can find our answers to each of the questions, containing not only the code outputs but also the explanation of *rationale* guiding our decisions and conclusions.

---
## Project Questions

### **Group I**

***a)*** Quality assessment of the raw sequencing data in the FASTQ files provided.

***b)*** Gene expression estimation in the sample using an aligner and an annotated human transcriptome. Compare (and comment) these estimates with those provided in the read count table for the same sample.

### **Group II**

***a)*** Sum up the distributions of read coverage and library complexity across TCGA samples based on the provided read count table.

***b)*** Transform and normalize the data to have gene expression profiles comparable among samples. Justify your procedure, showing that it worked, and analyzing the presence of alarming samples.

***c)*** Determine which phenotypic traits dominate data variance and which genes are associated with the main axes of variance. Analyze the presence of any non-biological batch effect worth acting on.

***d)*** Determine which are the main differences in expressed genes and activated pathways between primary tumors and normal breast samples. Investigate how the age of patients influences those differences.

### **Group III**

***a)*** Breast tumors are typically tested for some proteins:
+ Estrogen, encoded by genes ESR1 and ESR2;
+ Progesterone, encoded by gene PGR;
+ Human Epidermal growth factor Receptor 2 (HER2), encoded by gene ERBB2.

The applied treatment differs accordingly to the presence/absence of these proteins.

The results of such tests for the TCGA primary tumour samples are provided in a patient annotation table.

 Determine how good is the cognate genes' mRNA expression at recovering the binary classifications that result from the aforementioned (protein presence/absence) tests.

***b)*** *"PAM50 is a 50-gene signature that classifies breast cancer into five molecular intrinsic subtypes: Luminal A, Luminal B, HER2-enriched, Basal-like and Normal-like (v. corresponding column in the provided patient annotation table). This classification is also used in the clinic to support therapeutic decisions."*

Rank the five subtypes in terms of prognosis, according to an overall survival analysis.

***c)*** Find the gene signature that best classifies the molecular subtype, explicitly assessing its performance, based on gene expression. How many of the genes in the obtained signature coincide with the ones from PAM50*? Compare the performance of obtained gene signature with the one resulting from PAM50 genes.

\* - *"UBE2T, BIRC5, NUF2, CDC6, CCNB1, TYMS, MYBL2, CEP55, MELK, NDC80, RRM2, UBE2C, CENPF, PTTG1, EXO1, ORC6L, ANLN, CCNE1, CDC20, MKI67, KIF2C, ACTR3B, MYC, EGFR, KRT5, PHGDH, CDH3, MIA, KRT17, FOXC1, SFRP1, KRT14, ESR1, SLC39A6, BAG1, MAPT, PGR, CXXC5, MLPH, BCL2, MDM2, NAT1, FOXA1, BLVRA, MMP11, GPR160, FGFR4, GRB7, TMEM45B, ERBB2"*

***d)*** Many samples do not have assigned a PAM50 subtype. Using their gene expression and the signature from *c)*, classify them into their subtype. Use the information provided to give you an hint (as no ground truth is provided) on the classifier performance.
