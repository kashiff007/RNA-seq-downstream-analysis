# RNA-seq differential expression work flow

## Introduction

Here we analyzed RNA-seq reads for estimation of gene expression and differential gene expression. Transcriptomes were extracted from 2020 and 2050 environment condition from Berlin and Rabat.

This document presents an RNAseq differential expression workflow. We will start from the FASTQ files, align to the reference genome, prepare gene expression values as a count table by counting the sequenced fragments, perform differential gene expression analysis, and visually explore the results.


## Input data
We will use publicly available data from the article by Pranav Sahoo et al., 2020; Pecinka Lab. The purpose of the experiment was to investigate the role of the temperature and precipitation between 2020 and 2050 in *Arabidopsis thaliana*. The investigators derived *A. thaliana* somatic cells from 4 different conditions; i.e. Rabat-2020, Rabat-2050, Berlin-2020 and Berlin 2050. Total RNA was isolated from the plant harvested (stage-specific) using RNeasy Mini Kit (QIAGEN) according to manufacturer’s instructions. RNA quality was monitored on Bioanalyzer (Agilent) and the samples with RNA integrity number ≥8 were used for library preparation. A total of eight Illumina type sequencing libraries (four climate scenarios × two biological replicates each) were prepared, multiplexed and sequenced as 150 bp single-end reads on a HiSeq2500 sequencer (Illumina). At least 15 Million reads were sequenced per sample.


## Filtering low quality reads
What we get from the sequencing machine is a set of FASTQ files that contain the nucleotide sequence of each read and a quality score at each position. Before aligning, these reads must undergo quality checkup. [FastQC tool] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) were used to monitor the quality of each read from fastq files. Reads were trimmed/removed having quality less than 25 (25 is the Phread score). For trimming we used [FASTX-Toolkit] (http://hannonlab.cshl.edu/fastx_toolkit/). 

Example:
If the average per base read quality is bad at the 3' end.

![untitled presentation](https://www.researchgate.net/profile/Richard_Tennant/publication/312355161/figure/fig2/AS:450870568591361@1484507328098/Sequence-Quality-Per-base-Before-and-After-Trimming-and-Adapter-Removal-The-per-base.png)

we can trim 5 bases from 3' end by following command:
```
fastx_trimmer -f [read-size - 5] -i Input.fastq -o Input_filtered.fastq
```

Example:
If the average per base read quality is bad throughout the read.

![untitled presentation](https://dwheelerau.files.wordpress.com/2013/03/fastqc_stats.png)

we can filter reads for more the 25 Phread score by following command:
```
fastq_quality_filter -q 25 -p 75 -i Input.fastq -o Input_filtered.fastq
```

Note: Change the command parameters according to desired output.


## Aligning reads to a reference
After quality control reads were align to a reference genome. Here reference genome was [TAIR10](https://www.arabidopsis.org/download/index-auto.jsp%3Fdir%3D%252Fdownload_files%252FGenes%252FTAIR10_genome_release). It is important to know if the sequencing experiment was single-end or paired-end, as the alignment software will require the user to specify both FASTQ files for a paired-end experiment. In our case we had single-end fastq files. The output of this alignment step is commonly stored in a file format called BAM.

Here we used the TopHat2 alignment software for the alignment of reads on TAIR10 genome with the reference of TAIR10 gff file. GFF/GTF files are annotation files which directs the mapping to the gene location and further helps in estimation of fragment numbers for each gene.

For example, the paired-end RNA-Seq reads for the Rabat-2020 (each replicate) aligned using TopHat2 with 20 threads:
```
tophat2 -o Rabat-2020_rep1_tophat_out -p 8 -G path/to/gff/ path/to/genome Rabat-2020_rep1_forward.fastq Rabat-2020_rep1_reverse.fastq
```
This command will generate an accepted.bam file inside Rabat-2020_rep1_tophat_out folder. It also have a file align.summary which has the information about the percentage of reads mapped on the reference genome. 

Bam file is binary file (sam is human readable), which has the information about each read mapped-location on the reference genome. This file can be further used to extract the read counts under each gene for the particular condition. 


## Counting reads in genes
To count how many read map to each gene, we need bam file and transcript annotation. We used htseq-count to get the read counts.

Example:
```
htseq-count sam_file gff_file > count_file.txt
```
In the above sample command its showing sam file; it can interconvert from sam to bam using samtools. 

Read-count can also be estimated through samtools by the following command:
```
bedtools multicov -bams file1.BAM file2.BAM file3.BAM -bed <BED/GFF/VCF>
```
With bedtools one can use multiple bam files at once and a matrix is generated with corresponding multiple columns.


The final count matrix looks like following (example):

| Scaffold  | start  | stop   | geneID    |strand| Rabat-2020_rep1 | Rabat-2020_rep2 | Rabat-2050_rep1 | Rabat-2050_rep2 |
| --------- |:------:|:------:|:---------:|:----:|:---------------:|:---------------:|:---------------:|:---------------:|
| Chr1      | 23573  | 25552  | AT1G01060 |   +  |       234       |       343       |       8778      |      7667       |
| Chr1      | 453533 | 458347 | AT1G03243 |   -  |        2        |       11        |       2234      |      2455       |
| Chr1      | 767688 | 768898 | AT1G34510 |   +  |       4343      |      5454       |       45        |      55         |

Here, numeric value below each sample are read-count for corresponding genes.


## Estimation of normalised gene expression and differential gene expression.

### Reading the matrix file in R
Above read table were uploaded into R and select the matrix only for read-number (cloumn starting from 6 to last)

```
## Load read counts with replicates
CountTable <- as.data.frame(round(as.matrix(read.csv("replicate.read_counts", sep="\t", header=TRUE, row.names=4))))
CountTable <- CountTable[:c(6:)]
```
In the table two replicates of each CA20, CA50, CE20, CE50, TA20, TA50, TE20 and TA50 are present.

### Assigning the columns for respective samples replicate
Then we assign each column to each sample:

```
##Comparisons combination
CA20_CA50 <- (CountTable[,c(1,2,3,4)])
CA20_CE20 <- (CountTable[,c(1,2,5,6)])
CA20_CE50 <- (CountTable[,c(1,2,7,8)])
CA20_TA20 <- (CountTable[,c(1,2,9,10)])
CA20_TA50 <- (CountTable[,c(1,2,11,12)])
CA20_TE20 <- (CountTable[,c(1,2,13,14)])
CA20_TE50 <- (CountTable[,c(1,2,15,16)])
```

### Designing the experiment

We design each experiment by assigning the replicates condition to untreated and treated.

```
CA20_CA50_Design = data.frame(row.names = colnames(CA20_CA50 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_CE20_Design = data.frame(row.names = colnames(CA20_CE20 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_CE50_Design = data.frame(row.names = colnames(CA20_CE50 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_TA20_Design = data.frame(row.names = colnames(CA20_TA20 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_TA50_Design = data.frame(row.names = colnames(CA20_TA50 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_TE20_Design = data.frame(row.names = colnames(CA20_TE20 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
CA20_TE50_Design = data.frame(row.names = colnames(CA20_TE50 ),condition = c( "untreated", "untreated", "treated","treated"), libType = c( "single-end", "single-end", "single-end", "single-end"))
```


### Reading each experiment reads by DEseq2

Reading each experiment for camparion in the form of dataframe through DEseq2 package:

```
library(DESeq2)

CA20_CA50_dds <- DESeqDataSetFromMatrix(countData = CA20_CA50,colData=CA20_CA50_Design, design= ~ condition)
CA20_CE20_dds <- DESeqDataSetFromMatrix(countData = CA20_CE20,colData=CA20_CE20_Design, design= ~ condition)
CA20_CE50_dds <- DESeqDataSetFromMatrix(countData = CA20_CE50,colData=CA20_CE50_Design, design= ~ condition)
CA20_TA20_dds <- DESeqDataSetFromMatrix(countData = CA20_TA20,colData=CA20_TA20_Design, design= ~ condition)
CA20_TA50_dds <- DESeqDataSetFromMatrix(countData = CA20_TA50,colData=CA20_TA50_Design, design= ~ condition)
CA20_TE20_dds <- DESeqDataSetFromMatrix(countData = CA20_TE20,colData=CA20_TE20_Design, design= ~ condition)
CA20_TE50_dds <- DESeqDataSetFromMatrix(countData = CA20_TE50,colData=CA20_TE50_Design, design= ~ condition)

## setting refence Level
dds$condition <- relevel(dds$condition, ref="untreated")
```


### Saving normalized results before differential gene analysis.

DEseq2 uses [median of ratios](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) normalization method which means counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene. If you want get this value for each gene then use the following command: 
```
CA20_CA50_dds_factor <- estimateSizeFactors(dds)
CA20_CA50_dds_factor_normalize <- counts(CA20_CA50_dds_factor, normalized=TRUE)
*similary for other condition comparison*
```


### Run DESeq2 for the comparison
DRSeq command will compare the difference in gene expression for each gene, and for such comparison it calculates the p-value. It further adjust the p-value and hence also estimate the adjusted p-value (or q-value). It also report baseMean, log2FoldChange, lfcSE (lfc - standard error) and stat.

```
CA20_CA50_res <- results(DESeq(CA20_CA50_dds))
CA20_CE20_res <- results(DESeq(CA20_CE20_dds))
CA20_CE50_res <- results(DESeq(CA20_CE50_dds))
CA20_TA20_res <- results(DESeq(CA20_TA20_dds))
CA20_TA50_res <- results(DESeq(CA20_TA50_dds))
CA20_TE20_res <- results(DESeq(CA20_TE20_dds))
CA20_TE50_res <- results(DESeq(CA20_TE50_dds))
```

Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage [(Zhu, Ibrahim, and Love 2018)](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bty895), which improves on the previous estimator.

We provide the dds object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as it appears in resultsNames(dds).


```
CA20_CA50_res_LFC <- lfcShrink(DESeq(CA20_CA50_dds),coef=2, type="apeglm")
CA20_CE20_res_LFC <- lfcShrink(DESeq(CA20_CE20_dds),coef=2, type="apeglm")
CA20_CE50_res_LFC <- lfcShrink(DESeq(CA20_CE50_dds),coef=2, type="apeglm")
CA20_TA20_res_LFC <- lfcShrink(DESeq(CA20_TA20_dds),coef=2, type="apeglm")
CA20_TA50_res_LFC <- lfcShrink(DESeq(CA20_TA50_dds),coef=2, type="apeglm")
CA20_TE20_res_LFC <- lfcShrink(DESeq(CA20_TE20_dds),coef=2, type="apeglm")
CA20_TE50_res_LFC <- lfcShrink(DESeq(CA20_TE50_dds),coef=2, type="apeglm")
```


### For exporting the final processed file into text file:

```
write.csv(as.data.frame(CA20_CA50_res_LFC), 
          file="CA20_CA50_results.csv")

*similary for other condition comparison*      
```         

For much more detail explanation about DESeq2 package please refer to [DESeq2 tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#)


# edgeR workflow:
```
library(edgeR)
known <- read.table("known_miRNA.txt", sep="\t",header=T)
new <- read.table("new_miRNA.txt", sep="\t",header=T)

count_known_H <- known[,c(5,6,8,9)]
count_known_P <- known[,c(11,12,14,15)]
count_new_H <- new[,c(5,6,8,9)]
count_new_P <- new[,c(11,12,14,15)]

group <- c(1,1,2,2)
y_known_H <- DGEList(counts=count_known_H, group=group)
y_known_P <- DGEList(counts=count_known_P, group=group)
y_new_H <- DGEList(counts=count_new_H, group=group)
y_new_P <- DGEList(counts=count_new_P, group=group)

y_known_H <- estimateDisp(y_known_H)
y_known_H <- estimateTagwiseDisp(y_known_H)
y_known_P <- estimateDisp(y_known_P)
y_known_P <- estimateTagwiseDisp(y_known_P)
y_new_H <- estimateDisp(y_new_H)
y_new_H <- estimateTagwiseDisp(y_new_H)
y_new_P <- estimateDisp(y_new_P)
y_new_P <- estimateTagwiseDisp(y_new_P)

et_y_known_H <- exactTest(y_known_H)
et_y_known_P <- exactTest(y_known_P)
et_y_new_H <- exactTest(y_new_H)
et_y_new_P <- exactTest(y_new_P)

```
