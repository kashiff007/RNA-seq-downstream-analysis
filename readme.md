# RNA-Seq differential expression work flow

## Introduction

Here we analyzed RNA-seq reads for estimation of gene expression and differential gene expression. Transcriptomes were extracted from 2020 and 2050 environment condition from Berlin and Rabat.

This document presents an RNAseq differential expression workflow. We will start from the FASTQ files, align to the reference genome, prepare gene expression values as a count table by counting the sequenced fragments, perform differential gene expression analysis, and visually explore the results.

## Input data
We will use publicly available data from the article by Pranav Sahoo et al., 2020; Pecinka Lab. The purpose of the experiment was to investigate the role of the temperature and precipitation between 2020 and 2050 in *Arabidopsis thaliana*. The investigators derived *A. thaliana* somatic cells from 4 different conditions; i.e. Rabat-2020, Rabat-2050, Berlin-2020 and Berlin 2050. Total RNA was isolated from the plant harvested (stage-specific) using RNeasy Mini Kit (QIAGEN) according to manufacturer’s instructions. RNA quality was monitored on Bioanalyzer (Agilent) and the samples with RNA integrity number ≥8 were used for library preparation. A total of eight Illumina type sequencing libraries (four climate scenarios × two biological replicates each) were prepared, multiplexed and sequenced as 150 bp single-end reads on a HiSeq2500 sequencer (Illumina). At least 15 Million reads were sequenced per sample.

## Filtering low quality reads
What we get from the sequencing machine is a set of FASTQ files that contain the nucleotide sequence of each read and a quality score at each position. Before aligning, these reads must undergo quality checkup. [FastQC tool] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) were used to monitor the quality of each read from fastq files. Reads were trimmed/removed having quality less than 25 (25 is the Phread score). For trimming we used [FASTX-Toolkit] (http://hannonlab.cshl.edu/fastx_toolkit/). 

Example:
If the average per base read quality is bad at the 3' end.

we can trim 5 bases from 3' end by following command:
```
fastx_trimmer -f 5 -i Input.fastq -o Input_filtered.fastq
```

## Aligning reads to a reference
These reads must first be aligned to a reference genome or transcriptome. It is important to know if the sequencing experiment was single-end or paired-end, as the alignment software will require the user to specify both FASTQ files for a paired-end experiment. The output of this alignment step is commonly stored in a file format called BAM.

Here we use the TopHat2 spliced alignment software in combination with the Bowtie index available at the Illumina iGenomes.

For example, the paired-end RNA-Seq reads for the parathyroidSE package were aligned using TopHat2 with 8 threads, with the call:

 tophat2 -o file_tophat_out -p 8 path/to/genome file_1.fastq file_2.fastq samtools sort -n file_tophat_out/accepted_hits.bam _sorted 

The second line sorts the reads by name rather than by genomic position, which is necessary for counting paired-end reads within Bioconductor. This command uses the SAMtools software.

The BAM files for a number of sequencing runs can then be used to generate count matrices, as described in the following section.

Example BAM files
