#RNA-Seq differential expression work flow

##Introduction

Here we analyzed RNA-seq reads for estimation of gene expression and differential gene expression.
esmation and One of the aim of RNAseq data analysis is the detection of differentially expressed genes. The package DESeq2 provides methods to test for differential expression analysis.

This document presents an RNAseq differential expression workflow. We will start from the FASTQ files, align to the reference genome, prepare gene expression values as a count table by counting the sequenced fragments, perform differential gene expression analysis, and visually explore the results.

##Input data
We will use publicly available data from the article by Felix Haglund et al., J Clin Endocrin Metab 2012. The purpose of the experiment was to investigate the role of the estrogen receptor in parathyroid tumors. The investigators derived primary cultures of parathyroid adenoma cells from 4 patients. These primary cultures were treated with diarylpropionitrile (DPN), an estrogen receptor beta agonist, or with 4-hydroxytamoxifen (OHT). RNA was extracted at 24 hours and 48 hours from cultures under treatment and control.

Part of the data from this experiment is provided in the Bioconductor data package parathyroidSE.

Aligning reads to a reference
What we get from the sequencing machine is a set of FASTQ files that contain the nucleotide sequence of each read and a quality score at each position. These reads must first be aligned to a reference genome or transcriptome. It is important to know if the sequencing experiment was single-end or paired-end, as the alignment software will require the user to specify both FASTQ files for a paired-end experiment. The output of this alignment step is commonly stored in a file format called BAM.

Here we use the TopHat2 spliced alignment software in combination with the Bowtie index available at the Illumina iGenomes.

For example, the paired-end RNA-Seq reads for the parathyroidSE package were aligned using TopHat2 with 8 threads, with the call:

 tophat2 -o file_tophat_out -p 8 path/to/genome file_1.fastq file_2.fastq samtools sort -n file_tophat_out/accepted_hits.bam _sorted 

The second line sorts the reads by name rather than by genomic position, which is necessary for counting paired-end reads within Bioconductor. This command uses the SAMtools software.

The BAM files for a number of sequencing runs can then be used to generate count matrices, as described in the following section.

Example BAM files
