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
After quality control reads were align to a reference genome. Here reference genome was [TAIR10](https://www.arabidopsis.org/download/index-auto.jsp%3Fdir%3D%252Fdownload_files%252FGenes%252FTAIR10_genome_release). It is important to know if the sequencing experiment was single-end or paired-end, as the alignment software will require the user to specify both FASTQ files for a paired-end experiment. In our case we had paired-end fastq files. The output of this alignment step is commonly stored in a file format called BAM.

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

