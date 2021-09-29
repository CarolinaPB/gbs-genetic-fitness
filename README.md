# Perform GBS and assess genetic fitness

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT
This pipeline starts by scaffolding your genome using long reads. Then it performs GBS using FastGBS2 and outputs stats and a PCA.  
In order to assess genetic fitness, it calculates a distance matrix, IBD and looks at runs of homozygosity

#### Tools used:
- Seqtk - convert fasta to one line fasta
- Longstitch (ntLink-arks) - scaffolding with long reads
- bwa and samtools - index genome
- python - get assembly statistics
- FastGBS2 - perform GBS
- Bcftools - get stats from vcf
- Plink - calculate distance matrix, IBD, homozygosity

| ![DAG](https://github.com/CarolinaPB/population-mapping/blob/wur/workflow.png) |
|:--:|
|*Pipeline workflow* |

## Set up
#### 1. Copy your long reads fq.gz file to the working directory (where the config.yaml and Snakefile are).
#### 2. Follow the **Preparation of the sequence data** from the [FastGBS2 wiki](https://bitbucket.org/jerlar73/fast-gbs_v2/wiki/Home) or see below

Move your sequence files fastq in the data directory. Each file must be in the following format:

Paired-end reads:
```
FLOWCELL_LANE#_R1.fq.gz
FLOWCELL_LANE#_R2.fq.gz
```
Single-end reads:
```
FLOWCELL_LANE#.fq.gz
```
FLOWCELL: Sequencing batch. Originally can have up to 8 lanes


LANES: Each Lane includes sequences of multiplexed DNA samples. Several lanes can be associated to a FLOWCELL.

Suppose you want to analyse 6 sequence files fastq from 3 different FLOWCELLS (HADES, HOMER and ULYSSE) and each FLOWCELL has 2 LANES like:
```
HADES: LANES 1 & 2
HOMER: LANES 3 & 4
ULYSSE: LANES 5 & 6
```
You will have to name your files in this way:
```
HADES_1.fq.gz  
HADES_2.fq.gz  
HOMER_3.fq.gz  
HOMER_4.fq.gz  
ULYSSE_5.fq.gz  
ULYSSE_6.fq.gz
```
In this example the reads are single-end. For PE reads you will have the same naming format plus R1 and R2 after LANE number.

**It is important to keep the extension .fq.gz because, for simplicity, it is hardcoded in the pipeline.**

#### 3. Follow the **Preparation of the barcode files** from the [FastGBS2 wiki](https://bitbucket.org/jerlar73/fast-gbs_v2/wiki/Home) or see below

You can use any spreadsheet software (openOffice, Excel) to format your barcodes files. If you do so, save the active sheet in text (.txt) format with tab separated columns.

The barcode file should be like:
```
CTATTA      BC01
TAACGA      BC02
CGTGTGGTB   BC03
```
Format all your barcode files in the same manner, in a tab-separated column. Move these files in the barcodes directory.

Continuing with the same example as above, you have to create one barcode file per each sequence file and name them as fastq files. See example below:
```
barcodes_HADES_1
barcodes_HADES_2
barcodes_HOMER_3
barcodes_HOMER_4
barcodes_ULYSSE_5
barcodes_ULYSSE_6
```
Do not put any extension, like .txt, to the filenames.

Warning: Do not use special characters (such as: * / \ | $ " etc.) in the names given to your samples in the barcode file.

### Edit config.yaml with the paths to your files
```
ASSEMBLY: /path/to/assembly
LONGREADS: longreads.fq 
PREFIX: <prefix>
GENOME_SIZE: <approximate genome size>
PARAMETERS_FILE: <parameters file>
```

- ASSEMBLY - path to the assembly file
- LONGREADS - name of file with long reads. This file should be in the working directory (where this config and the Snakefile are)
- PREFIX -  prefix for the created files
- GENOME_SIZE - approximate genome size ```haploid genome size (bp)(e.g. '3e9' for human genome)``` from [longstitch](https://github.com/bcgsc/longstitch#full-help-page)
- PARAMETERS_FILE - name of FastGBS2 parameters file
