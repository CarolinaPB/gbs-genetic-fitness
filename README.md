# Perform GBS and assess genetic fitness

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT
This pipeline starts by scaffolding your genome using long reads. Then it performs GBS using FastGBS2 (it skips the imputation step) and outputs stats and a PCA.  
In order to assess genetic fitness, it calculates a distance matrix, IBD and inbreeding.

#### Tools used:
- Seqtk - convert fasta to one line fasta
- Longstitch (ntLink-arks) - scaffolding with long reads
- bwa and samtools - index genome
- python - get assembly statistics, plot interactive PCA
- FastGBS2 - perform GBS
- Bcftools - get stats from vcf
- Plink - calculate distance matrix, IBD, PCA, inbreeding

#### If you want to skip the scaffolding step do this:
1. Add your genome to the `refgenome` directory and name it as <prefix>.fa (it's important to have the .fa extension). 
2. Index with bwa by running the `index_genome.sh` script in the `scripts` directory from the main directory as `./scripts/index_genome.sh`
You will get two files with assembly stats, but both will have the same information.

# REPLACE WORKFLOW
| ![DAG](https://github.com/CarolinaPB/gbs-genetic-fitness/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |

## Set up
#### 1. Copy your long reads fq.gz file to the working directory (where the config.yaml and Snakefile are).
The reads must have .fq.gz or .fa.gz extension
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

#### 4. Follow the **Preparation of the parameters file** from the [FastGBS2 wiki](https://bitbucket.org/jerlar73/fast-gbs_v2/wiki/Home) AND see specific instructions below
**You'll also need to do this in addition to what is written in the FastGBS wiki**  
Copy the contents of parameters_template.txt to a new file and name that file **parameters_<prefix>.txt**, where `<prefix>` is the name you'll use on the config file and that will be appended to the output. For example, it could be the name of your organism `parameters_chicken.txt`.   
In the text you copied from `parameters_template.txt` there are some variables named `<prefix>_othername`. Where you see `<prefix>` you'll need to change it to the prefix you've chosen.   
 From [FastGBS2 wiki](https://bitbucket.org/jerlar73/fast-gbs_v2/wiki/Home) -- If a step has **DIFFERENT**, follow the instructions in this page instead of the ones in the FastGBS wiki. These are changes specific to this pipeline  
  
**change only the fields in this section**  
Parameter file is text file including all information related to GBS experiment and required data to run Fast-GBS_V2

Below is an example of parameter file included in Fast-GBS pipeline. Here we separated each section with a description:

1- Naming log file, FLOWCELL and LANES. This part is been detailed in the Preparation of the sequence data section  
  **DIFFERENT -- edit the log file name with your prefix \<prefix>_logfile_fastgbs.log**
```
; NAME OF FILE CONTAINING LOGS
LOGFILE=<PREFIX>_logfile_fastgbs.log

; FLOWCELL INFORMATION 
FLOWCELL=ULYSSE HOMER HADES 

; LANES INFORMATION
ULYSSE_LANES=5 6
HOMER_LANES=3 4
HADES_LANES=1 2
```
2- The parameters used in Fast-GBS pipeline is different for each sequencing technology and sequence types (single- or paired-end). Please select the right technology and sequence type.
```
; SEQUENCING TECHNOLOGY: ILLUMINA IONTORRENT
TECHNOLOGY=ILLUMINA

; SEQUENCE TYPE: SE (Single End) or PE (Paired Ends)
SEQTYPE=SE
```
3- Define the reference genome. The name should be exactly the same. If you want to use SRG please visit here. This part is been detailed in the Preparation of the reference genome section.  
**DIFFERENT -- edit the REFGEN file name with your prefix**
```
; NAME OF THE REFERENCE GENOME FILES
REFGEN=<prefix>.fa
```
4- You should get the adaptor recognition sequence from your sequencing service platform. This sequence here is for Illumina and used by IBIS
```
; SEQUENCE USED FOR ADAPTOR RECOGNITION (THE COMPLETE ADAPTOR SEQUENCE IS
; AGATCGGAAGAGCGGG)
ADAPFOR=AGATCGGAA

; IF PAIRED END READS, ADD ADAPTER SEQUENCE IN REVERSE SEQUENCE (R2)
ADAPREV=AGATCGGAA
```
5- Sequencing reads length are different based on different technology, also dimers are short reads. Here you define the minimum length of the reads for mapping and variant calling.
```
; MINIMUM READ LENGTH TO KEEP
READLEN=50
```


10- If you want to call variant from a database write the name VCF file desired to be used. Either write None.   
**DIFFERENT - NOT TESTED IN THIS PIPELINE**
```
; VARIANT CALING FROM A DATABASE: None or the name of a vcf file
SOURCE=file_option.vcf
```
11- You can here define the minimum depth of coverage for variant calling.
  ```
; FOR PLATYPUS: Minimum number of supporting reads required before a
; variant candidate will be considered (equivalent to minDP in vcftools)
MINREADS=2
  ```
12- **DIFFERENT - edit the file name with your prefix**
```
; PLATYPUS LOGFILE
LOGPLAT=<prefix>_FastGBS_platypus_log.txt
```
13- Name your result VCF file name  
 **DIFFERENT - edit the file name with your prefix**
```
; PLATYPUS OUTPUT VCF FILE
OUTPLAT=<prefix>_FastGBS_platypus
```
14- Filtering steps
```
;PROPORTION OF MISSING DATA (0 to 1 (no missing data)) AND AND MAF (minimum allele frequency, 0 to 1)
MAX-MISSING=0.2
MAF=0.01
```
Verify the information in the file parameters.txt. If necessary, change the value given to the variables. It is very important not to change the words in capital letters (before the sign =). If you do this, you will get an error message because the pipeline will not find that variable. Respect what is an integer and a string.
 


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

 ## Working directory
 There are many files and directories that you don't need to pay attention to, but they need to be there for the pipeline to run.  
The directories and files that should not be touched beyond what is stated in the preparation steps are:
- **barcodes** - where you'll add your barcode files
- **data** - where you'll add your fq.gz files and where the demultiplexed fq files will be added, as well as the mapped .bam files
- **refgenome** where the scaffolded genome will be
- **reject**
- **rules** extra snakemake rules
- **scripts** scripts used in the pipeline
- **stat**  where the missing samples will be. Samples with less than 10% of the mean number of reads in the pool will be moved here and removed from next steps
- **wiki** wiki from FastGBS2
- FastGBS scripts:
  - count_nbseq_V2.sh
  - fastgbs_V2.sh
  - makeBarcodeSabre_V2.py
  - makeDir_V2.sh
  - SLURM_GBS.sh
  - Summary4VCF.py
  - txt2unix.sh
- Many files with the name starting in \<prefix>_oneline.fa will be created by the scaffolding tool. The scaffolded genome will be in the **refgenome** directory.  
 This pipeline saves the results to several folders in the running directory. Since this pipeline will create many files, some of them saved to the main directory, it will look quite messy.  
 Below you'll find an overview of the most important results
 
 
 ## RESULTS
The working directory will be messy with all the necessary files and results from the several pipeline steps.
The most important files are and directories are:  
- **<run_date>_files.txt** dated file with an overview of the files used to run the pipeline (for documentation purposes)
- **refgenome** directory that contains scaffolded genome
- **results** directory that contains
  - **assembly_stats_<prefix>_new.txt** and **assembly_stats_<prefix>_old.txt** files with assembly statistics for the new scaffolded assembly and for the assembly before the extra scaffolding
  - **{prefix}_FastGBS_platypus.vcf** raw VCF file including all variants with PASS or FILTERED flag
  - **{prefix}_FastGBS_platypus.recode.vcf** filtered VCF file including only variants with PASS flag and meeting filtering criteria defined in parameters file
  - **{prefix}_FastGBS_platypus.GT.FORMAT** VCF file with only genotypes
  - **Summary_By_Samples_python.txt** and **Summary_By_Sites_python.txt** GBS summary stats from FastGBS2
  - **{prefix}_pca.html** interactive PCA where you can zoom in and hover a certain point to get sample name (open with browser)
  - **{prefix}.eigenvec** and **{prefix}.eigenval** output of PCA
 - **fitness** directory that contains results from fitness analysis
   - **{prefix}.het.gz** inbreeding results: [Plink](https://www.cog-genomics.org/plink/1.9/basic_stats) computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates 
   - **{prefix}.genome.gz** IBD results (see [here](https://www.cog-genomics.org/plink/1.9/ibd))
   - **{prefix}.dist.gz** and **{prefix}.dist.id** distance matrix and sample ids as described [here](https://www.cog-genomics.org/plink/1.9/distance)

## OTHER
 If you need to install platypus yourself:
 1. Download Platypus from [here](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data)
 2. Load python 2.7 and htslib. If you're in Anunna:
 ```
 module load python/2.7.15
 module load htslib
 ```
 3. Build
 ```
tar -xvzf Platypus_x.x.x.tgz
cd Platypus_x.x.x
./buildPlatypus.sh
 ```
 4. Add Platypus to your PATH
 ```
 export PATH=/path/to/intallation/Platypus_x.x.x:$PATH
 ```
 
 If you use your own installation of Platypus, in the Snakefile, rule `run_fast_GBS2`, you should edit the line `export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/Platypus_0.8.1:$PATH` to the path to your installation, as above.
