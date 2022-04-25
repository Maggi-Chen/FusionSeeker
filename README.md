# FusionSeeker

A gene fusion caller for long-read single-molecular sequencing data.

Author: Maggi Chen

Email: maggic@uab.edu

Draft date: Apr. 23, 2022

## Quick Start
```sh
git clone https://github.com/Maggi-Chen/FusionSeeker.git
cd FusionSeeker/
./fusionseeker -h

# quick gene fusion calling with sorted bam
fusionseeker --bam merged.sort.bam --datatype isoseq --outpath fusionseeker_out/ --human38

# Gene fusion discovery with custom reference
fusionseeker --bam merged.sort.bam --datatype nanopore --outpath fusionseeker_out/ --ref hg38.fa --gtf Human_hg38.gtf 

```



## Description

FusionSeeker is a tool for gene fusion discovery with long-read transcriptome sequencing data. The input should be a sorted BAM (PacBio Iso-Seq, Nanopore, or mixed platform). The output is a list of confident gene fusions and their transcript sequences. By default, FusionSeeker uses Human GRCh38 and Ensembl annotation v104 (Homo_sapiens.GRCh38.104.chr.gtf.gz) as reference.<br />
When using custom reference, make sure the chromosome name in BAM, GTF, and reference_genome.fa are identical. By default, FusionSeeker only considers gene with valid "gene_name" in GTF and skips the remaining genes, unless --geneid is set.<br />
This program was tested on a x86_64 Linux system with a 128GB physical memory.


## Depencency

Dependencies for FusionSeeker:

* python3
* pysam  (tested with version 0.17.0)
* minimap2  (tested with version 2.24)
* samtools  (tested with version 1.9)
* bsalign  (tested with version 1.2.1)


## Installation

```
git clone https://github.com/Maggi-Chen/FusionSeeker.git
```
Then, please also add this directory to your PATH:
```
export PATH=$PWD/FusionSeeker/:$PATH
```


To simplify the environment setup process, Anaconda2 (https://www.anaconda.com/) is recommended:
```
conda create --name fusions
conda activate fusions
conda install -c bioconda minimap2=2.24 pysam=0.17 samtools=1.9 -y
git clone https://github.com/ruanjue/bsalign.git
cd bsalign && make
export PATH=$PWD:$PATH

```

A test dataset is available to verify successful installation:
```
fusionseeker --bam FusionSeeker/testdata/test.bam  -o test_out/ --datatype isoseq --ref FusionSeeker/testdata/test.fa.gz
```
Output should be identical to confident_genefusion.txt and confident_genefusion_transcript_sequence.fa in the testdata folder, with 1 gene fusion and its transcript sequence. 
(The FusionSeeker gene fusion discovery on test dataset should finish within several minutes with 4 CPUs and 2GB memory.)


## General usage


```
fusionseeker [-h] --bam <sort.bam>

Gene fusion caller for long-read sequencing data

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --bam BAM             Input sorted BAM. index required
  --datatype DATATYPE   Input read type (isoseq, nanopore) [nanopore]
  --gtf GTF             Genome annotation file
  --ref REF             Reference genome. Required for breakpoint polishing
  --geneid              Use Gene ID instead of Gene name [False]
  --human38             Use reference genome and GTF for Human GCRh38 (default)
  --human19             Use reference genome and GTF for Human GCRh37
  -o OUTPATH, --outpath OUTPATH
                        Output directory [./fusionseeker_out/]
  -s MINSUPP, --minsupp MINSUPP
                        Minimal reads supporting an event [auto]
  --maxdistance MAXDISTANCE
                        Maximal distance to cluster raw signals [20 for isoseq, 40 for nanopore]
  --keepfile            Keep intermediate files [False]
  --thread THREAD       Number of threads [8]


```

## Use cases
FusionSeeker requires a input of read alignment results in BAM format sorted by coordinates. If you start with sequencing reads (Fasta or Fastq format), you may use minimap2 and samtools to map them to a reference genome before you can apply FusionSeeker:
```
# PacBio Iso-Seq
minimap2 -ax splice:hq reference.fa  isoseq.fastq | samtools sort -o isoseq.bam
samtools index isoseq.bam
# Nanopore
minimap2 -ax splice reference.fa  nanopore.fastq | samtools sort -o nanopore.bam
samtools index nanopore.bam
```

FusionSeeker can be applied with built-in Human reference genome (hg38) and annotation (Ensembl v104):
```
fusionseeker --bam isoseq.bam --datatype isoseq -o fusionseeker_out/
```
Or with custom reference genome and annotation (Make sure the chromosome name in both files are identical to those in BAM file):
```
fusionseeker --bam nanopore.bam  --datatype nanopore  -o fusionseeker_out/ --gtf annotation.gtf --ref reference.fa
```

By default, FusionSeeker uses only gene records with valid "gene_name" in the GTF file. To include all genes in the GTF file, use Gene ID instead:
```
fusionseeker --bam isoseq.bam --datatype isoseq -o fusionseeker_out/ --geneid 
```


### Options of FusionSeeker
#### 1. --minsupp, minimal number of supporting reads
--min_supp is the most important argument for gene fusion candidate filtering of FusionSeeker. It is used to remove false-positive signals generated during sequencing or read alignment processes.
By default, FusionSeeker estimates the volumn of noise signals from input dataset to assign a resonable --minsupp. If you find number of gene fusions is too few under default settings, you can speficy a lower --minsupp cutoff to allow in more candidates:
```
fusionseeker --bam isoseq.bam --datatype isoseq -o test_out/ --minsupp 5 
```
It is not suggested to set a --minsupp below 3, unless the sequencing depth of input dataset is extremely low.

#### 2. --maxdistance, maxminal distance cutoff in density-based spatial clustering of applications with noise
This option adjusts max distance cutoff used for clustering gene fusion raw signals. By default, FusionSeeker sets 20 for highly accurate reads (IsoSeq) and 40 for noisy reads (Nanopore).
you can set a larger value of --maxdistance to tolerate more shifts in the breakpoint positions of raw signals:

```
fusionseeker --bam isoseq.bam --datatype isoseq -o test_out/ --maxdistance 100
```

#### 3. --ref, reference genome
Input reference genome allows FusionSeeker to align transcript sequences and refine breakpoint positions of confident gene fusion calls. Make sure to provide the same reference file used for read alignment.
By default, FusionSeeker does NOT refine breakpoint positions when no reference genome is provided. 
(minimap2(>=2.24) is required to map transcript sequences to the reference.)
```
fusionseeker --bam isoseq.bam --datatype isoseq -o test_out/  --ref reference.fa
```


## Output files
The output directory includes:
```
confident_genefusion.txt                       A list of confident gene fusion calls from input BAM file. Includes gene names, breakpoint positions, number and name of fusion-supporting reads.
confident_genefusion_transcript_sequence.fa    Transcript sequences of reported confident gene fusion calls. 
clustered_candidate.txt                        A full list of gene fusion candidates before applying filters.
rawsignal.txt                                  A list of all gene fusion raw signals.
log.txt                                        Log file for debug.
raw_signal/                                    Intermediate files during raw signal detection. Removed by default.
poa_workspace/                                 Intermediate files during transcript sequence generation with POA. Removed by default.
```



