# ChromE-seq-Alignment-Pipeline_V1
Pipeline for demultiplexing and processing paired ChromE-seq data generated using a two-round combinatorial indexing strategy.
This pipeline is adapted from the SHARE-seq alignment workflow and includes support for two-round combinatorial barcoding, with 192 barcodes in each round, as well as Hi-C processing and RNA-seq processing.

Original SHARE-seq pipeline author: Sai Ma
ChromE-seq author: Zhicong Liao
Modified for two-round Hi-C/RNA processing

# Overview
This pipeline processes sequencing data from combinatorial indexing experiments and supports multiple assay types, including:
- DipC2 — Hi-C data generated using the two-round barcoding strategy
- RNA2 — RNA data generated using the two-round barcoding strategy
- ATAC
- RNA

For the paired DipC2/RNA2 workflow described here, the typical configuration is:
Project=(Demo_DNA Demo_RNA)
Type=(DipC2 RNA2)
Genomes=(hg38 hg38)

# Installation
The pipeline requires a Linux environment with the required command-line tools installed and available in the system PATH.

Major dependencies include:
- GNU Parallel
- Python 3
- R
- fastp
- pigz
- STAR
- bowtie2
- bwa
- bwa-mem2
- biscuit
- samtools
- bedtools
- featureCounts
- umi_tools
- RSeQC
  - geneBody_coverage.py
  - read_distribution.py
- pairtools
- cooler
- bgzip/tabix
- Picard
- Juicer Tools
- igvtools
- wigToBigWig
- trim_galore
Additional Python and R packages required by the scripts in scripts/SHARE-seq-nova/ must also be installed.

We recommend creating and activating a dedicated Conda environment before running the pipeline.

For example:
conda activate splitseq
The pipeline uses executables from the active Conda environment or system PATH. It does not require modification of the user's .bashrc.

# Repository structure
The recommended repository structure is:
Split-seq/
├── config_demo.yaml
├── demo_fastq/
├── Genome/
│   ├── biscuit/
│   ├── bowtie2/
│   ├── bwa/
│   ├── bwa-mem2/
│   └── star/
├── output/
├── references/
│   ├── genomeBed/
│   ├── gtf/
│   └── TSSfiles/
├── scripts/
│   └── SHARE-seq-nova/
└── Split_seq_V7.7.nova.walk_bothgenome.2rounds2plates.sh
Genome indices can be stored elsewhere on the system. The paths can be specified using environment variables instead of placing large genome files directly in the GitHub repository.

# Reference files
By default, genome indices are expected under:
Genome/
├── biscuit/
├── bowtie2/
├── bwa/
├── bwa-mem2/
└── star/
For example, for an hg38 DipC2 experiment using bwa-mem2, the pipeline expects a reference corresponding to:
Genome/bwa-mem2/hg38.fa
RNA alignment with STAR expects the STAR index directory:
Genome/star/hg38/
Biscuit references used for methylation-aware processing are expected under:
Genome/biscuit/

# Annotation files
Reference annotations are organized under:
references/
├── genomeBed/
├── gtf/
└── TSSfiles/
The genomeBed directory contains chromosome-size and other genome annotation files used by the pipeline.

The gtf directory contains gene annotation files used by featureCounts.

The TSSfiles directory contains transcription start site annotations used for QC.

# Configuring paths
The pipeline is designed so that machine-specific paths do not need to be hard-coded into the script.

The default base directory is:
basedir="${BASEDIR:-/path/to/Split-seq}"
Before running the pipeline, set BASEDIR to the cloned repository:
export BASEDIR=/path/to/Split-seq
The pipeline then uses:
$BASEDIR/demo_fastq
$BASEDIR/output
$BASEDIR/config_demo.yaml
$BASEDIR/scripts/SHARE-seq-nova
$BASEDIR/references
$BASEDIR/Genome
Individual locations can also be overridden:
export RAWDIR=/path/to/fastq
export OUTDIR=/path/to/output
export CONFIG_YAML=/path/to/config.yaml
export PIPELINE_SCRIPTS=/path/to/SHARE-seq-nova
export REFERENCE_DIR=/path/to/references
export GENOME_DIR=/path/to/genome/indexes
This is particularly useful when genome indices are stored on a shared server.

For example:
export BASEDIR=/path/to/Split-seq
export GENOME_DIR=/mnt/Genome

# Input FASTQ files
The pipeline can start from:
Start=Bcl
Start=Fastq
or:
Start=Demultiplexed
For standard Illumina FASTQs containing four files per lane, the expected reads are:
R1 = biological read 1
R2 = biological read 2
I1 = index read 1
I2 = index read 2
Examples of supported naming patterns include:
*_S0_L001_R1_001.fastq.gz
*_S0_L001_R2_001.fastq.gz
*_S0_L001_I1_001.fastq.gz
*_S0_L001_I2_001.fastq.gz
and single-lane files such as:
*_S1_R1_001.fastq.gz
*_S1_R2_001.fastq.gz
*_S1_I1_001.fastq.gz
*_S1_I2_001.fastq.gz

# Sample configuration
The main sample configuration is located near the beginning of:
Split_seq_V7.7.nova.walk_bothgenome.2rounds2plates.sh
For a paired DipC2/RNA2 experiment:
Project=(Demo_DNA Demo_RNA)
Type=(DipC2 RNA2)

Genomes=(hg38 hg38)
ReadsPerBarcode=(100 100)
The entries in these arrays correspond by position.

Therefore:
Demo_DNA → DipC2 → hg38
Demo_RNA → RNA2  → hg38

# General run configuration
Start=Fastq
Runtype=full
chem=rev
cores=32
Start
Input type:
Bcl
Fastq
Demultiplexed
Runtype
QC
full
QC processes a subset of reads for rapid quality assessment.

full processes the complete dataset.
chem
Sequencing chemistry/orientation:
rev
fwd
For example:
chem=rev
is used for NovaSeq 1.5/NextSeq-style reverse-complement index orientation.

# YAML configuration
Sample-specific PCR/index information is specified in the YAML configuration file.

Example:
Project1:
  Name: Demo_DNA
  Type: DipC2
  P5: P1.29
  P7: P2.29

Project2:
  Name: Demo_RNA
  Type: RNA2
  P5: P1.30
  P7: P2.30
  N6: T
The project names and assay types in the YAML file should match those specified in the main Bash script:
Project=(Demo_DNA Demo_RNA)
Type=(DipC2 RNA2)
For RNA libraries using random priming, the YAML configuration can include:
N6: T

# DipC2 options
DipC2 data are processed through the DipC branch of the pipeline.

A typical configuration is:
Type=(DipC2 RNA2)
taps=(F F)
For standard DipC2 without TAPS/methylation processing:
taps=F
DipC2 reads are adapter-trimmed and aligned using bwa-mem2.

The expected reference is:
Genome/bwa-mem2/<genome>.fa
For example:
Genome/bwa-mem2/hg38.fa
The DipC workflow subsequently uses pairtools to generate and process chromatin contact pairs.

Typical contact-map outputs include:
*.pairs.gz
*.pairs.unique.gz
*.pairs.unique.final.gz
*.contacts.final.gz
*.cool
*.mcool
*.hic

# RNA2 options
RNA2 data are processed through the RNA branch of the pipeline.

RNA-specific settings include:
removeSingelReadUMI=F
keepIntron=T
cores=32
genename=gene_name
refgene=gencode
mode=fast
keepMultiMapping=(F F)
N6=(T T)
removeSingelReadUMI
F = retain single-read UMIs
T = remove UMIs represented by only one read
Default:
removeSingelReadUMI=F
keepIntron
T = include intronic reads
F = use exon-associated reads only
Default:
keepIntron=T
genename
Gene identifier used for feature assignment:
gene_name
gene_id
mode
Two RNA UMI-processing modes are supported:
fast
regular
fast uses the custom deduplication workflow.

regular uses umi_tools.
N6
Controls processing of random-primer-derived RNA reads.

For example:
N6=(T T)
The corresponding RNA project can additionally specify:
N6: T

# Running the pipeline
Clone or download the repository and enter the repository directory.

Set the repository location:
export BASEDIR=/path/to/Split-seq
If the genome indices are stored elsewhere:
export GENOME_DIR=/path/to/Genome
Activate the software environment:
conda activate splitseq
Then run:
./Split_seq_V7.7.nova.walk_bothgenome.2rounds2plates.sh
Alternatively:
bash Split_seq_V7.7.nova.walk_bothgenome.2rounds2plates.sh
The default configuration uses:
demo_fastq/
as input,
config_demo.yaml
as the sample configuration, and
output/
as the output directory.

# Example
A minimal DipC2/RNA2 configuration is:
Project=(Demo_DNA Demo_RNA)
Type=(DipC2 RNA2)

Genomes=(hg38 hg38)
ReadsPerBarcode=(100 100)

Start=Fastq
Runtype=full
chem=rev

cores=32
taps=(F F)
N6=(T T)
cleanup=F
with:
Project1:
  Name: Demo_DNA
  Type: DipC2
  P5: P1.29
  P7: P2.29

Project2:
  Name: Demo_RNA
  Type: RNA2
  P5: P1.30
  P7: P2.30
  N6: T
Run with:
export BASEDIR=/path/to/Split-seq
conda activate splitseq

bash Split_seq_V7.7.nova.walk_bothgenome.2rounds2plates.sh

# Pipeline output
After successful processing, results for each project are moved into a directory corresponding to the project name.

For example:
output/
├── Demo_DNA/
├── Demo_RNA/
└── Useful/

DipC2
Important DipC2 outputs include processed chromatin-contact files and contact maps such as:
*.pairs.unique.final.gz
*.contacts.final.gz
*.cool
*.mcool
*.hic
*.filtered.counts.csv.gz
The pipeline also generates contact statistics and QC outputs.

RNA2
Important RNA2 outputs include:
*.gene.bc.matrices.h5
*.filtered.counts.csv.gz
along with RNA QC plots and intermediate UMI/gene-count files.

Useful
Selected plots and summary files are collected into:
output/Useful/

# Cleanup
Intermediate files can consume substantial disk space.

To retain intermediate files for troubleshooting:
cleanup=F
After confirming that the pipeline runs correctly, cleanup can be enabled:
cleanup=T
For a first run on a new system, cleanup=F is recommended.

# Notes
Paths containing spaces are not recommended.
Genome indices can be large and should generally not be committed to the GitHub repository.
Users should activate the appropriate Conda/software environment before starting the pipeline.
Project names and assay types in the Bash configuration and YAML file must agree.
DipC2 and RNA2 refer to the two-round barcoding workflow implemented in this version of the pipeline.

# Citation
If you use this pipeline, please cite the relevant ChromE-seq and experimental-method publications associated with your dataset and analysis.

# Contact
For questions about the ChromE-seq alignment pipeline:
Sai Ma (masai.zju@gmail.com)
Zhicong Liao (zhicongliao@gmail.com)
For questions about modifications in this repository, please open a GitHub issue.









