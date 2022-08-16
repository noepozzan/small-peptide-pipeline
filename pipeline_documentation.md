# SMAPP: workflow documentation

This document describes the individual steps of the workflow. For instructions
on installation and usage please see [here](README.md).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Third-party software used](#third-party-software-used)
- [Description of workflow steps](#description-of-workflow-steps)
  - [Rule graph](#rule-graph)
  - [Understanding the config files](#understanding-the-config-files)
  - [Parameters](#parameters)
    - [Input Files](#input-files)
    - [Parameter table](#parameter-table)
  - [Profiles](#profiles)
    - [`full`](#full)
    - [`test`](#test)
    - [`proteomics`](#proteomics)
    - [`prepare`](#prepare)
    - [`genome_map`](#genome)
    - [`qc`](#qc)
    - [`ribotish`](#ribotish)
  - [Subworkflows](#subworkflows)
    - [`PREPARE`](#PREPARE)
- [FAQ](#faq)

## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **bedtools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][code-bedtools] |
| **cutadapt** | [MIT][license-mit] | _"[...] finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads"_ | [code][code-cutadapt] / [manual][docs-cutadapt] / [publication][pub-cutadapt] |
| **gffread** | [MIT][license-mit] | _"[...] validate, filter, convert and perform various other operations on GFF files"_ | [code][code-gffread] / [manual][docs-gffread] |
| **FastQC** | [GPLv3][license-gpl3] | _"A quality control analysis tool for high throughput sequencing data"_ | [code][code-fastqc] / [manual][docs-fastqc] |
| **MultiQC** | [GPLv3][license-gpl3] | _"Aggregate results from bioinformatics analyses across many samples into a single report"_ | [code][code-multiqc] / [manual][docs-multiqc] / [publication][pub-multiqc] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |
| **STAR** | [MIT][license-mit] | _"**S**pliced **T**ranscripts **A**lignment to a **R**eference"_ - _"RNA-seq aligner"_ | [code][code-star] / [manual][docs-star] / [publication][pub-star] |
| **Ribo-TISH** | [GPLv3][license-gpl3] | _"Ribo-TISH: Ribo-seq data-driven Translation Initiation Sites Hunter"_ | [code][code-ribotish] / [manual][docs-ribotish] / [publication][pub-ribotish] |
| **segemehl** | [GPLv3][license-gpl3] | _"segemehl is a software to map short sequencer reads to reference genomes."_ | [code][code-segemehl] / [manual][docs-segemehl] / [publication][pub-segemehl] |
| **Philosopher** | [GPLv3][license-gpl3] | _"A complete toolkit for shotgun proteomics data analysis"_ | [code][code-philosopher] / [manual][docs-philosopher] / [publication][pub-philosopher] |
| **MSFragger** | [msfragger licencse][license-msfragger] | _"Ultrafast, comprehensive peptide identification for mass spectrometry–based proteomics"_ | [code][code-msfragger] / [manual][docs-msfragger] / [publication][pub-msfragger] |
| **FASTX** | [GPLv3][license-gpl3] | _"The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing."_ | [code][code-fastx] / [manual][docs-fastx] / [publication][pub-fastx] |


## Description of workflow steps

> The workflow consists of multiple Nextflow files: A main `Nextflow` and
> individual subworkflow files for each configuration mode 
> The `main.nf` file contains the general steps .... Individual steps of the workflow are described briefly, and
> links to the respective software manuals are given. Parameters that can be
> modified by the user (via the samples table) are also described. Descriptions
> for steps for which individual "rules" exist for single- and paired-end
> sequencing libraries are combined, and only differences between the modes are
> highlighted.

### Rule graph

![rule_graph][rule-graph]

Visual representation of workflow. Automatically prepared with [Nextflow][docs-nextflow].

## Understanding the `.config` files

Since this workflow can be run in different configurations, the requirements differ.

After choosing the `profile` that fits your analysis needs, you need to do _2_ things:
- Most importantly: change the paths to your own files under `conf/params/`
- Change the environment settings under `conf/envs/`

#### The files under `conf/params/` look like this:

- The first `params` part of the files is made up of the file path specifications.
- The second `params` part is made up of the [parameters](#parameters) the tools of the pipeline use.
- The third `params` part specifies the output paths of your results.
- The last `params` part contains variables important for the workflow that you probably shouldn't change.

https://github.com/noepozzan/small-peptide-pipeline/blob/4af1760caecbbdf82f1168000edad05b5a7b1003/conf/params/qc.config#L1-L30

#### The files under `conf/envs/` look like this:

- The first part of the files tell the virtualization systems what to do. You probably don't have to change anything about this.
- The `params` part contains [parameters](#parameters) that will be highlighted below.
- The first `process` part defines your system specifications of `memory` and `cpus`.
- The second `process` part contains the specifications for the Docker images used. Please do **not** change this.

https://github.com/noepozzan/small-peptide-pipeline/blob/4af1760caecbbdf82f1168000edad05b5a7b1003/conf/envs/docker.config#L1-L108

## Parameters

Parameters are found in the `conf` and the `nextflow.config`.
As explained earlier, each `profile` file contains the parameters necessary to run the workflow.
The `nextflow.config` contains very few parameters that should probably not be changed.

The table below will hopefully give you some deeper understanding of all parameters used:

### Input Files

Parameter name | Description | Data type(s)
--- | --- | ---
riboseq_reads | Main ribosome sequencing fastq.gz files | `str`
proteomics_reads | Main proteomics mzML files | `str`
gtf | Main gene annotation file in gene transfer format | `str`
other_RNAs_sequence | Fasta file of rRNA library to be used for filtering out riboseq noise | `str`
genome | Fasta genome file to be used as your reference genome | `str`

The rest of the parameters are either important for the logic of the workflow or are input parameters to tools run by the pipeline.
These are the default parameters that are sure to work with the pipeline. The tools used for analying the data our way obviously contain way more parameters. If you are familiar with the tools, feel free to change these parameters or pass them other parameters, too!
If you want to know more about the tools (and the specific version you are running thereof), you may do the following:

```bash
singularity shell ~/.singularity/cache/library/<tool you are running>
<tool you are running> -h
```

### Parameter table

Parameter name | Description | Process | Data type(s)
--- | --- | --- | ---
singularity_store | Specifies path of the folder where singularity images get stored | General | `str`
cut | Required for [Cutadapt](#third-party-software-used), _Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning._ | from `TRIM_FIRST_BASES` process in `subworkflows/prepare_reads.nf` | `int`
minimum_length | Required for [Cutadapt](#third-party-software-used), _Discard reads shorter than LENGTH. Default: 0_ | from `TRIM_FIRST_BASES` process in `subworkflows/prepare_reads.nf` | `int`
clip_reads_args | Required for [FASTX](#third-party-software-used), contains parameters described right below (tool used: fastx_clipper)| from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `multiple`
-v | Required for [FASTX](#third-party-software-used), _Verbose - report number of sequences._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `none`
-n | Required for [FASTX](#third-party-software-used), _discard sequences shorter than N nucleotides. default is 5._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `int`
-l | Required for [FASTX](#third-party-software-used), _keep sequences with unknown (N) nucleotides. default is to discard such sequences._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `none`
-c | Required for [FASTX](#third-party-software-used), _Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter)._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `none`
-z | Required for [FASTX](#third-party-software-used), _Compress output with GZIP._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `none`
-a | Required for [FASTX](#third-party-software-used), _ADAPTER string. default is CCTTAAGG (dummy adapter)._ | from `CLIP_READS` process in `subworkflows/prepare_reads.n` | `str`
trim_reads_args | Required for [FASTX](#third_party_software_used), contains parameters described right below (tool used: fastq_quality_trimmer) | from `TRIM_READS` process in `subworkflows/prepare_reads.nf` | `multiple`
-v | Required for [FASTX](#third-party-software-used), _Verbose - report number of sequences._ | from `TRIM_READS` process in `subworkflows/prepare_reads.n` | `none`
-l | Required for [FASTX](#third-party-software-used), _Minimum length - sequences shorter than this (after trimming) will be discarded. Default = 0 = no minimum length._ | from `TRIM_READS` process in `subworkflows/prepare_reads.n` | `int`
-t | Required for [FASTX](#third-party-software-used), _Quality threshold - nucleotides with lower quality will be trimmed (from the end of the sequence)._ | from `TRIM_READS` process in `subworkflows/prepare_reads.n` | `int`
-z | Required for [FASTX](#third-party-software-used), _Compress output with GZIP._ | from `TRIM_READS` process in `subworkflows/prepare_reads.n` | `none`
filter_reads_args | Required for [FASTX](#third_party_software_used), contains parameters described right below (tool used: fastq_quality_filter) | from `FILTER_READS` process in `subworkflows/prepare_reads.nf` | `multiple`
-v | Required for [FASTX](#third-party-software-used), _Verbose - report number of sequences._ | from `FILTER_READS` process in `subworkflows/prepare_reads.n` | `none`
-q | Required for [FASTX](#third-party-software-used), _Minimum quality score to keep._ | from `FILTER_READS` process in `subworkflows/prepare_reads.n` | `int`
-p | Required for [FASTX](#third-party-software-used), _Minimum percent of bases that must have [-q] quality._ | from `FILTER_READS` process in `subworkflows/prepare_reads.n` | `int`
-z | Required for [FASTX](#third-party-software-used), _Compress output with GZIP._ | from `FILTER_READS` process in `subworkflows/prepare_reads.n` | `none`
fastq_to_fasta_args | Required for [FASTX](#third_party_software_used), contains parameters described right below (tool used: fastq_to_fasta) | from `FASTQ_TO_FASTA` process in `subworkflows/prepare_reads.nf` | `multiple`
-v | Required for [FASTX](#third-party-software-used), _Verbose - report number of sequences._ | from `FASTQ_TO_FASTA` process in `subworkflows/prepare_reads.n` | `none`
-n | Required for [FASTX](#third-party-software-used), _keep sequences with unknown (N) nucleotides. Default is to discard such sequences._ | from `FASTQ_TO_FASTA` process in `subworkflows/prepare_reads.n` | `none`
-r | Required for [FASTX](#third-party-software-used), _Rename sequence identifiers to numbers._ | from `FASTQ_TO_FASTA` process in `subworkflows/prepare_reads.n` | `none`
segemehl_args | Required for [segemehl](#third_party_software_used), contains parameters described right below (tool used: segemehl.x) | from `MAP_TRANSCRIPTOME_SEGEMEHL` process in `subworkflows/map_transcriptome.nf` and `MAP_rRNA_SEGEMEHL` in `subworkflows/map_rrna.nf` | `multiple`
--silent | Required for [segemehl](#third_party_software_used), _shut up!_ | from `MAP_TRANSCRIPTOME_SEGEMEHL` process in `subworkflows/map_transcriptome.nf` and `MAP_rRNA_SEGEMEHL` in `subworkflows/map_rrna.nf` | `none`
--accuracy | Required for [segemehl](#third_party_software_used), _min percentage of matches per read in semi-global alignment (default:90)_ | from `MAP_TRANSCRIPTOME_SEGEMEHL` process in `subworkflows/map_transcriptome.nf` and `MAP_rRNA_SEGEMEHL` in `subworkflows/map_rrna.nf` | `int`
--threads | Required for [segemehl](#third_party_software_used), _start <n> threads (default:1)_ | from `MAP_TRANSCRIPTOME_SEGEMEHL` process in `subworkflows/map_transcriptome.nf` and `MAP_rRNA_SEGEMEHL` in `subworkflows/map_rrna.nf` | `int`
star_map_threads | Required for [STAR](#third_party_software_used), _int: number of threads to run STAR_ | from `MAP_GENOME_STAR` process in `subworkflows/map_genome.nf` | `int`
check_peridocitiy_codnum | Required for [python script](#third_party_software_used), _Input codon coverage_ | from `CHECK_PERIODICITY` process in `subworkflows/qc.nf` | `int`
riboseq_mode | Required for [Ribo-TISH](#third_party_software_used), _choose `regular` if you have ordinary riboseq bam files, choose `TI` if you have TIS enriched riboseq bam files_ | from `RIBOTISH_QUALITY` and `RIBOTISH_PREDICT` in `subworkflows/ribotish.nf` | `str`
ribotish_quality_th | Required for [Ribo-TISH](#third_party_software_used), _Threshold for quality (default: 0.5)_ | from `RIBOTISH_QUALITY` in `subworkflows/ribotish.nf` | `float`
ribotish_predict_mode | Required for [Ribo-TISH](#third_party_software_used), _if `longest` chosen: Only report longest possible ORF results_ | from `RIBOTISH_PREDICT` in `subworkflows/ribotish.nf` | `str`
workspace | Required for [Philosopher](#third_party_software_used), _path where you want [Philosopher](#third_party_software_used) to do its tasks_ | from multiple processes in `subworkflows/ribotish.nf` | `str`
peptideprophet_args | Required for [Philosopher](#third_party_software_used), contains parameters described right below (tool used: philosopher peptideprophet) | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `multiple`
--combine | Required for [Philosopher](#third_party_software_used), _combine the results from PeptideProphet into a single result file_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
--decoy | Required for [Philosopher](#third_party_software_used), _semi-supervised mode, protein name prefix to identify decoy entries_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `str`
--ppm | Required for [Philosopher](#third_party_software_used), _use ppm mass error instead of Daltons for mass modeling_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
--accmass | Required for [Philosopher](#third_party_software_used), _use accurate mass model binning_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
--expectscore | Required for [Philosopher](#third_party_software_used), _use expectation value as the only contributor to the f-value for modeling_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
--decoyprobs | Required for [Philosopher](#third_party_software_used), _compute possible non-zero probabilities for decoy entries on the last iteration_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
--nonparam | Required for [Philosopher](#third_party_software_used), _use semi-parametric modeling, must be used in conjunction with --decoy option_ | from `PEPTIDEPROPHET` process in `subworkflows/philosopher.nf` | `none`
philosopher_filter_args | Required for [Philosopher](#third_party_software_used), contains parameters described right below (tool used: philosopher filter) | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `multiple`
--psm | Required for [Philosopher](#third_party_software_used), _psm FDR level (default 0.01)_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `float`
--ion | Required for [Philosopher](#third_party_software_used), _peptide ion FDR level (default 0.01)_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `float`
--pep | Required for [Philosopher](#third_party_software_used), _peptide FDR level (default 0.01)_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `float`
--prot | Required for [Philosopher](#third_party_software_used), _protein FDR level (default 0.01)_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `float`
--picked | Required for [Philosopher](#third_party_software_used), _apply the picked FDR algorithm before the protein scoring_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `none`
--tag | Required for [Philosopher](#third_party_software_used), _decoy tag (default "rev_")_ | from `FILTER_FDR` process in `subworkflows/philosopher.nf` | `str`
  
## Profiles

Different combinations of subworkflows make up the profiles. Profiles are the units of the pipeline that are reasonable to run by themselves since they do some bioinformatic process of interest, such as preparing reads or performing quality control.
  
### `full`
This profile contains all of the subworkflows in the order they are listed, beginning from [`PREPARE`](#prepare).

### `test`
This profile contains all of the subworkflows in the order they are listed, beginning from [`PREPARE`](#prepare), but **WITHOUT** the [**RIBOTISH**](#ribotish) subworkflow, since this tool cannot be tested with small data.

### `proteomics`
This profile only contains the [`PHILOSOPHER`](#philosopher) subworkflow.

### `prepare`
This profile only contains the [`PREPARE`](#prepare) subworkflow.

### `ribotish`
This profile only contains the [`RIBOTISH`](#ribotish) subworkflow.

### `genome_map`
This profile contains the [`PREPARE`](#prepare), [`rRNA`](#rrna) and [`GENOME`](#genome) subworkflows.

### `qc`
This profile contains the [`PREPARE`](#prepare), [`rRNA`](#rrna), [`ANNOTATE`](#annotate), [`TRANSCRIPTOME`](#transcriptome) and [`QC`](#qc) subworkflows.

## Subworkflows

Subworkflows are parts of the larger, complete workflow that do one of the tasks of the larger worfklow.
They may be considered as being a self-contained unit of the larger workflow and may thus be used independently of it.
  
### `PREPARE`
This subworkflow is made up of 5 [Nextflow](#third-party-software-used) processes.
It takes the `.fastq.gz` files you've gotten back from your ribosome sequencing experiments as input and returns trimmed, clipped and filtered `.fasta` files.
  
#### `TRIM_FIRST_BASES`
_"[cutadapt](#third-party-software-used) removes adapter sequences from high-throughput sequencing reads."_
- **Input**
  - Ribosome sequencing (`.fastq.gz`) files
- **Parameters**
  - `--cut`: _Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning._
  - `--minimum_length`: _Discard reads shorter than LENGTH. Default: 0_
- **Output**
  - Ribosome sequencing files with first bases trimmed; used in [**CLIP_READS**](#clip_reads)

#### `CLIP_READS`
_"Removing sequencing adapters / linkers"_
- **Input**
  - riboseq files with first bases trimmed; from [**TRIM_FIRST_BASES**](#trim_first_bases)
- **Parameters**
  - `-v`: _Verbose - report number of sequences._
  - `-n`: _discard sequences shorter than N nucleotides. default is 5._
  - `-l`: _keep sequences with unknown (N) nucleotides. default is to discard such sequences._
  - `-c`: _Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter)._
  - `-a`: _ADAPTER string. default is CCTTAAGG (dummy adapter)._
- **Output**
  - riboseq files with adapters removed; used in [**TRIM_READS**](#trim_reads)
- **Non-configurable & non-default**
  - `-z`: _Compress output with GZIP._

#### `TRIM_READS`
_"Trims (cuts) sequences based on quality"_
- **Input**
  - riboseq files with adapters removed; from [**CLIP_READS**](#clip_reads)
- **Parameters**
  - `-v`: _Verbose - report number of sequences._
  - `-l`: _Minimum length - sequences shorter than this (after trimming) will be discarded. Default = 0 = no minimum length._
  - `-t`: _Quality threshold - nucleotides with lower quality will be trimmed (from the end of the sequence)._
- **Output**
  - trimmed riboseq files with adapters removed; used in [**FILTER_READS**](#filter_reads)
- **Non-configurable & non-default**
  - `-z`: _Compress output with GZIP._
  
#### `FILTER_READS`
_"Filters sequences based on quality"_
- **Input**
  - trimmed riboseq files with adapters removed; from [**TRIM_READS**](#trim_reads)
- **Parameters
  - `-v`: _Verbose - report number of sequences._
  - `-q`: _Minimum quality score to keep._
  - `-p`: _Minimum percent of bases that must have [-q] quality._
- **Output**
  - trimmed and filtered riboseq files with adapters removed; used in [**FASTQ_TO_FASTA**](#fastq_to_fasta)
- **Non-configurable & non-default**
  - `-z`: _Compress output with GZIP._

#### `FASTQ_TO_FASTQ`
_"Convert FASTQ files to FASTA files."_
- **Input**
  - trimmed and filtered riboseq files with adapters removed; from [**FILTER_READS**](#filter_reads)
- **Parameters**
  - `-v`: _Verbose - report number of sequences._
  - `-n`: _keep sequences with unknown (N) nucleotides. Default is to discard such sequences._
  - `-r`: _Rename sequence identifiers to numbers._
- **Output**
  - trimmed and filtered riboseq fasta files without adapters; used in [**rRNA**](#rrna), 
  
### `rRNA`
This subworkflow is made up of 2 [Nextflow](#third-party-software-used) processes.
It takes the trimmed and filtered riboseq fasta file from [**PREPARE**](#prepare) as input and returns mapped and unmapped reads against a reference of rRNA.

#### `SEGEMEHL_INDEX_rRNA`
Index rRNA reference fasta file
- **Input**
  - rRNA reference fasta file
- **Output**
  - rRNA reference index; used in [**MAP_rRNA_SEGEMEHL**](#map_rrna_segemehl)

#### `MAP_rRNA_SEGEMEHL`
Map prepared riboseq reads to rRNA reference
- **Input**
  - prepared fasta reads; from [**PREPARE**](#prepare)
  - rRNA reference index; from [**SEGEMEHL_INDEX_rRNA**](#segemehl_index_rrna)
  - rRNA reference fasta file
- **Parameters**
  - `--silent`: _shut up!_
  - `--accuracy`: _min percentage of matches per read in semi-global alignment (default:90)_
  - `--threads`: _start <n> threads (default:1)_
- **Output**
  - Reads mapped to rRNA reference
  - Reads unmapped to rRNA reference; used in [**GENOME**](#genome), [`TRANSCRIPTOME`](#transcriptome) and [`QC`](#qc)

### `GENOME`
This subworkflow is made up of 3 [Nextflow](#third-party-software-used) processes.
It takes the unmapped reads from [**rRNA**](#rna), the reference genome and a genome annotation file as inputs and returns mapped and unmapped reads against the genome.
  
#### `STAR_INDEX_GENOME`
Index genome
- **Input**
  - reference genome
- **Parameters**
  - `star_map_threads`: _int: number of threads to run STAR_
- **Output**
  - genome index; used in [**MAP_GENOME_STAR**](#map_genome_star)

#### `MAP_GENOME_STAR`
Map the reads that didn't map to any rRNA to the genome
- **Input**
  - Reads unmapped to rRNA reference; from [**rRNA**](#rrna)
  - genome index; from [**STAR_INDEX_GENOME**](#star_index_genome)
  - genome annotation (`.gtf`)
- **Parameters**
  - `star_map_threads`: _int: number of threads to run STAR_
- **Output**
  - Reads mapped to genome; used in [**SAM_TO_BAM_SORT_AND_INDEX_STAR**](#sam_to_bam_sort_and_index_star)
- **Non-configurable & non-default**
  - `--outSAMattributes All`
  - `--quantMode GeneCounts`
  - `--outReadsUnmapped Fastx`

#### `SAM_TO_BAM_SORT_AND_INDEX_STAR`
Compress sam files to bam, then sort and index them and pack these sorted bam+index into folders for easier downstream handling
- **Input**
  - Reads mapped to genome; from [**MAP_GENOME_STAR**](#map_genome_star)
- **Output**
  - sorted bam and its index in folder; used in [**RIBOTISH**](#ribotish)

### `ANNOTATE`
This subworkflow is made up of 4 [Nextflow](#third-party-software-used) processes.
It takes the rRNA reference, the reference genome and a genome annotation file as inputs and returns a fasta file containing the longest transcripts and a tsv file containing "transcript id, gene id, CDS"(MERIC?)
  
#### `SELECT_LONGEST_CODING_TRANSCRIPT`
Python script that selects longest coding transcript per gene
- **Input**
  - genome annotation file (`.gtf`)
  - python script
- **Output**
  - gtf file containing only the longest coding transcript for each gene; used in [**EXTRACT_TRANSCRIPT_SEQUENCES**](#extract_transcript_sequences) and [**CREATE_TAB_DELIMITED_CDS_FILE**](#create_tab_delimited_cds_file)

#### `EXTRACT_TRANSCRIPT_SEQUENCES`
Extracts the longest transcripts for each gene into a fasta file
- **Input**
  - gtf file containing only the longest coding transcript for each gene; from [**SELECT_LONGEST_CODING_TRANSCRIPT**](#select_longest_coding_transcript)
  - reference genome
- **Output**
  - fasta file containing longest coding transcript for each gene; used in [**TRANSCRIPTOME**](#transcriptome) and [**CREATE_TAB_DELIMITED_CDS_FILE**](#create_tab_delimited_cds_file)
  
#### `CREATE_TAB_DELIMITED_CDS_FILE`
Generate transcript id, gene id, CDS table
- **Input**
  - gtf file containing only the longest coding transcript for each gene; from [**SELECT_LONGEST_CODING_TRANSCRIPT**](#select_longest_coding_transcript)
  - fasta file containing longest coding transcript for each gene; from [**EXTRACT_TRANSCRIPT_SEQUENCES**](#extract_transcript_sequences)
  - python script
- **Output**
  - tsv file containing transcript id, gene id, CDS; used in [**CREATE_BED_CDS_FILE**](#create_bed_cds_file) and [**QC**](#qc)

#### `CREATE_BED_CDS_FILE`
Reorder above tsv file to be in `.bed` format
- **Input**
  - tsv file containing transcript id, gene id, CDS; from [**CREATE_TAB_DELIMITED_CDS_FILE**](#create_tab_delimited_cds_file)
- **Output**
  - bed file containing transcript id, CDS, gene id
  
### `TRANSCRIPTOME`
  
### `QC`

### `RIBOTISH`

This subworfklow is made up of 6 nextflow processes.
As inputs, it takes a gtf file, a fasta genome file and folders containing both an indexed+sorted bam file and its index (from [**GENOME**](#genome)).
As output, the workflow returns a fasta file with the predicted peptides found by ribosome sequencing.
  
If you have the right data and want to predict translation inititation sites (TIS), please select "TI" for the parameter `riboseq_mode`.  
If you have regular ribosome sequencing data and want to predict open reading frames, please select "regular" for the parameter `riboseq_mode`.

#### `FASTA_INDEX`
Create genome index using [**SAMtools**](#third-party-software-used).
- **Input**
  - Genome sequence file (`.fa`)
- **Output**
  - Genome index file (`.fai`); used in [**GFFREAD**](#gffread)

#### `RIBOTISH_QUALITY`
Create the ribosome offsets to determine the proper position of the ribosome on the reads. Uses [**Ribo-TISH**](#third-party-software-used)
- **Input**
  - Folders containing both an indexed+sorted bam file and its index; from [**GENOME**](#genome)
  - GTF file (`.gtf`)
- **Parameters**
  - `ribotish_quality_th`: _Threshold for quality (default: 0.5)_
- **Output**
  - offsets; used in [**RIBOTISH_PREDICT**](#ribotish_predict)

#### `RIBOTISH_PREDICT`
Go from bam files of aligned ribosome sequencing reads to new, short peptides.
- **Input**
  - Folders containing both an indexed+sorted bam file and its index; from [**GENOME**](#genome)
  - offsets: from [**RIBOTISH_QUALITY**](#ribotish_quality)
  - GTF file (`.gtf`)
  - Genome sequence file (`.fa`)
- **Parameters**
  - `ribotish_predict_mode`: _if `longest` chosen: Only report longest possible ORF results_
- **Output**
  - predicted peptides in nucleotide fasta format; used in [**SORF_TO_PEPTIDE**](#sorf_to_peptide)

#### `GFFREAD`
Extract the transcriptome using [**GFFREAD**](#third-party-software-used)
- **Input**
  - Genome sequence file (`.fa`)
  - Genome index file (`.fai`); from [**FASTA_INDEX**](#fasta_index)
  - GTF file (`.gtf`)
- **Output**
  - Extracted transcripts in tabular format; used in [**SORF_TO_PEPTIDE**](#sorf_to_peptide)

#### `SORF_TO_PEPTIDE`
Translates the short open reading frames into peptide sequences.
- **Input**
  - predicted peptides in nucleotide fasta format; from [**RIBOTISH_PREDICT**](#ribotish_predict)
  - Extracted transcripts in fasta format; from [**GFFREAD**](#gffread)
  - script for sorf to peptide translation
- **Output**
  - predicted peptides in fasta aa format; used in [**COMBINE**](#combine)

#### `COMBINE`
Collect different files of predicted peptides into 1, using simple bash
- **Input**
  - predicted peptides in fasta aa format; from [**SORF_TO_PEPTIDE**](#sorf_to_peptide)
- **Output**
  - combined predicted peptides in aa format: used in [**PHILOSOPHER**](#philosopher)

### `PROTEOMICS`
This subworfklow is made up of 10 nextflow processes.
As inputs, it takes the predicted peptides in aa format (from [**COMBINE**](#combine)) and mzML proteomics spectra files.
As output, the workflow returns a csv file containing the proteomics-validated small peptides initially found by ribosome sequencing.

#### `WORKSPACE`
Initializes a physical directory with a hidden folder which is needed for [**Philosopher**](#third-party-software-used) to work.
- **Input**
  - predicted peptides; from [**COMBINE**](#combine), works as pseudo input to make pipeline run in the right order
- **Parameters**
  - `workspace`: _path where you want [Philosopher](#third_party_software_used) to do its tasks_
- **Output**
  - pseudo output; to be used in [**DATABASE**](#database)

#### `DATABASE`
Add decoys and contaminants to predicted peptides fasta file and format it for philosopher
- **Input**
  - Pseudo input from [**WORKSPACE**](#workspace)
  - predicted peptides; from [**COMBINE**](#combine)
- **Output**
  - philosopher database fasta file (`.fas`); used in [**GENERATE_CHANGE_PARAMS**](#generate_change_params), [**MSFRAGGER**](#msfragger) and [**PEPTIDEPROPHET**](#peptideprophet)

#### `GENERATE_CHANGE_PARAMS`
Generate a parameter file necessary for [**MSFRAGGER**](#msfragger) and change some parameters using a python script
- **Input**
  - database fasta file; from [**DATABASE**](#database)
  - python script to change params
- **Output**
  - msfragger parameter file with changed params; used in [**MSFRAGGER**](#msfragger)
- **Non-configurable & non-default**
  - `calibrate_mass = 0`: _"[...] speed up the search even more."_

#### `MSFRAGGER`
Search fasta database against mzML spectra for peptides that appear in both ("hits").
- **Input**
  - database fasta file; from [**DATABASE**](#database)
  - msfragger parameter file; from  [**GENERATE_CHANGE_PARAMS**](#generate_change_params)
  - mzML proteomics spectra files
- **Output**
  - (`.pepXML`) file containing hits; used in [**PEPTIDEPROPHET**](#peptideprophet) and [**IONQUANT**](#ionquant)
  
#### `PEPTIDEPROPHET`
Validates the peptide assignment
- **Input**
  - database fasta file; from [**DATABASE**](#database)
  - (`.pepXML`) file containing hits; from [**MSFRAGGER**](#msfragger)
- **Output**
  - `.xml` file containing validated peptides; used in [**PROTEINPROPHET**](#proteinprophet) and [**FILTER_FDR**](#filter_fdr)
- **Non-configurable & non-default**
  - `--combine`: _combine the results from PeptideProphet into a single result file_ 
  - `--decoy`: _semi-supervised mode, protein name prefix to identify decoy entries_
  - `--ppm`: _use ppm mass error instead of Daltons for mass modeling_
  - `--accmass`: _use accurate mass model binning_
  - `--expectscore`: _use expectation value as the only contributor to the f-value for modeling_
  - `--decoyprobs`: _compute possible non-zero probabilities for decoy entries on the last iteration_
  - `--nonparam`: _use semi-parametric modeling, must be used in conjunction with --decoy option_
  
#### `PROTEINPROPHET`
Performs protein inference (skipped for now in our case since we are only interested in peptides)
- **Input**
  - `.xml` file containing validated peptides; from [**PEPTIDEPROPHET**](#peptideprophet)
- **Output**
  - `.xml` file containing validated proteins; used in [**FILTER_FDR**](#filter_fdr)
#### `FILTER_FDR`
Filter the hits by specifying cutoff FDR value
- **Input**
  - `.xml` file containing validated peptides; from [**PEPTIDEPROPHET**](#peptideprophet)
  - `.xml` file containing validated proteins; from [**PROTEINPROPHET**](#proteinprophet)
- **Parameters**
  - `--psm`: _psm FDR level (default 0.01)_
  - `--ion`: _peptide ion FDR level (default 0.01)_
  - `--pep`: _peptide FDR level (default 0.01)_
  - `--prot`: _protein FDR level (default 0.01)_
- **Output**
  - pseudo output (since results of process get written into hidden dir); used in [**FREEQUANT**](#freequant)
- **Non-configurable & non-default**
  - `--picked`: _apply the picked FDR algorithm before the protein scoring_
  - `--tag`: _decoy tag (default "rev\_")_
#### `FREEQUANT`
Quantification using MS1 peak intensities
- **Input**
  - pseudo input; from [**FILTER_FDR**](#filter_fdr)
- **Output**
  - pseudo output; used in [**REPORT**](#report)
#### `REPORT`
Inspect the results
- **Input**
  - pseudo input; from [**FREEQUANT**](#freequant)
- **Output**
  - result files of philosopher (*.tsv, msstats.csv, ..)
#### `IONQUANT`
Quantification of MS1-precursor intensity
- **Input**
  - pseudo input; from [**REPORT**](#report)
  - `.xml` files containing hits; from [**MSFRAGGER**](#msfragger)
- **Output**
  - `.csv` files containing quantified hits

#### `map_genome_star`

Align short reads to reference genome and/or transcriptome with
[**STAR**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Index; from [**create_index_star**](#create_index_star)
- **Parameters**
  - **rule_config.yaml**
    - `--outFilterMultimapScoreRange=0`: the score range below the maximum score for multimapping alignments (default 1)
    - `--outFilterType=BySJout`: reduces the number of ”spurious” junctions
    - `--alignEndsType`: one of `Local` (standard local alignment with soft-clipping allowed) or `EndToEnd` (force end-to-end read alignment, do not soft-clip); specify in sample table column `soft_clip`
    - `--twopassMode`: one of `None` (1-pass mapping) or `Basic` (basic 2-pass mapping, with all 1st-pass junctions inserted into the genome indices on the fly); specify in sample table column `pass_mode`
    - `--outFilterMultimapNmax`: maximum number of multiple alignments allowed; if exceeded, read is considered unmapped; specify in sample table column `multimappers`
- **Output**
  - Aligned reads file (`.bam`); used in
    [**sort_genomic_alignment_samtools**](#sort_genomic_alignment_samtools),
  - STAR log file
- **Non-configurable & non-default**
  - `--outSAMattributes=All`: NH HI AS nM NM MD jM jI MC ch
  - `--outStd=BAM_Unsorted`: which output will be directed to `STDOUT` (default 'Log')
  - `--outSAMtype=BAM Unsorted`: type of SAM/BAM output (default SAM)
  - `--outSAMattrRGline`: ID:rnaseq_pipeline SM: *sampleID*
  
## FAQ

If you get an error message that looks like gibberish to you, have a look at the frequently asked questions below:
  
#### gzip: unexpected end of file
If your workflow run fails with the the error message containing this bit, please just restart your workflow run.
Nextflow seems to sometimes just be over-hasty and begin a process before having properly received the whole file to work on.

[code-bedtools]: <https://github.com/arq5x/bedtools2>
[code-cutadapt]: <https://github.com/marcelm/cutadapt>
[code-gffread]: <https://github.com/gpertea/gffread>
[code-fastqc]: <https://github.com/s-andrews/FastQC>
[code-multiqc]: <https://github.com/ewels/MultiQC>
[code-samtools]: <https://github.com/samtools/samtools>
[code-star]: <https://github.com/alexdobin/STAR>
[code-ribotish]: <https://github.com/zhpn1024/ribotish>
[code-segemehl]: <not known>
[code-philosopher]: <https://github.com/Nesvilab/philosopher>
[code-msfragger]: <https://github.com/Nesvilab/MSFragger>
[code-fastx]: <not known>

[docs-bedtools]: <https://bedtools.readthedocs.io/en/latest/>
[docs-cutadapt]: <https://cutadapt.readthedocs.io/en/stable/>
[docs-gffread]: <http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread>
[docs-fastqc]: <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/>
[docs-multiqc]: <https://multiqc.info/docs/>
[docs-samtools]: <http://www.htslib.org/doc/samtools.html>
[docs-nextflow]: <https://www.nextflow.io/docs/latest/index.html>
[docs-star]: <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
[docs-ribotish]: <https://github.com/zhpn1024/ribotish/blob/master/README.rst>
[docs-segemehl]: <http://www.bioinf.uni-leipzig.de/Software/segemehl/>
[docs-philosopher]: <https://github.com/Nesvilab/philosopher/wiki>
[docs-msfragger]: <https://github.com/Nesvilab/MSFragger/wiki>
[docs-fastx]: <http://hannonlab.cshl.edu/fastx_toolkit/>

[license-bsd2]: <https://opensource.org/licenses/BSD-2-Clause>
[license-gpl2]: <https://opensource.org/licenses/GPL-2.0>
[license-gpl3]: <https://opensource.org/licenses/GPL-3.0>
[license-mit]: <https://opensource.org/licenses/MIT>
[license-msfragger]: <http://msfragger-upgrader.nesvilab.org/upgrader/MSFragger-LICENSE.pdf>

[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-multiqc]: <https://doi.org/10.1093/bioinformatics/btw354>
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-star]: <https://doi.org/10.1093/bioinformatics/bts635>
[pub-ribotish]: <https://doi.org/10.1038/s41467-017-01981-8>
[pub-segemehl]: <https://doi.org/10.1371/journal.pcbi.1000502>
[pub-philosopher]: <https://doi.org/10.1038/s41592-020-0912-y>
[pub-msfragger]: <https://doi.org/10.1038/nmeth.4256>
[pub-fastx]: <not known>

[rule-graph]: images/flowchart.png
