# SMAPP: workflow documentation

This document describes the individual steps of the workflow. For instructions
on installation and usage please see [here](README.md).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Third-party software used](#third-party-software-used)
- [Description of workflow steps](#description-of-workflow-steps)
  - [Rule graph](#rule-graph)
  - [Preparatory](#preparatory)
    - [Understanding the config files](#understanding-the-config-files)
  - [Subworklows](#subworklows)
    - [`annotate`](#annotate)
    - [`qc`](#qc)
    - [`prepare`](#prepare)
    - [`map_genome`](#map_genome)
    - [`map_rrna`](#map_rrna)
    - [`map_transcriptome`](#map_transcriptome)
    - [`ribotish`](#ribotish)
    - [`philosopher`](#philosopher)

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


^ compatible with [GPLv3][license-gpl3]

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

### Preparatory

#### Understanding the `.config` files

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

### Parameters

Parameters are found in the `conf` and the `nextflow.config`.
As explained earlier, each `profile` file contains the parameters necessary to run the workflow.
The `nextflow.config` contains very few parameters that should probably not be changed.

The table below will hopefully give you some deeper understanding of all parameters used:

Parameter name | Description | Data type(s)
--- | --- | ---
riboseq_reads | Main ribosome sequencing fastq.gz files | `str`
proteomics_reads | Main proteomics mzML files | `str`
gtf | Main gene annotation file in gene transfer format | `str`
other_RNAs_sequence | Fasta file of rRNA library to be used for filtering out riboseq noise | `str`
genome | Fasta genome file to be used as your reference genome | `str`
  
fq1_3p | Required for [Cutadapt](#third-party-software-used). 3' adapter of mate 1. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. | `str`
fq1_5p | Required for [Cutadapt](#third-party-software-used). 5' adapter of mate 1. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. | `str`
fq2_5p | Required for [Cutadapt](#third-party-software-used). 5' adapter of mate 2. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. Value ignored for single-end libraries. | `str`
fq1_polya3p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 3' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. | `str`
fq1_polya5p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 5' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. | `str`
index_size | Required for [STAR](#third-party-software-used). Ideally the maximum read length minus 1. (`max(ReadLength)-1`). Values lower than maximum read length may result in lower mapping accuracy, while higher values may result in longer processing times. | `int`
kmer | Required for [Salmon](#third-party-software-used). Default value of 31 usually works fine for reads of 75 bp or longer. Consider using lower values if poor mapping is observed. | `int`
organism | Name or identifier of organism or organism-specific genome resource version. Has to correspond to the naming of provided genome and gene annotation files and directories, like "ORGANISM" in the path below. <br> **Example:** `GRCh38` | `str`
gtf | Required for [STAR](#third-party-software-used). Path to gene annotation `.gtf` file. File needs to be in subdirectory corresponding to `organism` field. <br> **Example:** `/path/to/GRCh38/gene_annotations.gtf` | `str`
genome | Required for [STAR](#third-party-software-used). Path to genome `.fa` file. File needs to be in subdirectory corresponding to `organism` field. <br> **Example:** `/path/to/GRCh38/genome.fa` | `str`
sd | Required for [kallisto](#third-party-software-used) and [Salmon](#third-party-software-used), but **only** for single-end libraries. Estimated standard deviation of fragment length distribution. Can be assessed from, e.g., BioAnalyzer profiles | `int`


### Sequencing mode-independent

#### `start`

Copy and rename read files.

> Local rule

- **Input**
  - Reads file (`.fastq.gz`)
- **Output**
  - Reads file, copied, renamed (`.fastq.gz`); used in [**fastqc**](#fastqc) and
    [**remove_adapters_cutadapt**](#remove_adapters_cutadapt)

#### `create_index_star`

Create index for [**STAR**](#third-party-software-used) short read aligner.

> Indices need only be generated once for each combination of genome, set of
> annotations and index size.

- **Input**
  - Genome sequence file (`.fasta`)
  - Gene annotation file (`.gtf`)
- **Parameters**
  - **samples.tsv**
    - `--sjdbOverhang`: maximum read length - 1; lower values may reduce accuracy, higher values may increase STAR runtime; specify in sample table column `index_size`
- **Output**
  - STAR index; used in [**map_genome_star**](#map_genome_star)
  - Index includes files:
    - Chromosome length table `chrNameLength.txt`; used in
      [**generate_alfa_index**](#generate_alfa_index) and
      [**prepare_bigWig**](#prepare_bigwig)
    - Chromosome name list `chrName.txt`; used in
      [**create_index_salmon**](#create_index_salmon)

#### `extract_transcriptome`

Create transcriptome from genome and gene annotations with
[**gffread**](#third-party-software-used).

- **Input**
  - Genome sequence file (`.fasta`)
  - Gene annotation file (`.gtf`)
- **Output**
  - Transcriptome sequence file (`.fasta`); used in
    [**concatenate_transcriptome_and_genome**](#concatenate_transcriptome_and_genome)
    and [**create_index_kallisto**](#create_index_kallisto)

#### `concatenate_transcriptome_and_genome`

Concatenate reference transcriptome and genome.

> Required by [Salmon][docs-salmon-selective-alignment]

- **Input**
  - Genome sequence file (`.fasta`)
  - Transcriptome sequence file (`.fasta`); from
    [**extract_transcriptome**](#extract_transcriptome)
- **Output**
  - Transcriptome genome reference file (`.fasta`); used in
    [**create_index_salmon**](#create_index_salmon)

#### `create_index_salmon`

Create index for [**Salmon**](#third-party-software-used) quantification.

> Required if Salmon is to be used in "mapping-based" mode.
>create_index_salmon
> Index is built using an auxiliary k-mer hash over k-mers of a specified
> length (default: 31). While the mapping algorithms will make use of
> arbitrarily long matches between the query and reference, the selected k-mer
> size will act as the minimum acceptable length for a valid match. Thus, a
> smaller value of k may slightly improve sensitivty. Empirically, the default
> value of k = 31 seems to work well for reads of 75 bp or longer. For shorter
> reads, consider using a smaller k.

- **Input**
  - Transcriptome genome reference file (`.fasta`); from
    [**concatenate_transcriptome_and_genome**](#concatenate_transcriptome_and_genome)
  - Chromosome name list `chrName.txt`; from
    [**create_index_star**](#create_index_star)
- **Parameters**
  - **samples.tsv**
    - `--kmerLen`: k-mer length; specify in sample table column `kmer`
- **Output**
  - Salmon index; used in [**quantification_salmon**](#quantification_salmon)

#### `create_index_kallisto`

Create index for [**kallisto**](#third-party-software-used) quantification.

> The default kmer size of 31 is used in this workflow and is not configurable
> by the user.

- **Input**
  - Transcriptome sequence file (`.fasta`); from
    [**extract_transcriptome**](#extract_transcriptome)
- **Output**
  - kallisto index; used in
    [**genome_quantification_kallisto**](#genome_quantification_kallisto)

#### `extract_transcripts_as_bed12`

Convert transcripts from `.gtf` to extended 12-column `.bed` format with
[custom-script][custom-script-gtf-to-bed12]. Note that the default transcript type setting is used, which is "protein_coding".

- **Input**
  - Gene annotation file (`.gtf`)
- **Output**
  - Transcript annotations file (12-column `.bed`); used in
    [**calculate_TIN_scores**](#calculate_tin_scores)

#### `fastqc`

Prepare quality control report for reads library with
[**FastQC**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from [**start**](#start)
- **Output**
  - FastQC output directory with report (`.txt`) and figures (`.png`); used in
    [**multiqc_report**](#multiqc_report)

#### `fastqc_trimmed`

Prepare quality control report for trimmed reads (after adapter and poly(A)-tail removal) with [**FastQC**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from [**remove_polya_cutadapt**](#remove_polya_cutadapt)
- **Output**
  - FastQC output directory with report (`.txt`) and figures (`.png`); used in
    [**multiqc_report**](#multiqc_report)

#### `sort_genomic_alignment_samtools`

Sort BAM file with [**SAMtools**](#third-party-software-used).

> Sort a genome aligned BAM file.

- **Input**
  - Alignemnts file (`.bam`); from [**map_genome_star**](#map_genome_star)
- **Output**
  - Alignemnts file (`.bam`); used in [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools) & [**star_rpm**](#star_rpm) & [**calculate_TIN_scores**](#calculate_tin_scores)

#### `index_genomic_alignment_samtools`

Index BAM file with [**SAMtools**](#third-party-software-used).

> Indexing a genome sorted BAM file allows one to quickly extract alignments
> overlapping particular genomic regions. Moreover, indexing is required by
> genome viewers such as IGV so that the viewers can quickly display alignments
> in a genomic region of interest.

- **Input**
  - Alignemnts file (`.bam`); from [**sort_genomic_alignment_samtools**](#sort_genomic_alignment_samtools)
- **Output**
  - BAM index file (`.bam.bai`); used in [**star_rpm**](#star_rpm) &
    [**calculate_TIN_scores**](#calculate_tin_scores)

#### `star_rpm`

Create stranded bedGraph coverage (`.bg`) with
[**STAR**](#third-party-software-used).

> Uses STAR's [RPM normalization][docs-star-rpm-norm] functionality.
>
> STAR RPM uses SAM flags to correctly tell where the read and its mate mapped
> to. That is, if mate 1 is mapped to the plus strand, then mate 2 is mapped to
> the minus strand, and STAR will count mate 1 and mate 2 to the plus strand.
> This is in contrast to `bedtools genomecov -bg -split`, where a mate is
> assigned to a strand irrespective of its corresponding mate.

- **Input**
  - Alignments file (`.bam`); from [**sort_genomic_alignment_samtools**](#sort_genomic_alignment_samtools)
  - BAM index file (`.bam.bai`); from
    [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
- **Output**
  - Coverage file (`.bg`); used in [**multiqc_report**](#multiqc_report) and
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)

#### `rename_star_rpm_for_alfa`

Rename and copy stranded bedGraph coverage tracks such that they comply with
[**ALFA**](#third-party-software-used).

> Local rule
>
> Renaming to `plus.bg` and `minus.bg` depends on library orientation, which is
> provided by user in sample table column `libtype`.

- **Input**
  - Coverage file (`.bg`); from [**star_rpm**](#star_rpm)
- **Output**
  - Coverage file, renamed (`.bg`); used in [**alfa_qc**](#alfa_qc) and
    [**sort_bed_4_big**](#sort_bed_4_big)

#### `sort_bed_4_big`

Sort bedGraph files with [**bedtools**](#third-party-software-used).

- **Input**
  - Coverage files, renamed (`.bg`); from
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
- **Output**
  - Coverage files, sorted (`.bg`); used in [**prepare_bigWig**](#prepare_bigwig)

#### `prepare_bigWig`

Generate bigWig from bedGraph with
[**bedGraphToBigWig**](#third-party-software-used).

- **Input**
  - Coverage files, sorted (`.bg`); from [**sort_bed_4_big**](#sort_bed_4_big)
  - Chromosome length table `chrNameLength.txt`; from
    [**create_index_star**](#create_index_star)
- **Output**
  - Coverage files, one per strand and sample (`.bw`); used in
    [**finish**](#finish)

#### `calculate_TIN_scores`

Calculates the Transcript Integrity Number (TIN) for each transcript with
[custom script][custom-script-tin] based on
[**RSeqC**](#third-party-software-used).

> TIN is conceptually similar to RIN (RNA integrity number) but provides
> transcript-level measurement of RNA quality and is more sensitive for
> low-quality RNA samples.
>
> - TIN score of a transcript reflects RNA integrity of the transcript
> - Median TIN score across all transcripts reflects RNA integrity of library
> - TIN ranges from 0 (the worst) to 100 (the best)
> - A TIN of 60 means that 60% of the transcript would have been covered if
>   the read coverage were uniform
> - A TIN of 0 will be assigned if the transcript has no coverage or coverage
>   is below threshold

- **Input**
  - Alignments file (`.bam`); from [**sort_genomic_alignment_samtools**](#sort_genomic_alignment_samtools)
  - BAM index file (`.bam.bai`); from
    [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
  - Transcript annotations file (12-column `.bed`); from
    [**extract_transcripts_as_bed12**](#extract_transcripts_as_bed12)
- **Parameters**
  - **rule_config.yaml**
    - `-c 0`: minimum number of read mapped to a transcript (default 10)
- **Output**
  - TIN score table (custom `tsv`); used in
    [**merge_TIN_scores**](#merge_tin_scores)
  

#### `salmon_quantmerge_genes`

Merge gene-level expression estimates for all samples with
[**Salmon**](#third-party-software-used).

> Rule is run once per sequencing mode

- **Input**
  - Gene expression tables (custom `.tsv`) for samples of same sequencing mode;
    from [**quantification_salmon**](#quantification_salmon)
- **Output**
  - Gene TPM table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)
  - Gene read count table (custom `.tsv`); used in
    [**pca_salmon**](#pca_salmon) and [**pca_kallisto**](#pca_kallisto)

#### `salmon_quantmerge_transcripts`

Merge transcript-level expression estimates for all samples with
[**Salmon**](#third-party-software-used).

> Rule is run once per sequencing mode

- **Input**
  - Transcript expression tables (custom `.tsv`) for samples of same sequencing
    mode; from [**quantification_salmon**](#quantification_salmon)
- **Output**
  - Transcript TPM table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)
  - Transcript read count table (custom `.tsv`); used in
    [**pca_salmon**](#pca_salmon) and [**pca_kallisto**](#pca_kallisto)

#### `kallisto_merge_genes`

Merge gene-level expression estimates for all samples with 
[custom script][custom-script-merge-kallisto].

> Rule is run once per sequencing mode

- **Input**
  - Transcript expression tables (custom `.h5`) for samples of same sequencing
    mode; from [**genome_quantification_kallisto**](#genome_quantification_kallisto) 
  - Gene annotation file (custom `.gtf`)
- **Output**
  - Gene TPM table (custom `.tsv`)
  - Gene read count table (custom `.tsv`)
  - Mapping gene/transcript IDs table (custom `.tsv`)
- **Non-configurable & non-default**
  - `-txOut FALSE`: gene-level summarization (default would be transcript level) 

#### `kallisto_merge_transcripts`

Merge transcript-level expression estimates for all samples with 
[custom script][custom-script-merge-kallisto].

> Rule is run once per sequencing mode

- **Input**
  - Transcript expression tables (custom `.h5`) for samples of same sequencing
    mode; from [**genome_quantification_kallisto**](#genome_quantification_kallisto) 
- **Output**
  - Transcript TPM table (custom `.tsv`)
  - Transcript read count table (custom `.tsv`)

#### `pca_kallisto`

Run PCA analysis on kallisto genes and transcripts with [custom script][custom-script-zpca].

> Rule is run one time for transcript estimates and one time for genes estimates

- **Input**
  - Transcript/Genes TPM table (custom `.tsv`)
- **Output**
  - Directory with PCA plots, scree plot and top loading scores.

#### `pca_salmon`

Run PCA analysis on salmon genes and transcripts with [custom script][custom-script-zpca].

> Rule is run one time for transcript estimates and one time for genes estimates

- **Input**
  - Transcript/Genes TPM table (custom `.tsv`)
- **Output**
  - Directory with PCA plots, scree plot and top loading scores.


#### `generate_alfa_index`

Create index for [**ALFA**](#third-party-software-used).

- **Input**
  - Gene annotation file (`.gtf`)
  - Chromosome length table `chrNameLength.txt`; from
    [**create_index_star**](#create_index_star)
- **Output**
  - ALFA index; stranded and unstranded; used in
    [**alfa_qc**](#alfa_qc)

#### `alfa_qc`

Annotate alignments with [**ALFA**](#third-party-software-used).

> For details on output plots, see [ALFA documentation][docs-alfa].   
> Note: the read orientation of a sample will be inferred from salmon `libtype` specified in `samples.tsv`

- **Input**
  - Coverage files, renamed (`.bg`); from
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
  - ALFA index, stranded; from [**generate_alfa_index**](#generate_alfa_index)
- **Parameters**

- **Output**
  - Figures for biotypes and feature categories (`.pdf`)
  - Feature counts table (custom `.tsv`); used in
    [**alfa_qc_all_samples**](#alfa_qc_all_samples)

#### `prepare_multiqc_config`

Prepare config file for [**MultiQC**](#third-party-software-used).

> Local rule

- **Input**
  - Directories created during
    [**prepare_files_for_report**](#prepare_files_for_report)
- **Parameters**
  All parameters for this rule have to be specified in main `config.yaml`
  - `--intro-text`
  - `--custom-logo`
  - `--url`
- **Output**
  - Config file (`.yaml`); used in [**multiqc_report**](#multiqc_report)

#### `multiqc_report`

Prepare interactive report from results and logs with
[**MultiQC**](#third-party-software-used).

- **Input**
  - Config file (`.yaml`); from
    [**prepare_multiqc_config**](#prepare_multiqc_config)
  - ALFA plot, combined (`.png`); from
    [**alfa_concat_results**](#alfa_concat_results)
  - Coverage file (`.bg`); from [**star_rpm**](#star_rpm)
  - FastQC output directory with report (`.txt`) and figures (`.png`); from
    [**fastqc**](#fastqc)
  - Gene TPM table (custom `.tsv`); from
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Gene read count table (custom `.tsv`); from
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Pseudoalignments file (`.sam`); from
    [**genome_quantification_kallisto**](#genome_quantification_kallisto)
  - TIN score box plots (`.pdf` and `.png`); from
    [**plot_TIN_scores**](#plot_tin_scores)
  - Transcript TPM table (custom `.tsv`); from
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)
  - Transcript read count table (custom `.tsv`); from
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)

- **Output**
  - Directory with automatically generated `.html` report

#### `finish`

Target rule as required by [Snakemake][docs-snakemake-target-rule].

> Local rule

- **Input**
  - MultiQC report; from [**multiqc_report**](#multiqc_report)
  - Coverage files, one per strand and sample (`.bw`); used in
    [**prepare_bigWig**](#prepare_bigwig)



### Sequencing mode-specific

> Steps described here have two variants, one with the specified names for
> samples prepared with a single-end sequencing protocol, one with `pe_`
> prepended to the specified name for samples prepared with a paired-end
> sequencing protocol.

#### `remove_adapters_cutadapt`

Remove adapter sequences from reads with
[**Cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from [**start**](#start)
- **Parameters**
  - **samples.tsv**
    - Adapters to be removed; specify in sample table columns `fq1_3p`, `fq1_5p`,
    `fq2_3p`, `fq2_5p`
  - **rule_config.yaml:**
    - `-m 10`: Discard processed reads that are shorter than 10 nt. If
    specified in `rule_config.yaml`, it will override ZARP's default value of
    `m=1` for this parameter. Note that this is different from `cutadapt`'s
    default behavior (`m=0`), which leads to empty reads being retained,
    causing problems in downstream applications in ZARP. We thus strongly
    recommend to **not** set the value of `m` to `0`! Refer to `cutadapt`'s
    [documentation][docs-cutadapt-m] for more information on the `m`
    parameter.
    - `-n 2`: search for all the given adapter sequences repeatedly, either until
    no adapter match was found or until 2 rounds have been performed. (default 1)

- **Output**
  - Reads file (`.fastq.gz`); used in
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)



#### `remove_polya_cutadapt`

Remove poly(A) tails from reads with
[**Cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_adapters_cutadapt**](#remove_adapters_cutadapt)
- **Parameters**
  - **samples.tsv**
    - Poly(A) stretches to be removed; specify in sample table columns `fq1_polya` and `fq2_polya`
  - **rule_config.yaml**
    - `-m 10`: Discard processed reads that are shorter than 10 nt. If
    specified in `rule_config.yaml`, it will override ZARP's default value of
    `m=1` for this parameter. Note that this is different from `cutadapt`'s
    default behavior (`m=0`), which leads to empty reads being retained,
    causing problems in downstream applications in ZARP. We thus strongly
    recommend to **not** set the value of `m` to `0`! Refer to `cutadapt`'s
    [documentation][docs-cutadapt-m] for more information on the `m`
    parameter.
    - `-O 1`: minimal overlap of 1 (default: 3)
- **Output**
  - Reads file (`.fastq.gz`); used in
    [**genome_quantification_kallisto**](#genome_quantification_kallisto),
    [**map_genome_star**](#map_genome_star) and
    [**quantification_salmon**](#quantification_salmon)



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

#### `quantification_salmon`

Estimate transcript- and gene-level expression with
[**Salmon**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Filtered annotation file (`.gtf`)
  - Index; from [**create_index_salmon**](#create_index_salmon)
- **Parameters**
  - **samples.tsv**
    - `libType`: see [Salmon manual][docs-salmon] for allowed values; specify in sample table column `libtype`
    - `--fldMean`: mean of distribution of fragment lengths; specify in sample table column `mean` **(single-end only)**
    - `--fldSD`: standard deviation of distribution of fragment lengths; specify in sample table column `sd` **(single-end only)**
  - **rule_config.yaml**
    - `--seqBias`: [correct for sequence specific
    biases](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias)
    - `--validateMappings`: enables selective alignment of the sequencing reads when mapping them to the transcriptome; this can improve both the sensitivity and specificity of mapping and, as a result, can [improve quantification accuracy](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings).
    - `--writeUnmappedNames`: write out the names of reads (or mates in paired-end reads) that do not map to the transcriptome. For paired-end this gives flags that indicate how a read failed to map **(paired-end only)**
- **Output**
  - Gene expression table (`quant.sf`); used in
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Transcript expression table ( `quant.sf`); used in
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)
  - `meta_info.json`
  - `flenDist.txt`
  

#### `genome_quantification_kallisto`

Generate pseudoalignments of reads to transcripts with
[**kallisto**](#third-party-software-used).
> Note: the kallisto strandedness parameter will be inferred from salmon `libtype` specified in `samples.tsv`

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Index; from [**create_index_kallisto**](#create_index_kallisto)
- **Parameters**
  - **samples.tsv**
    - `-l`: mean of distribution of fragment lengths; specify in sample table column `mean` **(single-end only)**
    - `-s`: standard deviation of distribution of fragment lengths; specify in sample table column `sd` **(single-end only)**
- **Output**
  - Pseudoalignments file (`.sam`) and
  - abundance (`.h5`) 
  used in [**kallisto_merge_genes**](#kallisto_merge_genes)
- **Non-configurable & non-default**
  - `--single`: Quantify single-end reads **(single-end only)**
  - `--pseudobam`: Save pseudoalignments to transcriptome to BAM file

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

[rule-graph]: images/flowchart.png
