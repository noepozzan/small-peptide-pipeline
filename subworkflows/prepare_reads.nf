#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRIM_FIRST_BASES {

    label "cutadapt"

    publishDir "${params.reads_dir}/trim_first_bases", mode: 'copy', pattern: '*.trimmed_first_bases'
    publishDir "${params.log_dir}/trim_first_bases", mode: 'copy', pattern: '*.log'

    input:
    path reads

    output:
    path '*.trimmed_first_bases', emit: trimmed_first_bases
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    (cutadapt \
		--cut ${params.cut} \
		--minimum-length ${params.minimum_length} \
		${reads} | gzip > \
		\${prefix}.trimmed_first_bases) \
		&> \${prefix}_trim_first_bases.log

    """

}

process CLIP_READS {
    
    label "fastx"

    publishDir "${params.reads_dir}/clip_reads", mode: 'copy', pattern: '*.pro_clipped'
    publishDir "${params.log_dir}/clip_reads", mode: 'copy', pattern: '*.log'

    input:
    path reads

    output:
    path '*.pro_clipped', emit: pro_clipped
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastx_clipper \
		${params.clip_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_clipped \
		&> \${prefix}_clip_reads.log
    
    """

}


process TRIM_READS {
    
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 5

    label "fastx"

    publishDir "${params.reads_dir}/trim_reads", mode: 'copy', pattern: '*.pro_trimmed'
    publishDir "${params.log_dir}/trim_reads", mode: 'copy', pattern: '*.log'

    input:
    path reads

    output:
    path '*.pro_trimmed', emit: pro_trimmed
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_trimmer \
		${params.trim_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_trimmed \
		&> \${prefix}_trim_reads.log

    """

}

process FILTER_READS {
    
    label "fastx"

    publishDir "${params.reads_dir}/filter_reads", mode: 'copy', pattern: '*.pro_filtered'
    publishDir "${params.log_dir}/filter_reads", mode: 'copy', pattern: '*.log'

    input:
    path reads

    output:
    path '*.pro_filtered', emit: pro_filtered
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_filter \
		${params.filter_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_filtered \
		&> \${prefix}_filter_reads.log

    """

}

process FASTQ_TO_FASTA {
    
    label "fastx"

    publishDir "${params.reads_dir}/fastq_to_fasta", mode: 'copy', pattern: '*.pro_filtered_fasta'
    publishDir "${params.log_dir}/fastq_to_fasta", mode: 'copy', pattern: '*.log'

    input:
    path reads

    output:
    path '*.pro_filtered_fasta', emit: fasta
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_to_fasta \
		${params.fastq_to_fasta_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_filtered_fasta \
		&> \${prefix}_fastq_to_fasta.log

    """

}


workflow READS_PIPE {

    take:
    riboseq_reads_ch

    main:   
	// prepare riboseq reads
	TRIM_FIRST_BASES(
		riboseq_reads_ch
	)
	
	CLIP_READS(
		TRIM_FIRST_BASES.out.trimmed_first_bases
	)
	
	TRIM_READS(
		CLIP_READS.out.pro_clipped
	)
	
	FILTER_READS(
		TRIM_READS.out.pro_trimmed
	)
	
	FASTQ_TO_FASTA(
		FILTER_READS.out.pro_filtered
	)

    emit:
	FASTQ_TO_FASTA.out.fasta

}


