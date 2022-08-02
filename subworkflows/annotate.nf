#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SELECT_LONGEST_CODING_TRANSCRIPT {
    
    label "htseq"

    publishDir "${params.annotate_dir}/select_longest_coding_transcript", mode: 'copy', pattern: 'longest_coding_transcript_per_gene.gtf'
    publishDir "${params.log_dir}/select_longest_coding_transcript", mode: 'copy', pattern: '*.log'

    input:
    path input_gtf
    path select_longest_ct_py

    output:
    path 'longest_coding_transcript_per_gene.gtf', emit: gtf
	path '*.log', emit: log

    script:
    """
    python ${select_longest_ct_py} \
		--gtf ${input_gtf} \
		--out longest_coding_transcript_per_gene.gtf \
		&> select_longest_coding_transcript.log

    """

}


process EXTRACT_TRANSCRIPT_SEQUENCES {
    
    label "cufflinks"

    publishDir "${params.annotate_dir}/extract_transcript_sequences", mode: 'copy', pattern: 'transcripts_sequences.out'
    publishDir "${params.log_dir}/extract_transcript_sequences", mode: 'copy', pattern: '*.log'

    input:
    path gtf
    path genome

    output:
    path 'transcripts_sequences.out', emit: fasta
	path '*.log', emit: log

    script:
    """
    gffread \
		${gtf} \
		-g ${genome} \
		-w transcripts_sequences.out \
		&> extract_transcript_sequences.log

    """

}

process CREATE_TAB_DELIMITED_CDS_FILE {
    
    echo true

    label "htseq_biopython"

    publishDir "${params.annotate_dir}/create_tab_delimited_CDS_file", mode: 'copy', pattern: 'CDS.tsv'
    publishDir "${params.log_dir}/create_tab_delimited_CDS_file", mode: 'copy', pattern: '*.log'

    input:
    path gtf
    path transcripts
    path td_CDS_script_py

    output:
    path 'CDS.tsv', emit: tsv
	path '*.log', emit: log

    script:
    """
    python ${td_CDS_script_py} \
		--gtf ${gtf} \
		--fasta ${transcripts} \
		--out CDS.tsv \
		&> create_tab_delimited_CDS_file.log

    """

}

process CREATE_BED_CDS_FILE {
    
    label "htseq_biopython"

    publishDir "${params.annotate_dir}/create_bed_cds_file", mode: 'copy', pattern: 'CDS.bed'
    publishDir "${params.log_dir}/create_bed_cds_file", mode: 'copy', pattern: '*.log'

    input:
    path tsv

    output:
    path 'CDS.bed', emit: bed
	path '*.log', emit: log

    script:
    """
    tail -n+2 ${tsv} \
		| awk '{print \$1 "\t" \$3-1 "\t" \$4 "\t" \$2 }' > CDS.bed \
		&> create_bed_cds_file.log

    """

}


workflow ANNOTATE_PIPE {

    take:
    gtf_ch
    other_RNAs_sequence_ch
    genome_ch

    main:
	SELECT_LONGEST_CODING_TRANSCRIPT(
		gtf_ch,
		params.lct_script
	)
	longest_pc_transcript_per_gene_gtf = SELECT_LONGEST_CODING_TRANSCRIPT.out.gtf
	
	EXTRACT_TRANSCRIPT_SEQUENCES(
		longest_pc_transcript_per_gene_gtf,
		genome_ch
	)
	longest_pc_transcript_per_gene_fa = EXTRACT_TRANSCRIPT_SEQUENCES.out.fasta

    CREATE_TAB_DELIMITED_CDS_FILE(
		longest_pc_transcript_per_gene_gtf,
		longest_pc_transcript_per_gene_fa,
		params.ctdCDS_script
	)
	transcript_id_gene_id_CDS_tsv = CREATE_TAB_DELIMITED_CDS_FILE.out.tsv

	CREATE_BED_CDS_FILE(
		transcript_id_gene_id_CDS_tsv
	)
    transcript_id_gene_id_CDS_bed = CREATE_BED_CDS_FILE.out.bed

    emit:
    longest_pc_transcript_per_gene_gtf
    longest_pc_transcript_per_gene_fa
    transcript_id_gene_id_CDS_tsv
    transcript_id_gene_id_CDS_bed

}













