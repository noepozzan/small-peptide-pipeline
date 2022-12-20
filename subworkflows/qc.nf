#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process COUNT_OVERREPRESENTED_SEQUENCES_OTHER {
    
    label "pysam"

    publishDir "${params.qc_dir}/count_overrepresented_sequences_other", mode: 'copy', pattern: '*.overrepresented_sequences_counts'
    publishDir "${params.log_dir}/count_overrepresented_sequences_other", mode: 'copy', pattern: '*.log'

    input:
    each(path(sam))
    path script_py

    output:
    path '*.overrepresented_sequences_counts', emit: sequences
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    python ${script_py} \
		--sam ${sam} \
		--out \${prefix}.overrepresented_sequences_counts \
		&> \${prefix}_count_overrepresented_sequences_other.log
    
    : '
    grep -P -v \"^@\" ${sam} \
	| cut -f 10 | sort | uniq -c \
	| sort -n -r > \
	\${prefix}.overrepresented_sequences_counts \
	2> \${prefix}_count_overrepresented_sequences_other.log
    '

    """

}

process READ_LENGTH_HISTOGRAM {
    
    label "rcrunch_python"

    publishDir "${params.qc_dir}/read_length_histogram", mode: 'copy', pattern: '*.dir'
    publishDir "${params.log_dir}/read_length_histogram", mode: 'copy', pattern: '*.log'

    input:
    each(path(sam))
    path script_py

    output:
    path '*.dir', emit: hist
	path '*.log', emit: log

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    python ${script_py} \
		--sam ${sam} \
		--outdir \${prefix}.dir \
		&> \${prefix}_read_length_histogram.log

    """

}

process DETERMINE_P_SITE_OFFSET {
    
    label "pysam"

    publishDir "${params.qc_dir}/determine_p_site_offset", mode: 'copy', pattern: '*.alignment_offset.json'
    publishDir "${params.log_dir}/determine_p_site_offset", mode: 'copy', pattern: '*.log'

    input:
    each(path(bam_folder))
    path transcript_id_gene_id_CDS
    path script_py

    output:
    path '*.alignment_offset.json', emit: offsets
    path '*.log', emit: log

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${bam_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    echo \$workd
    echo \${prefix}

    cp ${bam_folder}/* .

    python ${script_py} \
		--bam \${prefix}.*.bam \
		--cds_coordinates ${transcript_id_gene_id_CDS} \
		--outdir \${prefix}_p_site_offset \
		&> \${prefix}_p_site_offset.log

    cp \${prefix}_p_site_offset/* .
    mv alignment_offset.json \${prefix}.alignment_offset.json

    """

}

process COUNT_READS {

    label "pysam"

    publishDir "${params.qc_dir}/count_reads", mode: 'copy', pattern: '*.counts.tsv'
    publishDir "${params.log_dir}/count_reads", mode: 'copy', pattern: '*.log'

    input:
    each(path(bam_folder_offsets))
    path transcript_id_gene_id_CDS
    path script_py

    output:
    path '*.counts.tsv', emit: counts, optional: true
	path '*.log', emit: log

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${bam_folder_offsets[0]})
    prefix=\$(echo \$input | cut -d '.' -f 1)

	cp ${bam_folder_offsets[0]}/* .

    mkdir \${prefix}.count_reads
	
	python ${script_py} \
    	--bam *.bam \
        --tsv ${transcript_id_gene_id_CDS} \
        --json ${bam_folder_offsets[1]} \
        --outdir \${prefix}.count_reads \
		&> \${prefix}.count_reads.log

    cp \${prefix}.count_reads/* .
	mv counts.tsv \${prefix}.counts.tsv

    """

}

process CHECK_PERIODICITY {

    label "rcrunch_python"

    publishDir "${params.qc_dir}/check_periodicity", mode: 'copy', pattern: '*.periodicity_start.pdf'
	publishDir "${params.qc_dir}/check_periodicity", mode: 'copy', pattern: '*.periodicity_stop.pdf'
	publishDir "${params.qc_dir}/check_periodicity", mode: 'copy', pattern: '*.Periodicity_Analysis_Start_Ribo_Seq.txt'
    publishDir "${params.log_dir}/check_periodicity", mode: 'copy', pattern: '*.log'

    input:
    each(path(bam_folder_offsets))
    path transcript_id_gene_id_CDS
    path script_py

    output:
    path '*.periodicity_start.pdf', emit: start, optional: true
    path '*.periodicity_stop.pdf', emit: stop, optional: true
    path '*.Periodicity_Analysis_Start_Ribo_Seq.txt', emit: txt, optional: true

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${bam_folder_offsets[0]})
    prefix=\$(echo \$input | cut -d '.' -f 1)

	cp ${bam_folder_offsets[0]}/* .

	mkdir \${prefix}.check_periodicity

	python ${script_py} \
    	--bam *.bam \
      	--tsv ${transcript_id_gene_id_CDS} \
       	--json ${bam_folder_offsets[1]} \
       	--outdir \${prefix}.check_periodicity \
       	--codnum ${params.check_periodicity_codnum} \
        &> \${prefix}_check_periodicity.log

    	cp \${prefix}.check_periodicity/* .
    	mv periodicity_start.pdf \${prefix}.periodicity_start.pdf
		mv periodicity_stop.pdf \${prefix}.periodicity_stop.pdf
		mv Periodicity_Analysis_Start_Ribo_Seq.txt \${prefix}.Periodicity_Analysis_Start_Ribo_Seq.txt

    """

}

process MULTIQC {

    label "multiqc"

    publishDir "${params.qc_dir}/multiqc", mode: 'copy', pattern: '*.html'
    publishDir "${params.log_dir}/multiqc", mode: 'copy', pattern: '*.html'

    input:
    path upstream_process_pseudo

    output:
    path '*.html', emit: multiqc, optional: true

    script:
    """
    workd=\$(pwd)

    cd ${projectDir}
    multiqc . \
        -n SMAPP_report.html

    mv SMAPP_report_data SMAPP_report.html \$workd
    """

}

workflow QC {

    take:
	transcript_id_gene_id_CDS
	transcripts_mapped_unique_sam
	bam_bai_folder
	other_genes_mapped_sam

    main:   
    COUNT_OVERREPRESENTED_SEQUENCES_OTHER(
        other_genes_mapped_sam,
        params.find_overrepresented_sequences_script
    )

	READ_LENGTH_HISTOGRAM(
		transcripts_mapped_unique_sam,
		params.plot_read_lengths_script
	)
      
	DETERMINE_P_SITE_OFFSET(
		bam_bai_folder,
		transcript_id_gene_id_CDS,
		params.determine_p_site_offsets_script
	)
    alignment_offset_json = DETERMINE_P_SITE_OFFSET.out.offsets

	bam_bai_folder
        .combine(alignment_offset_json)
        .filter{ it[0].baseName.split("\\.")[0] == it[1].baseName.split("\\.")[0] }
        .set{ bam_folder_offsets }

	COUNT_READS(
		bam_folder_offsets,
		transcript_id_gene_id_CDS,
		params.count_reads_script
	)
      
	CHECK_PERIODICITY(
		bam_folder_offsets,
		transcript_id_gene_id_CDS,
		params.check_periodicity_script
	)

    MULTIQC(
        CHECK_PERIODICITY.out.start.collect()
    )
}


