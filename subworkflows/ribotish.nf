nextflow.enable.dsl=2

// https://github.com/zhpn1024/ribotish

process FASTA_INDEX {

	label "samtools"

	publishDir "${params.ribotish_dir}/fasta_index", mode: 'copy', pattern: '*.fai'
	publishDir "${params.log_dir}/fasta_index", mode: 'copy', pattern: '*.log'

	input:
	path genome_fasta

	output:
	path "*.fai", emit: fai
	path "*.log", emit: log

	script:
	"""
	samtools faidx \
		${genome_fasta} \
		&> fasta_index.log

	"""

}

process RIBOTISH_QUALITY {

    label "ribotish"

    publishDir "${params.ribotish_dir}/ribotish_quality", mode: 'copy', pattern: "*.para.py"
	publishDir "${params.ribotish_dir}/ribotish_quality", mode: 'copy', pattern: "*.pdf"
	publishDir "${params.ribotish_dir}/ribotish_quality", mode: 'copy', pattern: "*.txt"
    publishDir "${params.log_dir}/ribotish_quality", mode: 'copy', pattern: "*.log"

    input:
    each(path(bam_sort_index_folder))
    path gtf_file

    output:
	path '*.para.py', emit: offset
    path '*.pdf', emit: ribo_pdf
	path '*.txt', emit: ribo_txt
	path '*.log', emit: log

    script:
    if( params.riboseq_mode == "regular" )
	"""
    input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

    ribotish quality \
		-b \${prefix}.*.bam \
		-g ${gtf_file} \
		--th ${params.ribotish_quality_th} \
		&> \${prefix}_ribotish_quality.log

    """
	else if( params.riboseq_mode == "TI" )
	"""
	input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

	ribotish quality -t \
        -b \${prefix}.*.bam \
        -g ${gtf_file} \
        --th ${params.ribotish_quality_th} \
		&> \${prefix}_ribotish_quality.log

	"""

}

process RIBOTISH_PREDICT {

	echo true

    label "ribotish"
	label "predicting"

    publishDir "${params.ribotish_dir}/ribotish_predict", mode: 'copy', pattern: "*ribotish_pred_all.txt"
    publishDir "${params.log_dir}/ribotish_predict", mode: 'copy', pattern: "*.log"

    input:
	each(path(bam_folder_offsets))
	path gtf
	path genome
 
    output:
    path '*.ribotish_pred_all.txt', emit: ribo_pred
	path '*.log', emit: log

    script:
	if( params.riboseq_mode == "regular" )
    """
	input=\$(basename ${bam_folder_offsets[0]})
    prefix=\$(echo \$input | cut -d '.' -f 1)

	cp ${bam_folder_offsets[0]}/* .

	ribotish predict \
		-b *.bam \
		-g ${gtf} \
		-f ${genome} \
		--ribopara ${bam_folder_offsets[1]} \
		${params.ribotish_predict_mode} \
		-o \${prefix}.ribotish_pred.txt \
		&> \${prefix}_ribotish_predict.log

	"""
	else if( params.riboseq_mode == "TI" )
	"""
    input=\$(basename ${bam_folder_offsets[0]})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_folder_offsets[0]}/* .

    ribotish predict \
        -t *.bam \
        -g ${gtf} \
        -f ${genome} \
        --ribopara ${bam_folder_offsets[1]} \
        ${params.ribotish_predict_mode} \
        -o \${prefix}.ribotish_pred.txt \
        &> \${prefix}_ribotish_predict.log

    """

}

process GFFREAD {

	label "gffread"

	publishDir "${params.ribotish_dir}/gffread", mode: 'copy', pattern: 'transcripts.fa'
	publishDir "${params.log_dir}/gffread", mode: 'copy', pattern: '*.log'

	input:
	path genome
	path genome_fai
	path gtf

	output:
	path 'transcripts.fa', emit: fasta	
	path '*.log', emit: log

	script:
	"""

	gffread \
        -w transcripts.fa \
        -g ${genome} \
        ${gtf} \
		&> gffread.log

	"""

}

process SORF_TO_PEPTIDE {

	label "sorf_to_speptide"

    publishDir "${params.ribotish_dir}/sorf_to_peptide", mode: 'copy', pattern: '*.speptide'
    publishDir "${params.log_dir}/sorf_to_peptide", mode: 'copy', pattern: '*.log'

    input:
    each(path(ribo_pred))
	path fasta
	path python_script

    output:
	path '*.speptide', emit: prediction
	path '*.log', emit: log

    script:
    """
	input=\$(basename ${ribo_pred})
    prefix=\$(echo \$input | cut -d '.' -f 1)

	python3 ${python_script} \
		--ribo_pred ${ribo_pred} \
		--fasta ${fasta} \
		--out \${prefix}.speptide \
		&> \${prefix}_sorf_to_peptide.log
	
	"""

}

process COMBINE {
   
	label "sorf_to_speptide"

    publishDir "${params.ribotish_dir}/combine", mode: 'copy', pattern: 'combined_speptide*'

    input:
    path ribo_pred

    output:
	path 'combined_speptide*', emit: combined_prediction

    script:
    """
    : '
	WRITE_HEADER="true"
    for VAR in ${ribo_pred}
    do
		if [ "\$WRITE_HEADER" == "true" ]; then
	    	head -1 \$VAR > combined_speptide.fasta
	    	WRITE_HEADER="false"
		fi
		tail -n +2 -q \$VAR >> combined_speptide.fasta
    done
	'
    for VAR in ${ribo_pred}
    do
        cat \$VAR >> combined_speptide.fasta
    done
	"""

}


workflow RIBOTISH {

    take:
    gtf_ch
    bam_sort_index_folder_ch
	genome_ch
	
    main:

    FASTA_INDEX(
		genome_ch
	)
	
    RIBOTISH_QUALITY(
		bam_sort_index_folder_ch,
		gtf_ch
	)
	offsets = RIBOTISH_QUALITY.out.offset

	bam_sort_index_folder_ch
		.combine(offsets)
		.filter{ it[0].baseName.split("\\.")[0] == it[1].baseName.split("\\.")[0] }
		.set{ bam_folder_offsets }

	RIBOTISH_PREDICT(
        bam_folder_offsets,
        gtf_ch,
        genome_ch
    )

	GFFREAD(
		genome_ch,
        FASTA_INDEX.out.fai,
		gtf_ch
	)
	transcripts_fa = GFFREAD.out.fasta

	SORF_TO_PEPTIDE(
		RIBOTISH_PREDICT.out.ribo_pred,
		transcripts_fa,
		params.sorf_peptide_script
	)
	speptide = SORF_TO_PEPTIDE.out.prediction

	COMBINE(
		speptide.collect()
	)
	speptide_combined = COMBINE.out.combined_prediction

	emit:
	speptide_combined

}


