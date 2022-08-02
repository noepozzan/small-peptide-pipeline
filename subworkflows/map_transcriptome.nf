#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SEGEMEHL_INDEX_TRANSCRIPTOME {

    label "segemehl"
	label 'indexing'

    publishDir "${params.map_dir}/segemehl_index_transcriptome", mode: 'copy', pattern: 'segemehl_index_transcriptome'
    publishDir "${params.log_dir}/segemehl_index_transcriptome", mode: 'copy', pattern: '*.log'

    input:
    path sequence

    output:
    path 'segemehl_index_transcriptome', emit: index
	path '*.log', emit: log

    script:
    """
    segemehl.x \
    	-x segemehl_index_transcriptome \
    	-d ${sequence} \
    	&> segemehl_index_transcriptome.log

    """

}

process MAP_TRANSCRIPTOME_SEGEMEHL {
    
    label "segemehl"
	label 'mapping'

    publishDir "${params.map_dir}/map_transcriptome_segemehl", mode: 'copy', pattern: "*.transcripts_mapped_sam"
	publishDir "${params.map_dir}/map_transcriptome_segemehl", mode: 'copy', pattern: "*.transcripts_unmapped_fa"
    publishDir "${params.log_dir}/map_transcriptome_segemehl", mode: 'copy', pattern: '*.log'

    input:
    each(path(reads))
    path index
    path sequence

    output:
    path '*.transcripts_mapped_sam', emit: mapped
    path '*.transcripts_unmapped_fa', emit: unmapped
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    segemehl.x \
		-i ${index} \
		-d ${sequence} \
		-q ${reads} \
		${params.segemehl_args} \
		-o \${prefix}.transcripts_mapped_sam \
		-u \${prefix}.transcripts_unmapped_fa \
		&> \${prefix}_map_to_transcripts.log
	
	"""

}

process REMOVE_MULTIMAPPERS {
    
    label "pysam"

    publishDir "${params.map_dir}/remove_multimappers", mode: 'copy', pattern: '*.transcripts_mapped_unique_sam'
    publishDir "${params.log_dir}/remove_multimappers", mode: 'copy', pattern: '*.log'

    input:
    path transcripts_mapped_sam

    output:
    path '*.transcripts_mapped_unique_sam', emit: unique
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${transcripts_mapped_sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    (grep -P \"^@|\tNH:i:1\t\" \
		${transcripts_mapped_sam} \
		> \${prefix}.transcripts_mapped_unique_sam) \
		&> \${prefix}_remove_multimappers.log

    """

}

process SAM_TO_BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.map_dir}/sam_to_bam_sort_and_index", mode: 'copy', pattern: "*.sorted_indexed.bam"
	publishDir "${params.map_dir}/sam_to_bam_sort_and_index", mode: 'copy', pattern: "*.sorted_indexed.bam.bai"
	publishDir "${params.map_dir}/sam_to_bam_sort_and_index", mode: 'copy', pattern: "*.folder_sorted_indexed_bam"
    publishDir "${params.log_dir}/sam_to_bam_sort_and_index", mode: 'copy', pattern: '*.log'

    input:
    path sam

    output:
    path '*.sorted_indexed.bam', emit: bam
    path '*.sorted_indexed.bam.bai', emit: bai
    path '*.folder_sorted_indexed_bam', emit: folder
	path '*.log', emit: log

    script:
    """
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    samtools view -bS ${sam} \
		| samtools sort - > \${prefix}.sorted_indexed.bam; \
		samtools index \${prefix}.sorted_indexed.bam; \
		&> \${prefix}_bam_sorted_and_indexed.log

    mkdir \${prefix}.folder_sorted_indexed_bam
    cp \${prefix}.sorted* \${prefix}.folder_sorted_indexed_bam

    """

}


workflow TRANSCRIPTOME_PIPE {

    take:
    longest_pc_transcript_per_gene_fa  
    other_genes_unmapped_fasta

    main:   
	// generate index for later mapping
	SEGEMEHL_INDEX_TRANSCRIPTOME(
		longest_pc_transcript_per_gene_fa
	)
	transcriptome_index = SEGEMEHL_INDEX_TRANSCRIPTOME.out.index

	// map the filtered reads to the transcriptome to then do QC
	MAP_TRANSCRIPTOME_SEGEMEHL(
		other_genes_unmapped_fasta,
		transcriptome_index,
		longest_pc_transcript_per_gene_fa
	)
	transcripts_mapped_sam = MAP_TRANSCRIPTOME_SEGEMEHL.out.mapped
	transcripts_unmapped_fasta = MAP_TRANSCRIPTOME_SEGEMEHL.out.unmapped

	REMOVE_MULTIMAPPERS(
		transcripts_mapped_sam
	)
	transcripts_mapped_unique_sam = REMOVE_MULTIMAPPERS.out.unique
	
	SAM_TO_BAM_SORT_AND_INDEX(
		transcripts_mapped_unique_sam
	)
	transcripts_mapped_unique_sorted_bam = SAM_TO_BAM_SORT_AND_INDEX.out.bam
	transcripts_mapped_unique_sorted_bam_bai = SAM_TO_BAM_SORT_AND_INDEX.out.bai
	bam_bai_folder = SAM_TO_BAM_SORT_AND_INDEX.out.folder

	emit:
	transcripts_mapped_unique_sam
	bam_bai_folder
}


