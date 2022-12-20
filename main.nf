nextflow.enable.dsl=2

include { PHILOSOPHER } from './subworkflows/philosopher.nf'
include { ANNOTATE } from './subworkflows/annotate.nf'
include { RIBOTISH } from './subworkflows/ribotish.nf'
include { PREPARE } from './subworkflows/prepare_reads.nf'
include { GENOME } from './subworkflows/map_genome.nf'
include { rRNA } from './subworkflows/map_rrna.nf'
include { TRANSCRIPTOME } from './subworkflows/map_transcriptome.nf'
include { QC } from './subworkflows/qc.nf'

// check if files are present by converting params to channels
if ( params.run_mode == "full" ) {
    riboseq_reads_ch = channel.fromPath(params.riboseq_reads, checkIfExists: true)
    proteomics_reads_ch = channel.fromPath(params.proteomics_reads, checkIfExists: true)
    other_RNAs_sequence_ch = channel.fromPath(params.other_RNAs_sequence, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
}
if ( params.run_mode == "map to genome" || params.run_mode == "qc" ) {
    riboseq_reads_ch = channel.fromPath(params.riboseq_reads, checkIfExists: true)
    other_RNAs_sequence_ch = channel.fromPath(params.other_RNAs_sequence, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
}
if ( params.run_mode == "full" || params.run_mode == "fasta" ) {
    swissprot_ch = channel.fromPath(params.swissprot, checkIfExists: true)
}
if ( params.special_run_mode == "test" || params.run_mode == "proteomics" ) {
    predicted_peptides = channel.fromPath(params.test_database, checkIfExists: true)
}
if ( params.run_mode == "proteomics" ) {
    proteomics_reads_ch = channel.fromPath(params.proteomics_reads, checkIfExists: true)
}
if ( params.run_mode == "ribotish" ) {
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
    bam_sort_index_folder_ch = channel.fromPath(params.bam_sort_index_folder, checkIfExists: true, type: 'dir')
    swissprot_ch = channel.fromPath(params.swissprot, checkIfExists: true)
}
if ( params.run_mode == "fasta" ) {
    proteomics_reads_ch = channel.fromPath(params.proteomics_reads, checkIfExists: true)
    other_RNAs_sequence_ch = channel.fromPath(params.other_RNAs_sequence, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)                  
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
    prepared_out = channel.fromPath(params.prepared_fasta, checkIfExists: true)
}
if ( params.run_mode == "prepare" ) {
    riboseq_reads_ch = channel.fromPath(params.riboseq_reads, checkIfExists: true)
}

// main workflow that calls all processes in the subworkflows dir
workflow {

    if ( params.run_mode == "pull" ) {
            PULL_FILES()
    }

    if ( params.run_mode == "full" || params.run_mode == "prepare" || params.run_mode == "qc" || params.run_mode == "map to genome") {
        PREPARE(
            riboseq_reads_ch
        )
        prepared_out = PREPARE.out
    }

    if ( params.run_mode == "full" || params.run_mode == "qc" || params.run_mode == "map to genome" || params.run_mode == "fasta" ) {
        rRNA(
            genome_ch,
            other_RNAs_sequence_ch,
            prepared_out
        )
    }

    if ( params.run_mode == "full" || params.run_mode == "map to genome" || params.run_mode == "fasta" ) {
        GENOME(
            genome_ch,
            gtf_ch,
            rRNA.out.other_genes_unmapped_fasta
        )
        bam_sort_index_folder_ch = GENOME.out.bam_sort_index_folder
    }

	if ( params.run_mode == "full" || params.run_mode == "qc" || params.run_mode == "fasta" ) {
		ANNOTATE(
			gtf_ch,
			other_RNAs_sequence_ch,
			genome_ch
		)

		TRANSCRIPTOME(
			ANNOTATE.out.longest_pc_transcript_per_gene_fa,
			rRNA.out.other_genes_unmapped_fasta
		)

		QC(
			ANNOTATE.out.transcript_id_gene_id_CDS_tsv,
			TRANSCRIPTOME.out.transcripts_mapped_unique_sam,
			TRANSCRIPTOME.out.bam_bai_folder,
			rRNA.out.other_genes_mapped_sam
		)
	}

	if ( params.run_mode == "full" || params.run_mode == "ribotish" || params.run_mode == "fasta" ) {
        if ( params.special_run_mode != "test" ) {
            RIBOTISH(
                gtf_ch,
                bam_sort_index_folder_ch,
                genome_ch,
                swissprot_ch
            )
            predicted_peptides = RIBOTISH.out.speptide_combined
	    }
    }
    
    if ( params.run_mode == "full" || params.run_mode == "proteomics" || params.run_mode == "fasta" ) {
        PHILOSOPHER(
            predicted_peptides,
            proteomics_reads_ch
        )
    }
    
}
