nextflow.enable.dsl=2

include { PHILOSOPHER } from './subworkflows/philosopher.nf'
include { PHILOSOPHER_PARALLEL } from './subworkflows/philosopher_parallel.nf'
include { ANNOTATE } from './subworkflows/annotate.nf'
include { RIBOTISH } from './subworkflows/ribotish.nf'
include { READS } from './subworkflows/prepare_reads.nf'
include { GENOME } from './subworkflows/map_genome.nf'
include { rRNA } from './subworkflows/map_rrna.nf'
include { TRANSCRIPTOME } from './subworkflows/map_transcriptome.nf'
include { QC } from './subworkflows/qc.nf'
include { FIX_NAMES } from './subworkflows/fix_names.nf'
include { PULL_FILES } from './subworkflows/pull_files.nf'
include { RAWMZML } from './subworkflows/raw_to_mzml.nf'

// check if files are present by converting params to channels
if ( params.run_mode == "test" || params.run_mode == "full" ) {
    riboseq_reads_ch = channel.fromPath(params.riboseq_reads, checkIfExists: true)
    proteomics_reads_ch = channel.fromPath(params.proteomics_reads, checkIfExists: true)
    other_RNAs_sequence_ch = channel.fromPath(params.other_RNAs_sequence, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
    genome_fai_ch = channel.fromPath(params.genome_fai, checkIfExists: true)
}
if ( params.run_mode == "ribotish" ) {
    bam_sort_index_folder_ch = channel.fromPath(params.bam_sort_index_folder, checkIfExists: true)
}
if ( params.run_mode == "test" || params.run_mode == "proteomics" ) {
    predicted_peptides = channel.fromPath(params.test_database, checkIfExists: true)
}
if ( params.run_mode == "proteomics" ) {
    proteomics_reads_ch = channel.fromPath(params.proteomics_reads, checkIfExists: true)
}
if ( params.run_mode == "ribotish" ) {
    genome_fai_ch = channel.fromPath(params.genome_fai, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch = channel.fromPath(params.genome, checkIfExists: true)
    // this input is very special: It is a folder containing a sorted+indexed bam file and its index
    bam_sort_index_folder_ch = channel.fromPath(params.bam_sort_index_folder, checkIfExists: true)
}
if ( params.run_mode == "pull" ) {
    directory_path = channel.fromPath(params.directory_path, checkIfExists: true)
    PULL_FILES(
        directory_path
    )
}


// this should be replaced by Meric's image
workflow RAW_TO_MZML {

	RAWMZML()

}

// kind of unncecessary?
workflow MAP_NAMES {

	csv_file = channel.fromPath("${projectDir}/data/experimental_conditions.csv")
	fix_script = channel.fromPath("${projectDir}/data/python_scripts/fix_names.py")

	FIX_NAMES(
		csv_file,
		fix_script
	)

}

// main workflow that calls all processes in the subworkflows dir
workflow {

    if ( params.run_mode == "full" || params.run_mode == "test" || params.run_mode == "prepare" || params.run_mode == "qc" || params.run_mode == "map to genome") {
        READS(
            riboseq_reads_ch
        )
    }

    if ( params.run_mode == "full" || params.run_mode == "test" || params.run_mode == "qc" || params.run_mode == "map to genome") {
        rRNA(
            genome_ch,
            other_RNAs_sequence_ch,
            READS.out
        )
    }

    if ( params.run_mode == "full" || params.run_mode == "test" || params.run_mode == "map to genome" ) {
        GENOME(
            genome_ch,
            gtf_ch,
            rRNA.out.other_genes_unmapped_fasta
        )
        bam_sort_index_folder_ch = GENOME.out.bam_sort_index_folder
    }

	if ( params.run_mode == "full" || params.run_mode == "test" || params.run_mode == "qc" ) {
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

	if ( params.run_mode == "full" || params.run_mode == "ribotish" ) {
		RIBOTISH(
			gtf_ch,
			bam_sort_index_folder_ch,
			genome_ch,
			genome_fai_ch,
		)
		predicted_peptides = RIBOTISH.out.speptide_combined
	}
    
    if ( params.run_mode == "full" || params.run_mode == "test" || params.run_mode == "proteomics") {
        PHILOSOPHER(
            predicted_peptides,
            proteomics_reads_ch
        )
    }
    
/*
	PHILOSOPHER_PARALLEL(
		//PHILOSOPHER.out.ionquant,
		PHILOSOPHER.out.report,
		PHILOSOPHER.out.msfragger_params,
		predicted_speptides,
        proteomics_reads_ch		
	)
*/

}
