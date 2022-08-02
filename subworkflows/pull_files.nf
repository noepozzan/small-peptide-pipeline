nextflow.enable.dsl=2

process PULL {

	label "rcrunch_python"

	publishDir "${params.pull_dir}/pull", mode: 'copy', pattern: "*.fa"
    publishDir "${params.pull_dir}/pull", mode: 'copy', pattern: "*.gtf"
    publishDir "${params.pull_dir}/pull", mode: 'copy', pattern: "*.gff3"
	publishDir "${params.log_dir}/pull_files/pull", mode: 'copy', pattern: '*.log'

	input:
	path dir

	output:
	path '*.fa', emit: genome_fasta
	path '*.gtf', emit: gtf
	path '*.gff3', emit: gff3
	path '*.log', emit: log

	script:
	"""
	# better idea:
	# because the versions change all the time, just fix 1 version from the archives

	# whole genome DNA fasta file
	wget http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
		&> genome_fasta.log
	# gff3 annotation files of noncoding RNA and mRNA
	wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/20.0/genome_coordinates/gff3/mus_musculus.GRCm39.gff3.gz \
		&> ncRNA_annotation.log
	wget http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.chr.gtf.gz \
		&> mRNA_annotation.log

	# decompress files
	zcat Mus_musculus.GRCm39.dna.primary_assembly.fa.gz > Mus_musculus.GRCm39.dna.primary_assembly.fa
	zcat mus_musculus.GRCm39.gff3.gz > mus_musculus.GRCm39.gff3
	zcat Mus_musculus.GRCm39.105.chr.gtf.gz > Mus_musculus.GRCm39.105.chr.gtf

	# cp *.fa *.gtf *.gff3 ${dir}
	"""


}

process FASTA_INDEX {

	label "samtools"

	publishDir "${params.pull_dir}/fasta_index", mode: 'copy', pattern: '*.fai'
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

process GFFREAD {

	label "gffread"

	publishDir "${params.pull_dir}/gffread", mode: 'copy', pattern: '*.gtf'
	publishDir "${params.log_dir}/gffread", mode: 'copy', pattern: '*.log'

	input:
	path gff3

	output:
	path "*.gtf", emit: gtf
	path "*.log", emit: log

	script:
	"""
	input=\$(basename ${gff3})
    prefix=\$(echo \$input | cut -d '.' -f 1,2)

	gffread ${gff3} \
		-T \
		-o \${prefix}.gtf \
		&> gff3_to_gtf.log

	"""
}

process DUPLICATES {

    label "rduplicate"

    publishDir "${params.pull_dir}/duplicates", mode: 'copy', pattern: '*_out.gtf'
    publishDir "${params.log_dir}/duplicates", mode: 'copy', pattern: '*.log'

    input:
    path gtf
    path duplicates_r_script

    output:
    path "*_out.gtf", emit: deduplicated_gtf
    path "*.log", emit: log

    script:
    """
    input=\$(basename ${gtf})
    prefix=\${input%.*}

    Rscript ${duplicates_r_script} \
        --gtf ${gtf} \
        --out \${prefix}_out.gtf \
        --verbose \
        &> \${prefix}_duplicates.log

    """
    
}

process CONCAT {

	label "gffread"

	publishDir "${params.pull_dir}/concat", mode: 'copy', pattern: 'combined.gtf'
    publishDir "${params.log_dir}/concat", mode: 'copy', pattern: '*.log' 

	input:
	//path gene_gtf
	//path noncoding_gtf
    path gtfs

	output:
	path "combined.gtf", emit: gtf
	path "*.log", emit: log, optional: true

	script:
	"""
	cat ${gtfs} >> tmp_combined.gtf
    grep -Ev '^#' tmp_combined.gtf > combined.gtf
	"""
}

process MOVE {

	label "philosopher"

	publishDir "${params.pull_dir}/move", mode: 'copy', pattern: '*'

	input:
	path genome_fasta
	path genome_fai
	path gene_gtf
	path noncoding_gtf
	path combined_gtf

	script:
	"""
	cp ${genome_fasta} ${params.genome}
	cp ${genome_fai} ${params.genome_fai}
	cp ${gene_gtf} ${params.gtf}
	cp ${noncoding_gtf} ${params.rnacentral_gtf}
	cp ${combined_gtf} ${params.combined_gtf}
	"""

}

workflow PULL_FILES {

	take:
	directory_path

	main:
	// this subworkflow pulls gtf, gff3 and fasta files from the web
	// and converts gff3 files to gtf, then concatenates the gtf files
	// into 1 larger gtf file which will be used to map reads to.
	// the workflow also puts the files in the right dir
	// and writes this to the nextflow.config file

	PULL(
		directory_path
	)

	FASTA_INDEX(
		PULL.out.genome_fasta
	)

	GFFREAD(
		PULL.out.gff3
	)

    PULL.out.gtf
        .mix(GFFREAD.out.gtf)
        .set{gtfs_ch}

    DUPLICATES(
        gtfs_ch,
        params.duplicate_r_script
    )
	
    CONCAT(
        DUPLICATES.out.deduplicated_gtf.collect()
	)
	combined_gtf = CONCAT.out.gtf

	MOVE(
		PULL.out.genome_fasta,
        FASTA_INDEX.out.fai,
        PULL.out.gtf,
        GFFREAD.out.gtf,
		combined_gtf
	)

	emit:
	combined_gtf

}




