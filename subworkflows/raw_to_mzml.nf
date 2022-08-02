nextflow.enable.dsl=2

process RAW_TO_MZML {

	echo true
	//label ".."

	publishDir "${projectDir}/results/raw_to_mzML", mode: 'copy', pattern: '*.mzML'

	input:
	each(path(bash_script))
	path raw_files

	output:
	path "*.mzML"

	script:
	"""
	./${bash_script} ${raw_files}
	
	mkdir ${projectDir}/data/mzml_spectra
	cp *.mzML ${projectDir}/data/mzml_spectra
	"""

}

workflow RAWMZML {

	bash_script = channel.fromPath("${projectDir}/data/bash_scripts/raw_to_mzml.sh")
	raw_files = channel.fromPath("${projectDir}/data/raw_files/*.raw").collect()

	RAW_TO_MZML(
		bash_script,
		raw_files
	)

}
