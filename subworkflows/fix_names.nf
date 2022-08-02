nextflow.enable.dsl=2

process CONDITION {

	echo true

	label "philosopher"

	input:
	path csv
	path python_script

	script:
	"""
	dir=${params.proteomics_reads}
    parentdir="\$(dirname "\$dir")"

	python ${python_script} \
		--map_file ${csv} \
		--dir \${parentdir} \
		--mode rename # "rename" or "reverse"

	"""
}

workflow FIX_NAMES {

	take:
	csv_file
	py_script

	main:
	CONDITION(
		csv_file,
		py_script
	)

}
