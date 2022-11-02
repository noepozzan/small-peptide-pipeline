#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis


process WORKSPACE {

    label 'philosopher'

    input:
    path ribotish_speptide

	output:
	val 'workspace_initialized'

    script:
    """
	# idea: create a workspace to combine all results
	# and subworkspaces for each file
	rm -rf ${params.workspace}
	mkdir ${params.workspace}

	cp -R -n -p ${params.proteomics_reads} ${params.workspace}
	cp ${ribotish_speptide} ${params.workspace}
	
	cd ${params.workspace}


	philosopher workspace --clean
	philosopher workspace --init
	"""

}

process RAW_TO_MZML {

    label "thermorawfileparser"

    publishDir "${params.output_dir}/raw_to_mzml", mode: 'copy', pattern: '*.mzML'
    publishDir "${params.log_dir}/raw_to_mzml", mode: 'copy', pattern: '*.log'

    input:
    val workspace

    output:
	val 'raw_to_mzml_done', emit: done
    path '*.mzML*', emit: mzML, optional: true
    path '*.log*', emit: log, optional: true

    script:
    """
    workd=\$(pwd)
    cd ${params.workspace}

    count=`ls -1 *.raw *.RAW *.Raw 2>/dev/null | wc -l`
    if [ \$count != 0 ]
    then 
        ThermoRawFileParser.sh \
            --input_directory=./ \
            &> thermorawfileparser.log
    fi

    cp *.mzML thermorawfileparser.log \${workd}
    """
}
process DATABASE {

    label 'philosopher'
    
    publishDir "${params.philosopher_dir}/database", mode: 'copy', pattern: '*.fas*'
    publishDir "${params.log_dir}/database", mode: 'copy', pattern: '*.log'

    input:
    val raw_to_mzml
    path db

    output:
    path '*.fas*', emit: fas
	path '*.log', emit: log

    script:
    """
    workd=\$(pwd)
	# same pattern always:
	# change into working directory, execute command,
	# then copy the generated output files back to
	# nextflow directory to allow nextflow to track files
	cd ${params.workspace}
	#philosopher workspace --init
	philosopher database \
		--custom ${db} \
		--contam \
		&> database.log
	
    #cp -R -n -p *.fas database.log \$workd
	cp -R -n -p *-${db}* database.log \$workd
	"""

}

process GENERATE_CHANGE_PARAMS {

    label 'msfragger'

    publishDir "${params.philosopher_dir}/generate_change_params", mode: 'copy', pattern: 'msfragger.params'
    publishDir "${params.log_dir}/generate_change_params", mode: 'copy', pattern: '*.log'

    input:
    path db
    path change_file_script

    output:
    path 'msfragger.params', emit: params
	path '*.log', emit: log

    script:
    if (params.fragger_mode == "open")
    """
	# generate MSFRAGGER parameter file
    #java -jar /MSFragger.jar --config nonspecific

	# python script to change some parameters
	python ${change_file_script} \
		--database ${db} \
		--param "${projectDir}/data/open_fragger.params" \
		--out msfragger.params \
		&> generate_change_params.log

    cp msfragger.params ${params.workspace}
	"""
    else if (params.fragger_mode == "closed")
    """
	python ${change_file_script} \
		--database ${db} \
		--param "${projectDir}/data/open_fragger.params" \
		--out msfragger.params \
		&> generate_change_params.log

    cp msfragger.params ${params.workspace}
    """
    else if (params.fragger_mode == "nonspecific")
    """
	python ${change_file_script} \
		--database ${db} \
		--param "${projectDir}/data/nonspecific_fragger.params" \
		--out msfragger.params \
		&> generate_change_params.log

    cp msfragger.params ${params.workspace}
    """
    else if (params.fragger_mode == "glyco")
    """
	python ${change_file_script} \
		--database ${db} \
		--param "${projectDir}/data/Nglyco-HCD_fragger.params" \
		--out msfragger.params \
		&> generate_change_params.log

    cp msfragger.params ${params.workspace}
    """
}

process MSFRAGGER {

    label "msfragger"
    label "fragger_heavy"

    publishDir "${params.philosopher_dir}/msfragger/", mode: 'copy'
    publishDir "${params.log_dir}/msfragger/", mode: 'copy', pattern: '*.log'

    input:
    path fragger_params
    path db_file
	
    output:
    path '*.pepXML', emit: pepXML
	path '*.log', emit: log
    
    script:
	"""
	workd=\$(pwd)	

	cd ${params.workspace}
	java \
		-Xmx${params.fragger_ram}g \
		-jar /MSFragger.jar \
		${fragger_params} \
		*.mzML \
		&> msfragger.log

	cp *.pepXML msfragger.log \$workd

	"""

}

process PEPTIDEPROPHET {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/peptideprophet", mode: 'copy', pattern: 'interact.pep.xml'
    publishDir "${params.log_dir}/peptideprophet", mode: 'copy', pattern: '*.log'
       
    input:
    path db
    path pepXML
    
    output:
    path "interact.pep.xml", emit: combined_pepxml
    
	script:
	"""
	workd=\$(pwd)
	
	# peptide assignment validation on all files at once
	cd ${params.workspace}
	philosopher peptideprophet \
		${params.peptideprophet_args} \
		--database ${db} \
		*.pepXML \
		&> peptideprophet.log

	cp interact.pep.xml peptideprophet.log \$workd

	"""

}

process PROTEINPROPHET {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/proteinprophet", mode: 'copy', pattern: 'interact.prot.xml'
    publishDir "${params.log_dir}/proteinprophet", mode: 'copy', pattern: '*.log'
       
    input:
    path pepxml
    
    output:
    path "interact.prot.xml", emit: interact_prot_xml
	path "*.log", optional: true, emit: log

    script:
    """
    workd=\$(pwd)

    cd ${params.workspace}

    touch interact.prot.xml

    cp ${params.workspace}/interact.prot.xml \$workd
    """

}

process FILTER_FDR {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/filter_fdr", mode: 'copy', pattern: 'filter_done'
    publishDir "${params.log_dir}/filter_fdr", mode: 'copy', pattern: '*.log'

    input:
    path pepxml
    path protxml

    output:
    val "filter_done", emit: filter_done
	path "*.log", emit: log
 
    script:
    """
    workd=\$(pwd)
    cd ${params.workspace}

    philosopher filter \
        ${params.philosopher_filter_args} \
        --pepxml ${pepxml} \
        &> filter_fdr.log

    cp filter_fdr.log \$workd
    """

}

process FREEQUANT {

    label "philosopher"
    label "predicting"
    
	publishDir "${params.philosopher_dir}/freequant", mode: 'copy', pattern: 'freequant_done'
    publishDir "${params.log_dir}/freequant", mode: 'copy', pattern: '*.log'

    input:
    val filter_fdr

    output:
    val 'freequant_done', emit: freequant_done
	path '*.log', emit: log

    script:
    """
	workd=\$(pwd)
	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
    cd ${params.workspace}
    philosopher freequant \
		--dir . \
		2> freequant.log

	cp freequant.log \$workd
    """

}

process REPORT {

    label "philosopher"

	publishDir "${params.philosopher_dir}/report", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.log_dir}/report", mode: 'copy', pattern: '*.log'

    input:
    val freequant

    output:
    path 'msstats.tsv', emit: msstats
	path '*.tsv', emit: tsv
	path '*.log', emit: log
 
    script:
    """
	workd=\$(pwd)
    
	cd ${params.workspace}

	# reports about the findings for easy interpretation
	philosopher report \
		--msstats \
		--decoys \
		&> report.log

    mv msstats.csv msstats.tsv
	cp msstats.tsv peptide.tsv psm.tsv ion.tsv report.log \$workd

	"""

}


process IONQUANT {

    label "ionquant"

    publishDir "${params.philosopher_dir}/ionquant", mode: 'copy', pattern: '*_quant.csv'
    publishDir "${params.log_dir}/ionquant", mode: 'copy', pattern: '*.log'

    input:
    val report
	path pepXML

    output:
    path "*_quant.csv", emit: quant_csv
	path "*.log", emit: log

    script:
    """
	workd=\$(pwd)
    cd ${params.workspace}

	# extract parent dir
	dir=${params.proteomics_reads}
    parentdir="\$(dirname "\$dir")"
	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
	java -Xmx${params.ionquant_ram}G -jar /IonQuant.jar \
		--specdir \${parentdir} \
		--multidir . \
		*.pepXML \
		&> ionquant.log

	cp *_quant.csv ionquant.log \$workd

    """

}

workflow PHILOSOPHER {

    take:
    ribotish_predict_ch
    proteomics_reads

    main:
	WORKSPACE(
		ribotish_predict_ch,
	)

    RAW_TO_MZML(
        WORKSPACE.out.collect()
    )

	DATABASE(
		RAW_TO_MZML.out.done,
		ribotish_predict_ch
	)

	GENERATE_CHANGE_PARAMS(
		DATABASE.out.fas,
		params.change_params_script
	)

	MSFRAGGER(
		GENERATE_CHANGE_PARAMS.out.params,
		DATABASE.out.fas
	)
	
	PEPTIDEPROPHET(
		DATABASE.out.fas,
		MSFRAGGER.out.pepXML
	)
	
	PROTEINPROPHET(
		PEPTIDEPROPHET.out.combined_pepxml
	)

	FILTER_FDR(
		PEPTIDEPROPHET.out.combined_pepxml,
		PROTEINPROPHET.out.interact_prot_xml
	)
	
	FREEQUANT(
		FILTER_FDR.out.filter_done
	)
    
	REPORT(
		FREEQUANT.out.freequant_done
	)

	IONQUANT(
        REPORT.out.msstats,
        MSFRAGGER.out.pepXML
    )

    emit:
	report = REPORT.out.msstats
    ionquant = IONQUANT.out.quant_csv

}
