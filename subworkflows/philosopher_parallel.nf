nextflow.enable.dsl=2

//https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis

// run philosopher in parallel with subdirectories

process WORKSPACE_PARALLEL {

    label 'philosopher'
	echo true

    input:
    each(path(philosopher_done))
	each(path(ribotish_speptide))
	val dir

	output:
	val 'workspace_initialized'

    script:
    """
	workd=\$(pwd)

	rm -rf ${params.workspace}/${dir}
	mkdir -p ${params.workspace}/${dir}

	cd ${params.workspace}/${dir}
	philosopher workspace --clean
	philosopher workspace --init

	parent_dir=${params.proteomics_reads}
    parentdir="\$(dirname "\$parent_dir")"

	cp -R -n -p \${parentdir}/${dir}.mzML ${params.workspace}/${dir}
	cp -R -n -p \${workd}/${ribotish_speptide} ${params.workspace}/${dir}

	"""

}

process DATABASE_PARALLEL {

    label 'philosopher'
    
    publishDir "${params.philosopher_dir}/database_parallel", mode: 'copy', pattern: '*.fas*'
	publishDir "${params.log_dir}/database_parallel", mode: 'copy', pattern: "*.log"

    input:
    val workspace
	val dir
    each(path(db))

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
    cd ${params.workspace}/${dir}

	philosopher database \
		--custom ${db} \
		--contam
		&> ${dir}_database.log

    cp -R -n -p *.fas* *_database.log \$workd

	"""

}

process GENERATE_CHANGE_PARAMS_PARALLEL {

    label 'msfragger'

    publishDir "${params.philosopher_dir}/generate_change_params", mode: 'copy', pattern: 'msfragger.params'
    publishDir "${params.log_dir}/generate_change_params", mode: 'copy', pattern: '*.log'

    input:
    path db
    path change_file_script
	val dir

    output:
    path 'msfragger.params', emit: params
    path '*.log', emit: log

    script:
    """
    # generate MSFRAGGER parameter file
    java -jar /MSFragger.jar --config closed

    # python script to change some parameters
    python ${change_file_script} \
        --database ${db} \
        --param closed_fragger.params \
        --out msfragger.params \
        &> generate_change_params.log

    # copy params file to the working directory
    # since this process took place in some subdirectory
    cp msfragger.params ${params.workspace}/${dir}

    """

}

process MSFRAGGER_PARALLEL {

    label "msfragger"

    publishDir "${params.philosopher_dir}/msfragger_parallel/", mode: 'copy', pattern: '*.pepXML'
	publishDir "${params.log_dir}/msfragger_parallel/", mode: 'copy', pattern: '*.log'

    input:
    path closed_fragger
	path db_file
	val dir

    output:
    path '*.pepXML', emit: pepXML
    path '*.log', emit: log

    script:
	"""
	workd=\$(pwd)	

	cd ${params.workspace}/${dir}
	java \
		-Xmx${params.fragger_ram}g \
		-jar /MSFragger.jar \
		${closed_fragger} \
		${dir}.mzML \
		&> ${dir}_msfragger.log

	cp *.pepXML *_msfragger.log \$workd

	"""

}

process PEPTIDEPROPHET_PARALLEL {

	echo true
    label "philosopher"

    publishDir "${params.philosopher_dir}/peptideprophet_parallel", mode: 'copy', pattern: 'interact*.pep.xml'
    publishDir "${params.log_dir}/peptideprophet_parallel", mode: 'copy', pattern: '*.log'
    
	input:
    path db
	path pepXML
	val dir

    output:
	path "interact*.pep.xml", emit: pepxml
	path "*.log", emit: log

	script:
	"""
	workd=\$(pwd)
	cd ${params.workspace}/${dir}

	# peptide assignment validation on all files at once
	
	philosopher peptideprophet \
        ${params.peptideprophet_args} \
		--combine \
        --database ${db} \
        *.pepXML \
		--output "interact-${dir}" \
		&> ${dir}_peptideprophet.log

	cp interact*.pep.xml *_peptideprophet.log \$workd

	"""

}

process PROTEINPROPHET_PARALLEL {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/proteinprophet_parallel", mode: 'copy', pattern: 'interact.prot.xml'
    publishDir "${params.log_dir}/proteinprophet_parallel", mode: 'copy', pattern: '*.log'
       
    input:
    path pepxml
	val dir
    
    output:
    path "interact.prot.xml", emit: interact_prot_xml
	path "*.log", optional: true, emit: log

	script:
    if( params.skip_proteinprophet == true )
        """
		workd=\$(pwd)

		cd ${params.workspace}/$dir

        touch interact.prot.xml

		cp -n ${params.workspace}/${dir}/interact.prot.xml \$workd

	    """

    else if( params.skip_proteinprophet == false )
        """
        workd=\$(pwd)

        cd ${params.workspace}/${dir}

		# group peptides by their corresponding protein(s)
		# to compute probabilities that those proteins were
		# present in in the original sample
		philosopher proteinprophet \
			${dir}.mzML \
			&> ${dir}_proteinprophet.log

        cp interact.prot.xml *_proteinprophet.log \$workd

        """

}

process FILTER_FDR_PARALLEL {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/filter_fdr_parallel", mode: 'copy', pattern: 'filter_done'
    publishDir "${params.log_dir}/filter_fdr_parallel", mode: 'copy', pattern: '*.log'

    input:
    path pepxml
    path protxml
	val dir
    //tuple path(pepxml), val(dir)

    output:
    val "filter_done", emit: filter_done
	path "*.log", emit: log
 
    script:
    if( params.skip_proteinprophet == true )
        """
		workd=\$(pwd)
        cd ${params.workspace}/${dir}

        philosopher filter \
			${params.philosopher_filter_args} \
			--pepxml *.pep.xml \
			#--pepxml ${pepxml} \
            &> ${dir}_filter_fdr.log

		cp *_filter_fdr.log \$workd
        """
	else if( params.skip_proteinprophet == false && params.combine_peptideprophet == true )
		"""
		workd=\$(pwd)
		cd ${params.workspace}/${dir}

		philosopher filter \
			${params.philosopher_filter_args} \
			--pepxml interact-${dir}.pep.xml \
			--protxml interact.prot.xml \
			&> ${dir}_filter_fdr.log

		cp *_filter_fdr.log \$workd
		"""
    else if( params.skip_proteinprophet == false && params.combine_peptideprophet == false )
        """
		workd=\$(pwd)
        cd ${params.workspace}/${dir}

		# filter matches and estimate FDR
		# skip the --sequential parameter
		philosopher filter \
			--psm 0.05 \
			--ion 0.05 \
			--pep 0.05 \
			--prot 1 \
			--razor \
			--picked \
			--tag rev_ \
			--pepxml interact-0*.pep.xml \
			--protxml interact.prot.xml \
			&> filter_fdr.log

		cp *_filter_fdr.log \$workd
        """

}

process FREEQUANT_PARALLEL {

    label "philosopher"

    publishDir "${params.philosopher_dir}/freequant_parallel", mode: 'copy', pattern: 'freequant_done'
    publishDir "${params.log_dir}/freequant_parallel", mode: 'copy', pattern: '*.log'

    input:
    val filter_fdr
	val dir

    output:
    val 'freequant_done', emit: freequant_done
	path '*.log', emit: log

    script:
    """
	workd=\$(pwd)
	
	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
	
	cd ${params.workspace}/${dir}
    philosopher freequant \
		--dir \
		&> ${dir}_freequant.log

	cp *_freequant.log \$workd
    """

}

process REPORT_PARALLEL {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/report_parallel", mode: 'copy', pattern: '*msstats.csv'
	publishDir "${params.philosopher_dir}/report_parallel", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.log_dir}/report_parallel", mode: 'copy', pattern: '*.log'

    input:
    val freequant
	val dir

    output:
    path '*msstats.csv', emit: msstats
	path '*.tsv', emit: tsv
	path '*.log', emit: log
 
    script:
    """
	workd=\$(pwd)
	cd ${params.workspace}/${dir}

	# reports about the findings for easy interpretation
	philosopher report \
		--msstats \
		--decoys \
		&> ${dir}_report.log

	mv msstats.csv ${dir}_msstats.csv
    mv ion.tsv ${dir}_ion.tsv
    mv peptide.tsv ${dir}_peptide.tsv
    mv psm.tsv ${dir}_psm.tsv
	cp *msstats.csv *peptide.tsv *psm.tsv *ion.tsv *report.log \$workd

	"""

}


process IONQUANT_PARALLEL {

    label "ionquant"

	publishDir "${params.philosopher_dir}/ionquant_parallel", mode: 'copy', pattern: '*_quant.csv'
    publishDir "${params.log_dir}/ionquant_parallel", mode: 'copy', pattern: '*.log'

    input:
    val filter_fdr
	path pepXML
	val dir

    output:
    path "*_quant.csv", emit: quant_csv
	path "*.log", emit: log

    script:
    """
	workd=\$(pwd)
    cd ${params.workspace}/${dir}

	# extract parent dir
	dir=${params.proteomics_reads}
    parentdir="\$(dirname "\$dir")"
	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
	java -Xmx${params.ionquant_ram}G -jar /IonQuant.jar \
		--specdir \${parentdir} \
		--multidir . \
		*.pepXML \
		&> ${dir}_ionquant.log

	cp *_quant.csv *_ionquant.log \$workd

    """

}

process COMBINE {                                                               
                                                                                
    label "sorf_to_speptide"                                                    
                                                                                
    publishDir "${params.philosopher_dir}/combine", mode: 'copy', pattern: "msstats.csv"

    input:    
    path msstats                                                              
                                                                                
    output:                                                                     
    path "msstats.csv", emit: msstats_combined
                                                                                
    script:                                                                     
    """                                                                         
    WRITE_HEADER="true"                                                         
    for VAR in ${msstats}                                                     
    do                                                                          
        if [ "\$WRITE_HEADER" == "true" ]; then                                 
            head -1 \$VAR > msstats.csv
            WRITE_HEADER="false"                                                
        fi                                                                      
        tail -n +2 -q \$VAR >> msstats.csv
    done                                                                        
                                                                                
    """  

}

process CLEAN_UP_WORKSPACE {

	label "philosopher"

	input:
	path report

	script:
	"""

	cd ${params.workspace}
	cd ..
	rm -rf ${params.workspace}

	"""

}

workflow PHILOSOPHER_PARALLEL {

    take:
	philosopher_done
	msfragger_params
    ribotish_predict_ch
    proteomics_reads

    main:
	channel.fromPath("${params.proteomics_reads}")
        .map{ it.baseName.split("\\.")[0] }
        .set{ sub_dirs }

	WORKSPACE_PARALLEL(
		philosopher_done,
		ribotish_predict_ch,
		sub_dirs,
	)

    DATABASE_PARALLEL(
		WORKSPACE_PARALLEL.out.collect(),
		sub_dirs,
		ribotish_predict_ch
	)

	GENERATE_CHANGE_PARAMS_PARALLEL(
		DATABASE_PARALLEL.out.fas,
        params.change_params_script,
		sub_dirs
	)

	MSFRAGGER_PARALLEL(
        //msfragger_params,
        GENERATE_CHANGE_PARAMS_PARALLEL.out.params,
		DATABASE_PARALLEL.out.fas,
		sub_dirs
    )
	
	PEPTIDEPROPHET_PARALLEL(
		DATABASE_PARALLEL.out.fas,
		MSFRAGGER_PARALLEL.out.pepXML,
		sub_dirs
	)

	PROTEINPROPHET_PARALLEL(
		PEPTIDEPROPHET_PARALLEL.out.pepxml,
		sub_dirs
	)

    PEPTIDEPROPHET_PARALLEL.out.pepxml
        .combine(sub_dirs)                                                          
        .filter{ it[0].baseName.split("\\.")[0].replaceAll( /^interact-/, '' ) == it[1] }
        /*                                                                          
        .multiMap { it ->                                                           
                pepxml: it[0]                                                       
                dir: it[1]                                                          
        }                                                                           
        */                                                                          
        .set{ pepxml_dirs }

	FILTER_FDR_PARALLEL(
		PEPTIDEPROPHET_PARALLEL.out.pepxml,
		PROTEINPROPHET_PARALLEL.out.interact_prot_xml,
		sub_dirs
        //pepxml_dirs
	)
	FREEQUANT_PARALLEL(
		FILTER_FDR_PARALLEL.out.filter_done,
		sub_dirs
	)
    
	REPORT_PARALLEL(
		FREEQUANT_PARALLEL.out.freequant_done,
		sub_dirs
	)
    report = REPORT_PARALLEL.out.msstats

	IONQUANT_PARALLEL(
        FILTER_FDR_PARALLEL.out.filter_done,
        MSFRAGGER_PARALLEL.out.pepXML,
		sub_dirs
    )

    COMBINE(
        report.collect()
    )
    msstats_combined = COMBINE.out.msstats_combined

	/*
	CLEAN_UP_WORKSPACE(
		IONQUANT_PARALLEL.out.quant_csv.collect()
	)
	*/

    emit:
	msstats_combined

}
