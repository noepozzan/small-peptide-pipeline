[![ci](https://github.com/zavolanlab/zarp/workflows/CI/badge.svg?branch=dev)](https://github.com/zavolanlab/zarp/actions?query=workflow%3Aci)
[![GitHub license](https://img.shields.io/github/license/zavolanlab/zarp?color=orange)](https://github.com/zavolanlab/zarp/blob/dev/LICENSE)

<div align="left">
    <img width="20%" align="left" src=images/esel.webp>
</div> 

**Small Peptide Pipeline** ([Zavolan-Lab][zavolan-lab] whatever ... Pipeline) is a workflow that allows you
to analyze riboseq reads for the existence of small peptides and validating hits on peptidomics data.
The workflow relies on publicly available bioinformatics tools and currently handles (MERIC?)single-end stranded bulk RNA-seq and label-free peptidomics data.
The workflow is developed in [Nextflow][nextflow], a widely used workflow management system in the bioinformatics
community.

According to the current SMAPP implementation, reads are first pre-processed and then filtered against a library of rRNA.
Quality control with state-of-the-art tools gives you meaningful initial insights into the quality and composition of your ribo-Seq library.
After mapping of the ribo-Seq to your reference of choice, potential small peptides can be extracted, mainly using [Ribo-TISH][ribotish]. If you have experimental data that also comprises proteomics files, evaluation of your presumed small peptides is possible with [philosopher][philosopher].
Additional reports summarise the results of the individual steps and provide useful visualisations.

<div align="center">
    <img width="80%" src=images/flowchart.png>
</div> 


> **Note:** For a more detailed description of each step, please refer to the [workflow
> documentation][pipeline-documentation].


# Requirements

The workflow has been tested on:
- CentOS 7.5
- macOS 12.3.1

> **NOTE:**
> Currently, we only support **Mac & Linux** execution. 


# Installation

## 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone https://github.com/noepozzan/small_peptide_pipeline.git
cd small_peptide_pipeline
```

## 2. Conda and Mamba installation

Workflow dependencies can be conveniently installed with the [Conda][conda]
package manager. We recommend that you install [Miniconda][miniconda-installation] 
for your system (Linux). Be sure to select the Python 3 option. 
The workflow was built and tested with `miniconda 4.13.0`.
Other versions are not guaranteed to work as expected.

Given that Miniconda has been installed and is available in the current shell the first
dependency for ZARP is the [Mamba][mamba] package manager, which needs to be installed in
the `base` conda environment with:

```bash
conda install -y mamba -n base -c conda-forge
```

## 3. Dependencies installation

For improved reproducibility and reusability of the workflow,
each individual step of the workflow runs in its own [Singularity][singularity] or [Docker][docker]
container.
As a consequence, running this workflow has very few individual dependencies.
Since this pipeline depends on many different software tools, only **container execution** is possible. This requires Singularity or Docker to be installed on the system where the workflow is executed. 
As the functional installation of Singularity and Docker require root privilege, and Conda currently only provides Singularity for Linux architectures, the installation instructions are slightly different depending on your system/setup:

### For most users

If you do *not* have root privileges on the machine you want
to run the workflow on *or* if you do not have a Linux machine, please [install
Singularity][singularity-install] or [install Docker][docker-install] separately and in privileged mode, depending
on your system. You may have to ask an authorized person (e.g., a systems
administrator) to do that. This will almost certainly be required if you want
to run the workflow on a high-performance computing (HPC) cluster. 

> **NOTE:**
> The workflow has been tested with the following versions:  
>  * `Singularity v3.8.5-1.el7`
>  * `Docker 20.10.17`

After the installation has completed, install the remaining dependencies with:
```bash
mamba env create -f install/environment.yml
```

### As root user on Linux

(NOT RELEVANT RIGHT NOW. ASK MERIC)
If you have a Linux machine, as well as root privileges, (e.g., if you plan to
run the workflow on your own computer), you can execute the following command
to include Singularity in the Conda environment:

```bash
mamba env create -f install/environment.root.yml
```

## 4. Activate environment

Activate the Conda environment with:

```bash
conda activate small_peptides
```

# Extra installation steps (optional)

## 5. Non-essential dependencies installation

Most tests have additional dependencies. If you are planning to run tests, you
will need to install these by executing the following command _in your active
Conda environment_:

```bash
mamba env update -f install/environment.dev.yml
```

## 5. Before running the tests

It is important to know that this workflow relies on many external tools.
One of those is [MSFragger][msfragger].
Since MSFragger is only free for non-commercial use, you should run
`source data/scripts/echo_env.sh` from the main dir of this project.
This sets environment variables that allow you to pull the private MSFragger image from [noepozzan's dockerhub][dockerhub-np] repository.

## 6. Successful installation tests

I have prepared tests to check the integrity of the workflow and its components.
These can be run as follows:
	````bash
	cd <main directory of this project>
	nextflow run main.nf -profile test,<slurm,slurm_offline>
	```
ATTENTION:
Since even the testing files for this pipeline are quite large, I provide a github repo to pull from.
	````bash
	cd <main directory of this project>
	git clone git@github.com:noepozzan/small_peptide_pipeline_test_data.git
	mkdir -p data/tests/
	mv small_peptide_pipeline_test_data/* data/tests/
	```
Once you have done this, you can run the tests shown above.

The most critical of these tests enable you to execute the entire workflow on a 
set of small example input files. Note that for this and other tests to complete
successfully, [additional dependencies](#installing-non-essential-dependencies) need to be installed.

Remember to activate the conda environment!

Execute one of the following commands to run the test workflow on your local machine:
* Test workflow on local machine with **Singularity**:
(I need to test this if this works at all on my machine)
```bash
nextflow run main.nf -profile test,docker
```

Execute one of the following commands to run the test workflow 
on a [Slurm][slurm]-managed high-performance computing (HPC) cluster:
* Test workflow with **Singularity**:
```bash
nextflow run main.nf -profile test,<slurm,slurm_offline>
```


# Running the workflow on your own samples

If you want to run the workflow on your own files, running it is pretty straightforward:

	```bash
	cd <project's main directory>
	```
	````bash
	nextflow run main.nf -profile <profile of your choice>,<profile that fits your work environment>
	```

But before you start, you have to get the configuration right.
As you see above, this workflow needs 2 profiles:
	<profile of your choice>:  where you provide the files and parameters for the tools included in the workflow
	<profile that fits your work environment>: where you provide the details to the memory and the CPUs of your system

1. You have the choice of running the workflow in different configurations: (substitute one of the below options for the <profile of choice> above)
    full: to run the full pipeline (this is computationally quite heavy and should be done in cluster environment)
    test: to only run the test pipeline with small files
    qc: to only run the quality control part of the pipeline
    prepare: to prepare the reads
    ribotish: to only run ribotish
    proteomics: to quantify your proteomics files

While this looks quite straightforward up to this point, make sure to provide the right files for each of the run modes.
These files have to be provided, as follows:

In the project's root directory, there is a folder called `conf/`.
This folder houses all configuration files necessary to deal with the different run modes.
IMPORTANT: If you choose the profile <full>, you have to specify the paths to your files in the `conf/full.config` configuration file.
This works for all the other profiles in the same manner.
Use your editor of choice to populate these files with appropriate
paths. Every config files indicates the variables necessary to run this workflow in the way you want it to.

2. Have a look at the examples in the `conf/` directory to see what the
files should look like, specifically:

    - [slurm.config]
    - [docker.config]

    - For more details and explanations, refer to the [pipeline-documentation]

3. Pick one of the following choices for either local or cluster execution:

> **NOTE:** Depending on the configuration of your Slurm installation you may
> need to adapt the files under the `conf/` directory 
> and the arguments to options `memory` and `cpus`
> in the file `*.config` of the respective profile.
> Consult the manual of your workload manager as well as the section of the
> nextflow manual dealing with [profiles].

	slurm: for cluster execution (needs singularity installed)
	slurm_offline: for cluster execution
	(needs singularity installed and also needs you to first run:
	````bash
	cd <main directory of this project>
	bash data/scripts/pull_containers.sh
	```
	)
	docker: for local execution (needs docker installed and the daemon running)

4. Start your workflow run:

	Either, to view the output directly in your terminal:
    ```bash
	nextflow run main.nf -profile <profile of your choice>,<profile that fits your work environment>
    ```

	Or to have the workflow run in the background:
	(practical if you need to leave your computer while still running the pipeline)
	This option requires you to copy the exact nextflow command you intend to run into the `slurm.script`,
	which you'll find in the project's main directory.
	```bash
	sbatch slurm.script
	```

[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[profiles]: <https://www.nextflow.io/docs/latest/config.html#config-profiles>
[mamba]: <https://github.com/mamba-org/mamba>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[singularity]: <https://sylabs.io/singularity/>
[docker]: <https://docker.com/>
[msfragger]: <https://msfragger.nesvilab.org/>
[philosopher]: <https://github.com/Nesvilab/philosopher>
[nextflow]: <https://nextflow.io/>
[singularity-install]: <https://sylabs.io/guides/3.5/admin-guide/installation.html>
[docker-install]: <https://docs.docker.com/engine/install/>
[dockerhub-np]: <https://hub.docker.com/u/noepozzan>
[ribotish]: <https://bioinformatics.mdanderson.org/public-software/ribo-tish/>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[zavolan-lab]: <https://www.biozentrum.unibas.ch/research/researchgroups/overview/unit/zavolan/research-group-mihaela-zavolan/>
[slurm.config]: conf/slurm.config
[docker.config]: conf/docker.config
[pipeline-documentation]: pipeline_documentation.md

