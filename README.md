# cfRNA_pipeline
Pipeline to generate counts table from cfRNA data

# Overview
Full pipeline to process cfRNAseq data starting with fastq files

# Getting started

1. Install miniconda if you do not have it already - https://docs.conda.io/en/latest/miniconda.html
2. Git clone this repo
3. Change into this directory
4. Either:

	(A) Install snakemake into your base env
	> conda install -c bioconda snakemake=5.8.1=0 

	(B) Create a virtual environment with snakemake.
	> conda env create -f envs/smk.yaml 
	
	(C) Create a virtual environment with snakemake-minimal. Under 1GB in size for space issues.
	> conda env create -f envs/smk_minimal.yaml 

5. Modify "util/config/common_config.yaml" such that the paths specified work for you

# Setting up a new project
Typically for each new project, you'll have relevant config files and read data/output in 2 locations. 
* Configuration files or scripts required to run this pipeline which can be found under "run/" within this pipeline directory
	* See "Creating pipeline run specific folder and files" below 
* Sample metadata, read data and output produced by this pipeline however are saved in an independent location - usually a project-specific within OAK
	* See "Creating a project-specific folder within OAK" below

## Creating a project-specific folder and sample_metadata file within your file system
Within the root_dir you specified in common_config.yaml (typically your OAK dir if quake lab)

1. Create a project specific folder e.g. "root_dir/project_id/"
2. Change into that folder
3. Download raw reads from the cloud (eg AWS, GCloud, Etc) into a sub-directory e.g. "raw_data/"
4. Create an analysis specific directory e.g. "analysis/"
	* This is where all the pipeline specific files will be saved like the output from STAR, htseq, etc
5. In this same folder, create a sample_metadata file.

### Create a sample_metadata.csv file. 
---
This file can take 1 of 2 forms. For examples, see run/example

* If the input fastq files include the sample name in the path, take the lazy route and make a file like the example in run/example/sample_metadata.brief.csv
	* Example: If your fastqs are named like "S01_R001_L001.fastq.gz" and the sample name is "S01" then the pipeline is capable of using glob to find the relevant files without you specifying a full path
	* To see how this works, check out the fxn patt_match_reads in "rules/std/common.smk"
* If the input fastq files contain some random barcode that does not match your desired sample name, make a very specific file like the example in run/example/sample_metadata.full_file_paths.csv
	* Example: If your fastqs are named like "U1028985_R001_L001.fastq.gz" and the sample name is "S01" then the pipeline needs to know that sample S1 corresponds to fastq titled "U1028985_R001_L001"

## Creating pipeline run specific folder and files
For an example of both of these files that you can copy and modify where commented, see run/example

1. Change into this project repo directory.
	> cd path_prefix/cfRNA_pipeline
2. Create a new subdirectory titled "run/project_name" and a subsubdirectory titled "run/project_name/log"
	> mkdir run/project_name
	> mkdir run/project_name/log
3. Change into "run/project_name"
	> cd run/project_name
4. In this sub-directory, create the following files. For templates, copy from "run/example."

### Create a file titled config.yaml
---
This file has the following fields:

* "samples" - Path to sample_metadata file **relative** to the root_dir specified in your common_config.yaml - e.g. root_dir/project_id/sample_metadata.csv
* "reads" - Path to folder containing .fastqs **relative** to the root_dir specified in your common_config.yaml
* "data" - Path to folder where all analysis will be saved **relative** to the root_dir specified in your common_config.yaml
* "wildcard_constraints" - Regex that specifies how to identify sample name
* "merge_mult_lanes" - Either "0" or "1". Specifies whether reads for a single sample are split over multiple lanes and should be merged (1) or not
	* Usually on NovaSeq - you have R1 and R2 reads split over a few lanes like S1_R1_L001.fastq.gz and S1_R1_L002.fastq.gz
	* This ensures that reads from multiple lanes be merged into a single file prior to further analysis
	* If all your reads were already on a single lane (e.g. NextSeq results or just single NovaSeq lane) then pass 0 here to skip this step
* "is_two_pass" - Either "0" or "1". Specifies whether STAR should be run in 1-pass (0) or 2-pass mode (1)
	* 2-pass mode is only needed if calling SNPs or alternative splicing
* "seq_type" - Either "se" (single-end) or "pe" (paired-end). Specifies nature of sequencing.
* "stranded_type" - Either "reverse", "forward", or "unstranded". Specifies nature of cDNA construction.
* "adapters" - 3 acceptable values - "truseq2", "truseq3", "nextera"
	* Used to choose the correct adapter file to trim adapter seqs
	* "truseq3" is most common
	* "truseq2" only for older expts on like GAII
	* "nextera" only for lib preps that used those adapters. Note only PE fa file available
	* To see all available adapter files - look in "ref/adapters/" which are the files that come with trimmomatic
	* For more details - check out https://www.biostars.org/p/323087/
* "get_reads_from_sample_sheet" - Either "0" or "1". 
	* If 1, expects a sample_metadata sheet like "run/example/sample_metadata.full_file_paths.csv". See below for further details

### Create a run_snakemake.sh file
---
Change paths where necessary as specified by comments - namely 1st 4 variables defined

# Launching the pipeline
Once you've created all project-specific files and folders, you can launch the pipeline as follows

1. Change into this main directory
2. Perform a dry run to confirm everything works as expected for your project
	> bash run/project_name/run_snakemake.sh
3. If all goes as expected, launch snakemake by either:

	(A) If total sample number < 100
	> bash run/project_name/run_snakemake.sh unlock #Only run this if previously launched pipeline
	>
	> screen
	>
	> bash run/project_name/run_snakemake.sh snakemake
	>
	> Press Ctrl-A-D #Detachs from screen

	(B) If total sample number > 100
	> bash run/project_name/run_snakemake.sh unlock #Only run this if previously launched pipeline
	>
	> bash run/project_name/run_snakemake.sh sbatch #Submits a job to run snakemake continuously

4. Once all jobs are done, analyze data :) Maybe debug and check log files if jobs fail

# Output
Within the "data" directory you specified in the project-specific "config.yaml" file, you will find:
1. "htseq_merged.csv" - Full counts table merged across all samples
2. "exon_merged.csv" - Full exoncount table merged across all samples
3. "qc/multiqc_report.html" - MultiQC report
4. "qc/outlier_data.txt" - Aggregated QC metrics, plots are produced too in "plots/"
	* Lists raw value for each QC metric and whether the sample is an outlier for that metric. 	
	* Final two columns call whether sample is an outlier overall and what was one flag that failed for that sample
5. Sample-specific intermediate files in subdirectories
	* E.g. STAR-mapped bam files within folder "star"

# Notes
Sometimes you only want to rerun part of the pipeline (eg. rerun htseq and merge_htseq but not multiqc). You can do this by editing "run_snakemake.sh" as follows.

1. Navigate to the "CALL_SMK" fxn. Append the following to the end of the snakemake command in the body of "CALL_SMK"

	> --allowed-rules htseq merge_htseq all

2. Stanford Sherlock uses SLURM as job manager - yours may be different and require modification to the way job calls are made.
	
