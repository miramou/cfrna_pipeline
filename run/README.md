# Setting up a project-specific run - Details

# Setting up a new project
Typically for each new project, you'll have relevant config files and read data/output in 2 locations. 
* Configuration files or scripts required to run this pipeline which can be found under "run/" within this pipeline directory
	* See "Creating pipeline run specific folder and files" below 
* Sample metadata, read data and output produced by this pipeline however are saved in an independent location - usually a project-specific within OAK
	* See "Creating a project-specific folder within OAK" below

## Creating pipeline run specific folder and files
For an example of all of these files that you can copy and modify where commented, see run/example

1. Change into this project repo directory.
	> cd path_prefix/cfRNA_pipeline
2. Create a new subdirectory titled "run/project_name" and a subsubdirectory titled "run/project_name/log"
	> mkdir run/project_name
	> mkdir run/project_name/log
3. Change into "run/project_name"
	> cd run/project_name
4. In this sub-directory, create the following files. For templates, copy from "run/example."

### Configuring the config file
---
Fields and their accepted values:
1. "samples" - Path to sample_metadata file **relative** to the root_dir specified in your common_config.yaml
	* Example - "root_dir/project_id/sample_metadata.csv"
2. "reads" - Path to folder containing .fastqs **relative** to the root_dir specified in your common_config.yaml
	* Example - "root_dir/project_id/raw_data"
3. "data" - Path to folder where all analysis will be saved **relative** to the root_dir specified in your common_config.yaml
	* Example - "root_dir/project_id/analysis/"
4. "wildcard_constraints" - Regex that specifies how to identify sample name
	* See https://docs.python.org/3/library/re.html for further details on REGEX's in python
5. "merge_mult_lanes" - Either "0" or "1". Specifies whether reads for a single sample are split over multiple lanes and should be merged (1) or not
	* Usually on NovaSeq - you have R1 and R2 reads split over a few lanes like S1_R1_L001.fastq.gz and S1_R1_L002.fastq.gz
	* This ensures that reads from multiple lanes be merged into a single file prior to further analysis
	* If all your reads were already on a single lane (e.g. NextSeq results or just single NovaSeq lane) then pass 0 here to skip this step
6. "is_two_pass" - Either "0" or "1". Specifies whether STAR should be run in 1-pass (0) or 2-pass mode (1)
	* 2-pass mode is only needed if calling SNPs or alternative splicing
7. "seq_type" - Either "se" (single-end) or "pe" (paired-end). Specifies nature of sequencing.
8. "stranded_type" - Either "reverse", "forward", or "unstranded". Specifies nature of cDNA construction.
9. "adapters" - 3 acceptable values - "truseq2", "truseq3", "nextera"
	* Used to choose the correct adapter file to trim adapter seqs
	* "truseq3" is most common
	* "truseq2" only for older expts on like GAII
	* "nextera" only for lib preps that used those adapters. Note only PE fa file available
	* To see all available adapter files - look in "ref/adapters/" which are the files that come with trimmomatic
	* For more details - check out https://www.biostars.org/p/323087/
10. "get_reads_from_sample_sheet" - Either "0" or "1". 
	* If 1, expects a sample_metadata sheet like "run/example/sample_metadata.full_file_paths.csv". See below for further details

### Create a run_snakemake.sh file
---
Change paths where necessary as specified by comments - namely 1st 4 variables defined

### Create a sample_metadata.csv file. 
---
This file can take 1 of 2 forms. This file should be saved **relative** to the root_dir specified in common_config (e.g. root_dir/project_id/sample_metadata.csv). If you prefer to change this preference, modify "Snakefile_modular.py" accordingly - specifically where the global var, samples, is defined.

* If the input fastq files include the sample name in the path, take the lazy route and make a file like sample_metadata.brief.csv
	* Example: If your fastqs are named like "S01_R001_L001.fastq.gz" and the sample name is "S01" then the pipeline is capable of using glob to find the relevant files without you specifying a full path
	* To see how this works, check out the fxn patt_match_reads in "rules/std/common.smk"
* If the input fastq files contain some random barcode that does not match your desired sample name, make a very specific file like sample_metadata.full_file_paths.csv
	* Example: If your fastqs are named like "U1028985_R001_L001.fastq.gz" and the sample name is "S01" then the pipeline needs to know that sample S1 corresponds to fastq titled "U1028985_R001_L001"

	
