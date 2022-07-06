rule trim:
	input: ancient(unpack(get_fastq))
	output:
		r1 = temp(data_folder+"trimmed/{sample}_R1.fastq.gz"),
		r2 = temp(data_folder+"trimmed/{sample}_R2.fastq.gz"),
		u1 = data_folder+"trimmed/{sample}_R1.U.fastq.gz", 		
		u2 = data_folder+"trimmed/{sample}_R2.U.fastq.gz",
	threads: 2
	resources: mem_mb=10000
	params:
		trim_type=seq_type,
		trimmer="ILLUMINACLIP:{}:{} {}".format(get_adapter_file(), adapter_params, trim_params) #See common.smk for get_adapter_file() fxn
	conda:
		os.path.join(workflow.basedir, "envs/trim_align_sort.yaml")
	script:
		os.path.join(workflow.basedir, "util/rules/trim.py") #To combine SE and PE reads

