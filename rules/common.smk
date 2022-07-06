def patt_match_reads(sample, to_merge):
	suffix = "_001.fastq*"
	if merge_mult_lanes and to_merge:
		suffix = "_full.fastq.gz"
	reads_sample = {}

	folder_prefix = reads_folder + "/"
	for read in reads:
		key = str.lower(read)
		search_pattern = folder_prefix + "**/" + sample +"_*" + read + suffix
		if merge_mult_lanes and to_merge:
			val = folder_prefix + sample + "_" + read + suffix
		else:
			val = samples.loc[sample, key] if get_reads_from_sample_sheet else sorted(glob(search_pattern, recursive=True))
		reads_sample[key] = val
	return reads_sample

def get_reads(wildcards):
	dir_done = data_folder+"logs/make_log_folders.done"
	out = {"dirs":dir_done}
	out.update(patt_match_reads(wildcards["sample"], False))
	return out

def get_fastq(wildcards):
	dir_done = data_folder+"logs/make_log_folders.done"
	out = {"dirs":dir_done}
	out.update(patt_match_reads(wildcards["sample"], True))
	return out

def get_trimmed(wildcards):
	read_str = data_folder+"trimmed/{sample}_{read}.fastq.gz"
	return expand(read_str, read=reads, **wildcards)

def get_adapter_file():
	adapter_file = adapter_type.capitalize()

	if "seq" in adapter_file:
		idx_s = adapter_file.find("seq")
		adapter_file = adapter_file[:idx_s] + adapter_file[idx_s].capitalize() + adapter_file[(idx_s + 1):]
	if adapter_file == "Nextera":
		adapter_file += seq_type.upper()

	adapter_file += ("-" + seq_type.upper())

	if adapter_file == "TruSeq3-PE":
		# 2 adapter files for TruSeq3-PE should be second set since cfRNA usually isnt high quality and can contain adapters in strange places
		# See https://www.biostars.org/p/323087/ for further details
		adapter_file += "-2" 

	adapter_file += ".fa"
	full_file_path = os.path.join(workflow.basedir, "ref/adapters/" + adapter_file)

	# Ridiculously trimmomatic does not raise an error and fail if it cant find an adapter fa file. It just doesnt trim and returns an untrimmed file. 
	# Raise error and fail here instead to catch this case
	try:
		file = open(full_file_path, 'r')
	except FileNotFoundError as fnf_error: 
		print(fnf_error) 
		sys.exit(1)

	return full_file_path

rule make_dirs:
	input:
		dir_list = "ref/analysis_directory_struct.txt"
	output:
		touch(data_folder+"logs/make_log_folders.done")
	priority: 5
	threads: 1
	resources: mem_mb = 4000
	shell:
		'''
		sed -e "s#^#{data_folder}#" {input.dir_list} | xargs -d '\n' mkdir -p --
		'''

rule zcat:
	input: unpack(get_reads)
	output:
		r1 = protected(reads_folder+"/{sample}_R1_full.fastq.gz"),
		r2 = protected(reads_folder+"/{sample}_R2_full.fastq.gz"),
	threads: 1
	resources: mem_mb=5300
	run:
		read_cmd = "zcat " #"zcat " if "gz" in "{input.r1}" else "cat "
		
		if merge_mult_lanes:
			if seq_type == "pe":
				shell(read_cmd + "{input.r2} | gzip > {output.r2}")
			if seq_type == "se": #touch empty r2 so rule proceeds
				shell("touch {output.r2}")
			shell(read_cmd + "{input.r1} | gzip > {output.r1}")

rule sort_by_coord:
	input: "{prefix}.Aligned.out.{suffix}"
	output: "{prefix}.Aligned.sortedByCoord.out.{suffix,(bam|sam)}"
	threads: 1
	resources: mem_mb=lambda wildcards, attempt: attempt*20000
	conda:
		os.path.join(workflow.basedir, "envs/trim_align_sort.yaml")
	shell: "samtools sort -l 1 -o {output} -T {output}.tmp {input}"
		
