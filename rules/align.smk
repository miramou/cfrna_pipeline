rule star:
	input:
		sample = ancient(get_trimmed)
	output:
		aligned = data_folder+"star/{sample}.Aligned.{group}.bam",
	wildcard_constraints:
		group="sortedByCoord.out|out"
	threads: lambda wildcards, attempt: 6*attempt if attempt <=2 else 3
	resources: mem_mb=lambda wildcards, attempt: 32000*attempt if attempt<=2 else 32768
	params:
		sam_type=lambda wildcards, resources: "Unsorted" if resources.mem_mb==32768 else "SortedByCoordinate",
		#sj_filter=lambda config: star_sj_params if call_as else "",
		twopass=lambda config: star_2pass_params if is_two_pass else "",
		index=ref_index,
		out_prefix=data_folder+"star/{sample}.",
	conda:
		os.path.join(workflow.basedir, "envs/trim_align_sort.yaml")
	shell:
		"""
		STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.sample} --readFilesCommand zcat \
		--outSAMunmapped Within --outReadsUnmapped Fastx  --outFilterType BySJout --outSAMstrandField intronMotif \
		--outFileNamePrefix {params.out_prefix} --outSAMtype BAM {params.sam_type} {params.twopass} {star_params} 
		"""
		#{params.sj_filter}

