def get_aligned_coord(wildcards):
	return  expand([data_folder+"star/{sample}.Aligned.sortedByCoord.out.bam"], **wildcards)

rule dedup:
	input: get_aligned_coord
	output:
		bam=temp(data_folder+"dedup/{sample}.bam"),
		metrics=data_folder+"dedup/{sample}.metrics.txt"
	resources: mem_mb=lambda wildcards, attempt: 30680*attempt #10240*attempt
	threads: 1
	conda:
		os.path.join(workflow.basedir, "envs/dedup.yaml")
	shell:
		"gatk MarkDuplicates {dedup_params} -I {input} -O {output.bam} -M {output.metrics}"


rule sort:
	input: data_folder+"dedup/{sample}.bam"
	output:
		sorted_name = temp(data_folder+"dedup/{sample}.sortedByName.out.bam"),
		sorted_coord = data_folder+"dedup/{sample}.sortedByCoord.out.bam",
	threads: 1
	resources: mem_mb=lambda wildcards, attempt: attempt*21200 	
	params: 
		old_tmp_regex=data_folder+"dedup/{sample}*tmp*"
	conda:
		os.path.join(workflow.basedir, "envs/trim_align_sort.yaml")
	shell:
		#Delete prior tmp files otherwise cant relaunch job with more memory
		"""
		find {params.old_tmp_regex} -delete
		samtools index {input} && 
		samtools sort -n {input} -o {output.sorted_name} ; 
		samtools sort {input} -o {output.sorted_coord} && 
		samtools index {output.sorted_coord} 
		"""
