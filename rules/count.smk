def get_htseq_strandedness(wildcards):
	if stranded_type == "reverse":
		return htseq_params
	elif stranded_type == "forward":
		return htseq_params_forward
	return htseq_params_nostrand

rule htseq:
	input: data_folder+"dedup/{sample}.sortedByName.out.bam"
	output: 
		counts=temp(data_folder+"htseq/{sample}_counts.txt"),
		exon=temp(data_folder+"exoncount/{sample}_exon.txt"),
		intron=temp(data_folder+"introncount/{sample}_intron.txt")
	params:
		options=get_htseq_strandedness,
		genome_index=ref_anno,
		genome_index_exon=ref_anno_no_ERCC,
		genome_index_intron=ref_introns
	resources: mem_mb=lambda wildcards, attempt: attempt*21200
	threads: 1
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	shell:
		"""
		htseq-count {params.options} {input} {params.genome_index} > {output.counts} ;
		htseq-count {params.options} -t intron -i intron_id {input} {params.genome_index_intron} > {output.intron} ;
		htseq-count {params.options} -t exon -i exon_id --additional-attr=gene_id --additional-attr=exon_number {input} {params.genome_index_exon} > {output.exon}	
		"""
	
rule ribocount:
	input: data_folder+"dedup/{sample}.sortedByCoord.out.bam"
	output: 
		ribo=temp(data_folder+"ribocount/{sample}.RNA45SN5.counts.txt")
	resources: mem_mb=lambda wildcards, attempt: attempt*8192
	conda:
		os.path.join(workflow.basedir, "envs/trim_align_sort.yaml")
	#Issue https://github.com/conda-forge/r-base-feedstock/issues/67 produces error here sometimes by not being able to find ldpaths. To fix I commented out a line in the sh file that act R-base as sugg in issue
	shell:
		"""
		samtools view {input} {ribo_region} | cut -f 1 | sort | uniq | wc -l  > {output.ribo} &&
		samtools view {input} | cut -f 1 | sort | uniq | wc -l >> {output.ribo} &&
		Rscript util/qc/get_ribo_frac.R {output.ribo} ;
		"""
