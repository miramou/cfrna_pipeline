def get_rseqc_strandedness(wildcards):
	if stranded_type == "reverse":
		return '1+-,1-+,2++,2--'
	elif stranded_type == "forward":
		return '1++,1--,2+-,2-+'
	return None

rule blast_adapter_seq:
	input: unpack(get_fastq)
	output:
		out=data_folder+"blast_adapter_seq/{sample}.sam",
	threads: 4
	params:
		adapter=get_adapter_file(), #See common.smk for get_adapter_file() fxn
		extra="-no_unaligned"
	resources: mem_mb=10000
	run:
		if seq_type == "pe":
			shell("magicblast -query {input.reads[0]} -query_mate {input.reads[1]} -subject {params.adapter} -infmt fastq -out {output.out} {params.extra} -num_threads {threads};")
		

rule fastqc:
	input: unpack(get_fastq)
	output:
		temp(touch(data_folder+"qc/fastqc/{sample}.fastqc.done"))
	threads: 6
	params: 
		outdir=data_folder+"qc/fastqc",
		seq = seq_type
	resources: mem_mb=32000
	conda:
		os.path.join(workflow.basedir, "envs/fastqc.yaml")
	script:
		os.path.join(workflow.basedir, "util/rules/fastqc.py")

#RSeQC rules modified from Snakemake github - https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/rules/qc.smk

rule rseqc_gtf2bed:
	input: ref_anno	   
	output:
		bed=data_folder+"qc/rseqc/annotation.bed",
		db=temp(data_folder+"qc/rseqc/annotation.db")
	threads: 1
	resources: mem_mb=10240
	conda:
		os.path.join(workflow.basedir, "envs/gffutils.yaml")
	script:
		os.path.join(workflow.basedir, "util/qc/gtf2bed.py")

rule rseqc:
	input:
		bam=data_folder+"dedup/{sample}.sortedByCoord.out.bam",
		bed=ancient(data_folder+"qc/rseqc/annotation.bed")
	output:
		stats=data_folder+"qc/rseqc/{sample}/{sample}.stats.txt",
		read_dist=data_folder+"qc/rseqc/{sample}/{sample}.readdistribution.txt", 
		all_done=temp(data_folder+"qc/rseqc/{sample}/{sample}.rseqc.done"),
	log: 
		data_folder+"qc/rseqc/{sample}/{sample}.out"
	params:
		prefix=data_folder+"qc/rseqc/{sample}/{sample}",
		gene_body_prefix=data_folder+"qc/rseqc/{sample}/{sample}.gene_body_cov",
		strand=get_rseqc_strandedness
	resources: mem_mb=lambda wildcards, attempt: 10240*attempt if attempt<=2 else 32768
	conda:
		os.path.join(workflow.basedir, "envs/rseqc.yaml")
	shell:
		# Removed RPKM sat - Takes over a day to run and excessive memory consumption leading to OOM error for about ~5% of samples (Based on BC_DEN set)
		# RPKM_saturation.py -r {input.bed} -d {params.strand} -i {input.bam} -o {params.prefix} > {log} 2>&1 ;
		"""
		junction_annotation.py r'-q {uniq_read_val}' -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1 ;
		junction_saturation.py r'-q {uniq_read_val}' -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1 ;
		bam_stat.py -i {input} > {output.stats} 2> {log} ;
		inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1 ;
		read_distribution.py -r {input.bed} -i {input.bam} > {output.read_dist} 2> {log} ;
		geneBody_coverage.py -r {input.bed} -i {input.bam} -o {params.gene_body_prefix} > {log} 2>&1 && 
		touch {output.all_done} ;
		"""

rule multiqc:
	input: 
		htseq_done=expand(data_folder+"htseq/{sample}_counts.txt", sample=get_samples()),
		fastqc_done=expand(data_folder+"qc/fastqc/{sample}.fastqc.done", sample=get_samples()),
		#rseqc_done=expand(data_folder+"qc/rseqc/{sample}/{sample}.rseqc.done", sample=get_samples())
	output:
		multiqc_report = data_folder+"qc/multiqc_report.html",
		#rseqc_report = data_folder+"qc/rseqc/multiqc_report.html"
	threads: 1
	resources: mem_mb=lambda wildcards, attempt: 10240*attempt if attempt <=2 else 32768
	params:
		multiqc_dir=data_folder,
		rseqc_dir=data_folder+"qc/rseqc/"
	conda:
		os.path.join(workflow.basedir, "envs/multiqc.yaml")
	shell:
		"""
		multiqc -o {params.multiqc_dir} -d {params.multiqc_dir} {multiqc_params} &&
		mv {params.multiqc_dir}multiqc_report.html {output.multiqc_report} ;
		echo 'Created MultiQC report' ;
		"""
		#   multiqc -o {params.rseqc_dir} -d {params.rseqc_dir} {multiqc_rseqc} &&
		#   echo 'Created RSeQC report.' ;
		#   """

rule outlierqc:
	input: 
		dna=data_folder+"qc/intron_exon_ratios.txt",
		r2g=data_folder+"qc/reads_to_genes.txt",
		deg=data_folder+"qc/deg_3prime_bias_frac_1.txt",
	output:
		outlier_report = data_folder+"qc/outlier_data.txt",
	threads: 1
	resources: mem_mb=4000
	params:
		out_dir=data_folder+"qc/",
		cutoffs=os.path.join(workflow.basedir, "ref/outlier_cutoffs_95p_w_r2g0_r1.txt")
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	shell:
		"""
		python util/qc/id_outliers.py {input.deg} {input.dna} {input.r2g} {params.out_dir} --cutoffs_file {params.cutoffs}
		"""
