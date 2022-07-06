#Add local rules to aggregate data like htseq, ribo, exon, intron

rule merge_htseq:
	input:
		htseq=expand(data_folder+"htseq/{sample}_counts.txt", sample=get_samples()),
		exon=expand(data_folder+"exoncount/{sample}_exon.txt", sample=get_samples()),
		gene_ensg_table = os.path.join(workflow.basedir, "ref/human_ensemblId_geneName_table.tsv"),
	output:
		htseq=data_folder+"htseq_merged.csv",
		exon=data_folder+"exon_merged.csv",
		r2g=data_folder+"qc/reads_to_genes.txt"
	threads: 1
	resources: mem_mb=10000
	params:
		sample_regex=config["wildcard_constraints"],
		htseq_idx="0",
		htseq_colnames="gene_num",
		exon_idx="0,1,2",
		exon_colnames="exon_id,gene_num,exon_num"	
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	shell:
		"""
		python util/agg/merge_htseq.py {input.gene_ensg_table} '{input.htseq}' '{params.sample_regex}' {params.htseq_idx} {output.htseq} --reads_to_genes_path {output.r2g} --colnames {params.htseq_colnames} &&
		python util/agg/merge_htseq.py {input.gene_ensg_table} '{input.exon}' '{params.sample_regex}' {params.exon_idx} {output.exon} --colnames {params.exon_colnames} ;
		"""

rule merge_ribo:
	input:
		ribo=expand(data_folder+"ribocount/{sample}.RNA45SN5.counts.txt", sample=get_samples()),
	output:
		ribo=data_folder+"qc/ribo_frac.txt", 
	threads: 1
	resources: mem_mb=4000
	params:
		sample_regex=config["wildcard_constraints"],
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	script:
		os.path.join(workflow.basedir, "util/agg/summ_ribo_frac.py")

rule merge_dna:
	input:
		intron=expand(data_folder+"introncount/{sample}_intron.txt", sample=get_samples()),
		exon=expand(data_folder+"exoncount/{sample}_exon.txt", sample=get_samples()),
	output:
		dna=data_folder+"qc/intron_exon_ratios.txt", 
	threads: 1
	resources: mem_mb=4000
	params:
		sample_regex=config["wildcard_constraints"],
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	script:
		os.path.join(workflow.basedir, "util/agg/get_intron_exon_ratio.py")

rule merge_deg:
	input: expand(data_folder+"exoncount/{sample}_exon.txt", sample=get_samples())
	output:
		deg=data_folder+"qc/deg_3prime_bias_frac_1.txt"
	threads: 1
	resources: mem_mb=4000
	params:
		sample_regex=config["wildcard_constraints"],
	conda:
		os.path.join(workflow.basedir, "envs/count_qc.yaml")
	shell:
		"""
		python util/agg/summ_3prime_bias_ratio.py '{input}' {output.deg} --sample_regex '{params.sample_regex}'
		"""

