#Get correct output of trimmomatic rule
import subprocess


"""
		STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.sample} --readFilesCommand zcat \
		--outSAMunmapped Within --outReadsUnmapped Fastx  --outFilterType BySJout --outSAMstrandField intronMotif \
		--outFileNamePrefix {params.out_prefix} --outSAMtype BAM {params.sam_type} {params.twopass} {star_params} {params.sj_filter}
		"""

rule = "STAR --runThreadN " + str(snakemake.threads) + " --readFilesIn "
params = " --readFilesCommand zcat --outSAMunmapped Within --outReadsUnmapped Fastx  --outFilterType BySJout --outSAMstrandField intronMotif --outFileNamePrefix " + snakemake.params.out_prefix + " --outSAMtype BAM " + snakemake.params.sam_type + " " + snakemake.params.twopass + " " + snakemake.star_params " " + snakemake.params.sj_filter

both_reads = (" ").join(snakemake.input.sample) 
just_r1 = snakemake.input.sample[0]

if seq_type == "pe":
	just_r2 = snakemake.input.sample[1]

rule_both = rule + both_reads + params 
subprocess.run(rule_both, shell=True)

if call_as:
	rule_r1 = rule + just_r1 + params 
print(rule)
subprocess.run(rule, shell=True)
