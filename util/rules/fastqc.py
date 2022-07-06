#Get correct output of trimmomatic rule
import subprocess

rule = "fastqc -q -t " + str(snakemake.threads) + " --outdir " + snakemake.params.outdir + " " + str(snakemake.input.r1)

if snakemake.params.seq == 'pe':
	rule += " " + str(snakemake.input.r2)

print(rule)
subprocess.run(rule, shell = True)
