#Get correct output of trimmomatic rule
import subprocess

trim_type = snakemake.params.trim_type.upper()

rule = "trimmomatic " + trim_type + " -threads " + str(snakemake.threads) + " "

files = (" ").join([str(snakemake.input.r1), str(snakemake.input.r2), str(snakemake.output.r1)]) if trim_type == "PE" else (" ").join([str(snakemake.input.r1), str(snakemake.output.r1)])

if trim_type == "PE":
	 files += " " + str(snakemake.output.u1) + " " + str(snakemake.output.r2) + " " + str(snakemake.output.u2)

if trim_type == "SE": #Create empty files so that snakemake rule proceeds
	for f in [str(snakemake.output.u1), str(snakemake.output.r2), str(snakemake.output.u2)]:
		make_f = "touch " + f
		subprocess.run(make_f, shell=True)

rule += files + " " + snakemake.params.trimmer

print(rule)
subprocess.run(rule, shell=True)
