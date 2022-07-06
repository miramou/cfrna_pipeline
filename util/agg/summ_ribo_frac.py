import pandas as pd
import numpy as np
import sys
import re

def get_ribo_frac(ribo_file, regex_sample):
	sample = re.search(regex_sample, ribo_file)

	ribo = pd.read_csv(ribo_file, sep="\t", names=["count"])
	ribo_frac = ribo.iloc[2, 0]
	
	return pd.DataFrame([(sample.group(), ribo_frac)], columns=["sample", "ribo_frac"])

def main(ribo_files, regex_sample, out_path):
	ribo_files = sorted(ribo_files)
	
	for i in range(len(ribo_files)):
		nxt_line = get_ribo_frac(ribo_files[i], regex_sample)
		
		if i == 0:
			ratios = nxt_line
		else:
			ratios = pd.concat((ratios, nxt_line),axis=0)
	
	ratios.to_csv(out_path, sep="\t", index=False)

main(snakemake.input.ribo, snakemake.params.sample_regex, snakemake.output.ribo)

