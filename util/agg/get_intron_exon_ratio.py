import pandas as pd
import numpy as np
import sys
import re

def get_intron_exon_ratio(intron_file, exon_file, regex_sample):
	sample = re.search(regex_sample, intron_file)

	intron = pd.read_csv(intron_file, sep="\t", index_col=0, names=["intron_id", "intron_count"])
	intron = intron.loc[~intron.index.str.contains("__")]

	exon = pd.read_csv(exon_file, sep="\t", index_col=0, names=["exon_id", "exon_num", "exon_count"])
	exon = exon.loc[~exon.index.str.contains("__")]

	intron_exon_ratio = round((np.sum(intron.intron_count))/(np.sum(exon.exon_count)+0.0001), 4)
	
	return pd.DataFrame([(sample.group(), intron_exon_ratio)], columns=["sample", "intron_exon_ratio"])

def main(intron_files, exon_files, regex_sample, out_path):
	intron_files = sorted(intron_files)
	exon_files = sorted(exon_files)

	for i in range(len(intron_files)):
		nxt_line = get_intron_exon_ratio(intron_files[i], exon_files[i], regex_sample)
		
		if i == 0:
			ratios = nxt_line
		else:
			ratios = pd.concat((ratios, nxt_line),axis=0)
	
	ratios.to_csv(out_path, sep="\t", index=False)

main(snakemake.input.intron, snakemake.input.exon, snakemake.params.sample_regex, snakemake.output.dna)

