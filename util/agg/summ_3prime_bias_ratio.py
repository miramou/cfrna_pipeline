import pandas as pd
import numpy as np
import sys
import re
import os
import argparse

def argparser():
    parser = argparse.ArgumentParser(description = "Estimate RNA degradation: Obtain fraction of genes for which all reads map to 3' most exon")
    parser.add_argument("bam_files", type = str, help = "Space-separated list of all bam files to process. Make sure the list is in quotations! Eg 'bam1.bam bam2.bam bam3.bam'")
    parser.add_argument("out_path", type = str, help = "Path to tsv file name to save output")
    parser.add_argument("--sample_regex", "-sr", default = None, type = str, help = "Sample regex to pull sample name from bam file name. Make sure this is in quotations. If not provided will use entire bam_file name")
    return parser.parse_args()

def read_exon_htseq_file(file_path):
	df = pd.read_csv(file_path, sep="\t", header=None, names=["exon_id", "gene_id", "exon_num", "exon_count"])
	return df.loc[~df.exon_id.str.contains("__")].set_index("exon_id")

def calc_bias_frac(htseq_df):
	exon_sum = htseq_df.groupby(["gene_id"])['exon_count'].sum()
	has_det_genes = np.sum(exon_sum > 0) 
	det_genes = exon_sum.loc[(exon_sum > 0)].index #genes with counts > 0

	htseq_df = htseq_df.loc[htseq_df.gene_id.isin(det_genes.to_list())]
	exon_sum = exon_sum.loc[det_genes]

	frac_all_3prime = np.nan
	
	if has_det_genes: #At least 1 gene has count > 0

		idx_3prime_exon = htseq_df.groupby(["gene_id"])['exon_num'].idxmax() #the 3' most exon is that with the highest exon_num

		ratio_3prime = htseq_df.loc[idx_3prime_exon].reset_index().set_index('gene_id')
		ratio_3prime.insert(ratio_3prime.shape[1], "ratio", ratio_3prime.exon_count / exon_sum) #number of reads for 3' exon / total reads for gene = ratio

		frac_all_3prime = np.sum(ratio_3prime.ratio == 1.0)
		frac_all_3prime /= ratio_3prime.shape[0]
	
	return frac_all_3prime

def get_bias_frac(ratio_file, sample_regex):

	bias_frac = calc_bias_frac(read_exon_htseq_file(ratio_file))
	sample = re.search(sample_regex, ratio_file).group() if sample_regex is not None else ratio_file
	
	return pd.DataFrame([(sample, bias_frac)], columns=["sample", "bias_frac"])

def main():
	args = argparser()

	ratio_files = sorted(args.bam_files.split(" "))

	for i in range(len(ratio_files)):
		curr_df = get_bias_frac(ratio_files[i], args.sample_regex)
		ratios = curr_df if i == 0 else pd.concat((ratios, curr_df), axis=0)
		if i > 0 and i % 10 == 0:
			os.system("echo Now processed " + str(i) + " files.")
	
	ratios.to_csv(args.out_path, sep="\t", index=False)
	return


if __name__ == "__main__":
	main()

