import os
import sys
import pandas as pd
from glob import glob

#Common configuration
configfile: "util/config/common_config.yaml"

#Root dir
root_dir = config["root_dir"]

#Genome paths
ref_index = root_dir + config["ref"]["index"]
ref_anno = root_dir + config["ref"]["annotation"]
ref_anno_no_ERCC = root_dir + config["ref"]["annotation_no_ERCC"]
ref_introns = root_dir + config["ref"]["annotation_introns"]
ref_genes = root_dir + config["ref"]["annotation_genes"]

genome = root_dir + config["genome"]
ch_sizes = root_dir + config["chrom_size_file"]

#Adapter paths
adapter_params = config["adapter"]["params"]

#Tool params
trim_params = config["params"]["trim"]
star_params = config["params"]["star_all"] + config["params"]["star_all_tags"]
dedup_params = config["params"]["dedup"]
star_2pass_params = config["params"]["star_2pass"]
htseq_params = "-s reverse " + config["params"]["htseq"]
htseq_params_forward = "-s yes " + config["params"]["htseq"]
htseq_params_nostrand = "-s no " + config["params"]["htseq"]
ribo_region = config["params"]["ribo"]
multiqc_params = config["params"]["multiqc"]
multiqc_rseqc = config["params"]["multiqc_rseqc"]

#Unique read value
uniq_read_val = 60

#Common functions
def get_samples():
	return samples.index.tolist()

def get_read_names():
	reads = ["R1"]
	if seq_type == "pe":
		reads.append("R2")
	return reads

#For rule all
def get_output_files(wildcards):
	output = [data_folder+"qc/multiqc_report.html",
			data_folder+"htseq_merged.csv", 
			data_folder+"qc/outlier_data.txt"]

	return output
