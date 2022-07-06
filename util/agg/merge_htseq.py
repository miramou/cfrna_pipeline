import pandas as pd
import re
import argparse

def argparser():
	parser = argparse.ArgumentParser(description = "Merge htseq like output")
	parser.add_argument("ensg_gene_file", type = str, help = "Path to ref file that matches ENSG IDs with gene names"),
	parser.add_argument("file_list", type = str, help = "Comma sep string containing all file names to merge"),
	parser.add_argument("sample_regex", type = str, help = "Regex to ID sample name from file name"),
	parser.add_argument("index_cols", type = str, help = "Comma sep string containing integer locations of index columns"),
	parser.add_argument("--ensg_exon_file", type = str, default = "", help = "Path to file containing exon ensg id mapping"),
	parser.add_argument("--colnames", "-cnames", type = str, default = "", help = "Comma sep string containing column names"),
	parser.add_argument("--colselect", "-csel", type = str, default = "", help = "Comma sep string contianing column names to select"),
	parser.add_argument("out_path", type = str, help = "Path to save outputted csv file"),
	parser.add_argument("--reads_to_genes_path", type = str, help = "Path to save outputted txt file for reads to genes counts")
	return parser.parse_args()

def get_format_data(file, index_cols, regex_sample, colnames, col_to_select):
	sample = re.search(regex_sample, file).group()
	data = pd.read_csv(file, sep = "\s+", index_col = index_cols, names = colnames + [sample]) if colnames != [""] else pd.read_csv(file, sep = "\s+", index_col = index_cols)
	
	if col_to_select != "":
		data = pd.DataFrame(data.loc[:, col_to_select])
		data.rename(columns = {col_to_select : sample}, inplace = True) 
	return data

def merge_data(full_df, to_merge):
	return full_df.merge(to_merge, left_index=True, right_index=True, how='outer')

def add_gene_num(ensg_exon_file, full_df):
	exon_table = pd.read_csv(ensg_exon_file, sep="\s+", index_col=1)
	exon_table = exon_table.loc[~exon_table.index.duplicated(keep="first")]
	full_df.reset_index(inplace=True)
	full_df.set_index("exon_id", inplace=True)
	full_df.insert(0, "gene_num", exon_table)
	return full_df

def add_gene_name(ensg_gene_file, full_df, ensg_col_name):
	gene_names = pd.read_csv(ensg_gene_file, sep="\s+", index_col=0)
	full_df.reset_index(inplace=True)
	full_df.set_index(ensg_col_name, inplace=True)
	full_df.insert(0, "gene_name", gene_names)
	return full_df.reset_index()

def main():
	args = argparser()
	file_list = sorted(args.file_list.split(" "))
	icols = [int(col) for col in args.index_cols.split(",")]
	colnames = args.colnames.split(",")
	colsel = args.colselect.split(",")[0]

	first_one = True
	for file in file_list:
		tmp = get_format_data(file, icols, args.sample_regex, colnames, colsel)
		full_df = tmp if first_one else merge_data(full_df, tmp)
		first_one = False 
	
	full_df.fillna(value = 0, inplace = True)
	full_df = full_df.astype(int)
	
	#Ann file
	full_df = add_gene_num(args.ensg_exon_file, full_df) if args.ensg_exon_file != "" else full_df
	ensg_col_name = "gene_num" if "gene_num" in full_df.reset_index().columns.tolist() else "ensg" #To accnt for dif colnames for gene_num
	full_df = add_gene_name(args.ensg_gene_file, full_df, ensg_col_name = ensg_col_name)
	
	#Reorg columns
	cols = [col for col in full_df.columns.tolist() if col not in ["gene_name", ensg_col_name]]
	cols = ["gene_name", ensg_col_name] + cols

	full_df = full_df[cols]

	full_df.to_csv(args.out_path, index=False)

	if args.ensg_exon_file != "":
		return

	suffix_idx = args.out_path.find('.csv')
	addtional_out_path = args.out_path[:suffix_idx] + '_genes_only' + args.out_path[suffix_idx:] 
	full_df.iloc[:-6].to_csv(addtional_out_path, index = False)

	r2g = full_df.iloc[:-6, 2:].sum(axis = 0)
	r2g.index.rename('sample', inplace = True)
	r2g.rename('reads_to_genes', inplace = True)
	r2g.to_frame().to_csv(args.reads_to_genes_path, sep = '\t')

if __name__ == "__main__":
	main()


