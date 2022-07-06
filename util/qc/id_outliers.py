import argparse
import pandas as pd
import numpy as np
import re
import os
import matplotlib.pyplot as plt
import datetime

def argparser():
    parser = argparse.ArgumentParser(description = "ID Outliers")
    #parser.add_argument("ribo_paths", type = str, help = "Comma separated paths for all ribosomal QC files")
    parser.add_argument("deg_paths", type = str, help = "Comma separated paths for all intron exon ratio QC files")
    parser.add_argument("intron_exon_paths", type = str, help = "Comma separated paths for all RNA degradataion QC files")
    parser.add_argument("reads_to_genes_paths", type = str, help = "Comma separated paths for all sum reads that map to genes QC files")
    parser.add_argument("out_folder", type = str, help = "Path to parent folder within to save output")
    parser.add_argument("--cutoffs_file", type = str, default = '', help = "Path to cutoff file to use")
    parser.add_argument("--outlier_cutoff", "-oc", default = 95, type = int, help = "Specifies cutoff for calling samples outliers. Default = 95")
    parser.add_argument("--sample_regex", "-sr", default = None, type = str, help = "Sample regex to select outlier info for some of samples passed above")
    parser.add_argument("--outliers_to_include", '-out_incl', default = '', type = str, help = "Comma separated list of outlier types to include in final analysis. Accepted values are ribo,dna,deg,r2g")
    return parser.parse_args()

def get_outliers(df, threshold, key):
    df = df.assign(outlier = df.iloc[:,0] > threshold[0])
    outliers = df.loc[df.outlier == True]

    outliers = outliers.assign(cutoff = threshold[0], reason = key)
    col_name = outliers.columns[1]
    outliers.rename(index=str, columns={col_name:"value"}, inplace=True)
    return outliers

def read_qc_file(path):
    return pd.read_csv(path, sep = "\t", index_col = 0) #2 col file where first is sample name

def read_merge_qc_files(paths):
    df_made = 0

    for path in paths:
        df = read_qc_file(path)
        df_full = df if df_made == 0 else pd.concat((df_full, df), axis = 0)
        df_made += 1

    return df_full.round(2)

def get_percentiles(df, p):
    return np.round(np.percentile(df.iloc[:, 0].dropna(), q=[p]), 2)[0]

def plot_dist(df, percentile, figure_name, to_readjust = True):
    to_plot = df.iloc[:,0].copy()
    if to_readjust:
        to_plot.loc[to_plot > 1000] = 30 #Readjust big outliers
    to_plot.dropna(inplace = True)
    plt.figure()
    plt.hist(to_plot, bins=20, color='c', edgecolor='k')
    plt.axvline(percentile, color='k', linestyle='dashed', linewidth=1)
    plt.savefig(figure_name, bbox_inches='tight')

def check_make_folder(path):
    if not (os.path.exists(path)):
        os.makedirs(path)

def save_list(x, path):
    with open(path, "w") as f:
        for item in x:
            f.write("%s\n" % item)
    f.close()
    return

def main():
    args = argparser()
    paths = {#"ribo" : (args.ribo_paths.split(","), 'geq'),
            "r2g" : (args.reads_to_genes_paths.split(","), 'leq'),
            "deg" : (args.deg_paths.split(","), 'geq'),
            "dna" : (args.intron_exon_paths.split(","), 'geq'),
            }


    plot_folder = args.out_folder + "plots/"
    check_make_folder(plot_folder)

    outlier_df_exists = False

    p_row = str("percentile_" + str(args.outlier_cutoff))
    outlier_cutoffs = pd.DataFrame(columns = list(paths.keys()), index = [p_row])
    
    if args.cutoffs_file != '':
        outlier_cutoffs = pd.read_csv(args.cutoffs_file, sep = "\t", index_col = 0)
        print(outlier_cutoffs)
        p_row = outlier_cutoffs.index.to_list()[0]

    is_outlier_keys = []

    for qc_type, files_comparison_tuple in paths.items():
        files, to_compare = files_comparison_tuple
        #Load files, get cutoffs
        qc_df = read_merge_qc_files(files)

        if args.cutoffs_file == '':
            calc_p = get_percentiles(qc_df, args.outlier_cutoff) if to_compare == 'geq' else get_percentiles(qc_df, (100-args.outlier_cutoff))
            outlier_cutoffs.loc[p_row, qc_type] = calc_p

        is_outlier = qc_df.iloc[:, 0] >= outlier_cutoffs.loc[p_row, qc_type] if to_compare == 'geq' else qc_df.iloc[:, 0] <= outlier_cutoffs.loc[p_row, qc_type]
        
        #QC plots
        plot_dist(qc_df, outlier_cutoffs.loc[p_row, qc_type], plot_folder + qc_type + ".png", to_readjust = False if 'r2g' else True)

        #Make df
        outlier_df = qc_df if not outlier_df_exists else outlier_df.join(qc_df)
        is_outlier_keys.append(qc_type + "_outlier")
        outlier_df.insert(outlier_df.shape[1], is_outlier_keys[-1], is_outlier)
        outlier_df_exists = True

    is_outlier_keys_include = is_outlier_keys 
    if args.outliers_to_include != '': 
        outliers_to_include = args.outliers_to_include.split(',') 
        is_outlier_keys_include = [key for key in is_outlier_keys if key.split("_outlier")[0] in outliers_to_include]

    outlier_df.insert(outlier_df.shape[1], "is_outlier", (np.sum(outlier_df.loc[:, is_outlier_keys_include], axis = 1) > 0))
    outlier_df.insert(outlier_df.shape[1], 'outlier_type', 'none')
    for qc_type in paths.keys():
        outlier_df.loc[outlier_df.loc[:, qc_type + "_outlier"], "outlier_type"] = qc_type

    if args.sample_regex is not None:
        idx_matches_regex = [bool(re.search(args.sample_regex, idx))for idx in outlier_df.index.to_list()]
        outlier_df = outlier_df[idx_matches_regex]

    outlier_names = sorted(outlier_df.loc[outlier_df.is_outlier].index.to_list())
    save_list(outlier_names, args.out_folder + "outlier_names.txt")
    outlier_df.to_csv(args.out_folder + "outlier_data.txt", sep = "\t")

    if args.cutoffs_file == '':
        outlier_cutoffs.to_csv(args.out_folder + "outlier_cutoffs" + str(datetime.datetime.now()) + ".txt", sep = "\t")

    return

if __name__ == "__main__":
    main()

    
