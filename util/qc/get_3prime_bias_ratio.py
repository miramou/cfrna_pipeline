import sys
import numpy as np
import pandas as pd
from scipy import stats

sample_sheet = sys.argv[1]
csv_path = sys.argv[2]

sample_anno = pd.read_csv(sample_sheet, sep="\t", index_col=0, header=None, names=["exon_id", "gene_id", "exon_num", "exon_count"])
exon_sum = sample_anno.groupby(["gene_id"]).sum().exon_count
gene_idx = exon_sum > 0 #genes with counts > 0

if np.sum(gene_idx) > 0:
	ratio_3prime = sample_anno.groupby(["gene_id"]).max()[gene_idx].astype(float) #the 3' most exon is that with the highest exon_num
	ratio_3prime.loc[:, "exon_count"] /= exon_sum[gene_idx] #number of reads for 3' exon / total reads for gene = ratio

	frac_all_3prime = np.sum(ratio_3prime.loc[:, "exon_count"] == 1.0)
	frac_all_3prime /= ratio_3prime.shape[0]

	out = pd.DataFrame(list(set(ens_id.loc[:, "gene_id"])), columns=["gene_id"])
	out = out.set_index("gene_id")
	out = out.join(ratio_3prime.loc[:, "exon_count"], on="gene_id")
	out = out.append(pd.Series([frac_all_3prime], index=out.columns, name="frac_3prime"))
	out = out.append(pd.Series([np.mean(ratio_3prime.loc[:, "exon_count"])], index=out.columns, name="mean_3prime"))
	out = out.append(pd.Series([np.median(ratio_3prime.loc[:, "exon_count"])], index=out.columns, name="median_3prime"))
	out = out.append(pd.Series([stats.mode(ratio_3prime.loc[:, "exon_count"])[0][0]], index=out.columns, name="mode_3prime"))
	out = out.append(pd.Series([np.std(ratio_3prime.loc[:, "exon_count"])], index=out.columns, name="std_3prime"))

	out.to_csv(csv_path, sep="\t")
else:
	out = pd.DataFrame()
	out.to_csv(csv_path, sep="\t")

