include: "util/config/environ.Snakefile"

shell.prefix("set +euo pipefail;") #Force script to fail for debugging purposes

# From run config file, define folders
reads_folder = root_dir + config["reads"]
reads_folder = reads_folder[:-1] if reads_folder[-1] == "/" else reads_folder #Key assumption in rules is that this path does not ends with "/"

data_folder = root_dir + config["data"]
data_folder = data_folder + "/" if data_folder[-1] != "/" else data_folder #Key assumption in rules is that this path ends with "/"

# From run config file, define global vars
merge_mult_lanes = int(config["merge_mult_lanes"])
seq_type = config["seq_type"]

is_two_pass = int(config["is_two_pass"])

stranded_type = config["stranded_type"]
adapter_type = config["adapters"]
get_reads_from_sample_sheet = bool(int(config["get_reads_from_sample_sheet"]))

# Get samples and reads
samples = pd.read_csv(root_dir + config["samples"], index_col=0)
reads = get_read_names()

wildcard_constraints:
	sample=config["wildcard_constraints"]

# Common rules
ruleorder: star>sort_by_coord
localrules: all, make_dirs #, merge_htseq, merge_ribo, merge_dna, merge_deg

rule all:
	input:
		get_output_files

# Target rules
rule_dir = "rules/"

include: rule_dir + "common.smk"
include: rule_dir + "trim.smk"
include: rule_dir + "align.smk"
include: rule_dir + "dedup_sort.smk"
include: rule_dir + "count.smk"
include: rule_dir + "qc.smk"
include: rule_dir + "aggregate.smk"
