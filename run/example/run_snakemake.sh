#!/bin/bash
ROOT=/oak/stanford/groups/quake/miramou #CHANGE TO MATCH FOR YOU
LOG_DIR=$ROOT/example/analysis #CHANGE TO MATCH DATASET
FLOW_DIR=$ROOT/resources/comp_pipelines/cfRNA/cfRNA_pipeline #CHANGE TO MATCH FOR YOU
RUN_DIR=$FLOW_DIR/run/example #CHANGE TO MATCH DATASET
SNAKEFILE=$FLOW_DIR/Snakefile_modular.py
CONFIGFILE=$RUN_DIR/config.yaml
CLUSTERFILE=$FLOW_DIR/util/config/cluster_config.json
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$RUN_DIR/log/pipeline.$DATETIME.log

CALL_SMK () {
  snakemake all --snakefile $SNAKEFILE --directory=$FLOW_DIR --use-conda --configfile $CONFIGFILE --cluster-config $CLUSTERFILE --cluster "sbatch --ntasks=1 --job-name={cluster.name} --cpus-per-task={threads} --partition={cluster.partition} --mem={resources.mem_mb} -t {cluster.time} -o $LOG_DIR/{cluster.output} -e $LOG_DIR/{cluster.error}" --keep-target-files -j 100 -w 100 -p -r --rerun-incomplete --restart-times 1 $1
}

if [ $# -eq 0 ]
  then
    # dry run 
	CALL_SMK -n
elif [ $1 == "unlock" ]
  then
    # unlock
    CALL_SMK --unlock

elif [ $1 == "snakemake" ]
  then
    #Start snakemake and save output to log
    CALL_SMK "" 2>&1 | tee $LOGFILE

elif [ $1 == 'sbatch' ]
  then
    sbatch --ntasks=1 --job-name='snkamaka' --cpus-per-task=1 --time=1-0:00:00 -p quake,normal,owners --mem=2G -o $LOGFILE -e $LOGFILE $RUN_DIR/run_snakemake.sh snakemake

else
    echo 'error code 5402901u47 - sad snakamaka'
fi

echo snakamaka log is $LOGFILE
