#!/bin/bash
#$ -cwd
#$ -V
#$ -N anaylzetrait
#$ -l h_data=8G,time=10:00:00
#$ -o logs/
#$ -e logs/
#$ -t 1-45
#$ -m ea

set -e
set -u

. /u/local/Modules/default/init/modules.sh
module load R

array=($(ls -d data2/*.tsv.bgz))
file=${array["$SGE_TASK_ID"]}
bname=`basename $file`
t=$1
f=$2
Rscript scripts/analyze_trait_t_corrected.R $file $t $f > results_t_${t}_corrected_${f}/$bname.txt
