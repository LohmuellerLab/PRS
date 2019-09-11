#!/bin/bash
#$ -cwd
#$ -V
#$ -N anaylzetrait
#$ -l h_data=8G,time=10:00:00,highp
#$ -o logs/
#$ -e logs/
#$ -m ea

set -e
set -u

. /u/local/Modules/default/init/modules.sh
module load R

for file in `ls data/*.tsv`;
do
    Rscript scripts/analyze_trait.R $file >> results/private_va.txt
done