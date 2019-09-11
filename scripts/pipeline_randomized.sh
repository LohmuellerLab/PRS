#!/bin/bash
#$ -cwd
#$ -V
#$ -N anaylzetrait
#$ -l h_data=32G,time=23:00:00
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
python scripts/randomized_va.py -w 50000 -m 10000 -a 0.05 -f $file > results3/$bname.50000-10000-0.05.txt
python scripts/randomized_va.py -w 5000 -m 10000 -a 0.05 -f $file > results3/$bname.5000-10000-0.05.txt
python scripts/randomized_va.py -w 500 -m 10000 -a 0.05 -f $file > results3/$bname.500-10000-0.05.txt

python scripts/randomized_va.py -w 50000 -m 10000 -a 0.1 -f $file > results3/$bname.50000-10000-0.1.txt
python scripts/randomized_va.py -w 5000 -m 10000 -a 0.1 -f $file > results3/$bname.5000-10000-0.1.txt
python scripts/randomized_va.py -w 500 -m 10000 -a 0.1 -f $file > results3/$bname.500-10000-0.1.txt