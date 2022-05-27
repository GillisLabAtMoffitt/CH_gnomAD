#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1,pmem=100gb

cd /share/lab_gillis/Christelle/cosmic_raw_data/splitted_data

zcat CosmicMutant_w_chr.txt.gz |
awk '{ print | ("gzip -c > " $1 ".vcf.gz") }'
