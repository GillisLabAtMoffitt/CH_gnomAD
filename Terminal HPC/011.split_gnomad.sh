#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1,pmem=100gb

cd /share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data

zcat gnomad_variant_data.gz |
cut -f1,2-5,8- |
awk '{ print | ("gzip -c > " $1 ".vcf.gz") }'
