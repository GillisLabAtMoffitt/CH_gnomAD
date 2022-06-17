#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1,pmem=100gb
module load R/4.0.2
cd /share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/age_comparison_pval
Rscript age_comparison_pval.R
