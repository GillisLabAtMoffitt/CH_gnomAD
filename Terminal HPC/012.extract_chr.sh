#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1,pmem=100gb

cd /share/lab_gillis/Christelle/cosmic_raw_data/splitted_data

zcat CosmicMutantExport.tsv.gz | 
awk -v FS='\t' -v OFS='\t' -f 012.extract_chr.awk | 
awk '{if(NR==1){print $0} ;if(NR>1) {$0="chr"$0; print $0}}' |
gzip -c > CosmicMutant_w_chr.txt.gz
