#!/bin/bash
##Usage: bash match_gene.sh input_file ref_gene_file output_file

awk -v OFS='\t' 'NR==FNR {h[$4] = $10; next} {print $1,h[$1]}' $1 $2 > $3.temp
echo gene_id$'\t'FPKM | cat - $3.temp > $3
rm $3.temp
