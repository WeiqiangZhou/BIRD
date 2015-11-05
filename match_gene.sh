#!/bin/bash
##Usage: bash match_gene.sh input_file ref_gene_file output_file
##In input_file, column 4 is gene_id and column 10 is FPKM.

awk -v OFS='\t' 'NR==FNR {h[$4] = $10; next} {if(length(h[$1])!=0) print $1,h[$1];else print $1,"0";}' $1 $2 > $3.temp
echo gene_id$'\t'FPKM | cat - $3.temp > $3
rm $3.temp
