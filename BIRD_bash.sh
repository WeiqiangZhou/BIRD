#!/bin/bash
##Usage: bash BIRD_bash.sh BIRD_directory input_file output_name partition_size wig_mode
if [ $# -ne 5 ]
then echo "Usage: bash BIRD_bash.sh BIRD_directory input_file output_name partition_size wig_mode"
     echo "Example:"
     echo "WIG output: bash BIRD_bash.sh BIRD_v1.0 input.txt output_name 100 1"
     echo "Standard output: bash BIRD_bash.sh BIRD_v1.0 input.txt output_name 100 0"
else
mkdir ./BIRD_tmp
Rscript $1/R_script/split_input.r $2 $4

if [ -f "./BIRD_tmp/BIRD_input_part1" ]
then

if [ $5 -eq "1" ]
then
for file in ./BIRD_tmp/BIRD_input_part*
do
    echo "Processing ${file}"
    BIRD_predict -b $1/model/model_file.bin -i ${file} -o $3 -w
done
else
for file in ./BIRD_tmp/BIRD_input_part*
do
    echo "Processing ${file}"
    part_num=${file##*BIRD_input_}
    BIRD_predict -b $1/model/model_file.bin -i ${file} -o $3.${part_num}
done
fi

fi
rm -r ./BIRD_tmp
fi
