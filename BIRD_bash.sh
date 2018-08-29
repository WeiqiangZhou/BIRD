#!/bin/bash
##Usage: bash BIRD_bash.sh BIRD_directory model_file input_file output_name partition_size wig_mode
if [ $# -ne 6 ]
then echo "Usage: bash BIRD_bash.sh BIRD_directory model_file input_file output_name partition_size wig_mode"
     echo "Example:"
     echo "WIG output: bash BIRD_bash.sh path_to_BIRD path_to_BIRD/model/model_file.bin input.txt output_name 100 1"
     echo "Standard output: bash BIRD_bash.sh path_to_BIRD path_to_BIRD/model/model_file.bin input.txt output_name 100 0"
else

if [ ! -d "BIRD_tmp" ]
then mkdir ./BIRD_tmp
fi

Rscript $1/R_script/split_input.r $3 $5

if [ -f "./BIRD_tmp/BIRD_input_part1" ]
then

if [ $6 -eq "1" ]
then
for file in ./BIRD_tmp/BIRD_input_part*
do
    echo "Processing ${file}"
    $1/BIRD_predict -b $2 -i ${file} -o $4 -w
done
else
for file in ./BIRD_tmp/BIRD_input_part*
do
    echo "Processing ${file}"
    part_num=${file##*BIRD_input_}
    $1/BIRD_predict -b $2 -i ${file} -o $4_${part_num}.txt
done
fi

fi
rm -r ./BIRD_tmp
fi
