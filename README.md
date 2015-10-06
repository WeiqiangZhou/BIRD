###BIRD: Big data Regression for predicting DNase I hypersensitivity

####Overview
BIRD is a software to predict DNase I hypersensitivity (DNase-seq signal) based on gene expression data (support both human exon array and RNA-seq data). Using a pre-built model and input gene expression data, BIRD is capable to predict DNase-seq signal genome-wide (~1M genomic loci). BIRD provided two types of outputs: (1) data matrix format or (2) WIG format. Users can easily visualize the predicted DNase-seq signals in UCSC genome browser. 

####Installation
Currently, BIRD supports Linux/Unix system. Download the latest version of BIRD from: 

https://github.com/WeiqiangZhou/BIRD/releases/download/v1.0.2/BIRD_v1.0.2.zip
```
unzip BIRD_v1.0.2.zip
cd BIRD_v1.0.2
make
```
####How to use (for exon array)
BIRD accepts gene expression output file from **GeneBASE**.
If you have the raw exon array data (CEL file), use GeneBASE to get the gene expression. 

To download and install GeneBASE, see http://web.stanford.edu/group/wonglab/GeneBASE/

How to use GeneBASE, see https://github.com/WeiqiangZhou/BIRD/blob/master/exonarray_instruction.md

After running GeneBASE, you will get the gene expression data file (e.g. input_file.txt).

To get data matrix format output, run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i input_file.txt -o output_file.txt
```
To get WIG format output, run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i input_file.txt -o output_name -w
```
In this mode, BIRD will generate WIG file for each sample with prefix "output_name." followed by the column name in the input_file.txt.

WIG file can be visualized in UCSC genome browser by adding custom tracks:

http://genome.ucsc.edu/cgi-bin/hgGateway

For help information, run:
```
path_to_BIRD/BIRD_predict -h
```
```
Usage:                                                                                                      
Standard output: BIRD_predict -b model_file.bin -i input_file.txt -o output_file.txt                        
Standard output will save a matrix contained all predited value in log scale (log2(x+1) transformed).       
WIG output: BIRD_predict -b model_file.bin -i input_file.txt -o output_name -w                              
WIG output will save each sample as a WIG file.                                                             
Options:                                                                                                    
-b   Specify library file. If not sepecified,the program will use model_file.bin in the current directory.  
-i   Specify input file (gene expression obtained from GeneBASE).                                           
-o   Specify output file.                                                                                   
-w   Output WIG file for each sample.                                                                       
```
####For large dataset
It could be memory intensive to run BIRD on dataset with sample size larger than 100 (memory usage > 1G). In such cases, it is recommended to run BIRD via the bash script **BIRD_bash.sh**. It is required that **R** is installed and make sure the **Rscript** command is executable. The bash script will partition the input file into N files according to partition_size (e.g. 100, the number of samples in a partitioned input file). 

To get data matrix format output, run:
```
bash path_to_BIRD/BIRD_bash.sh path_to_BIRD input_file.txt output_name partition_size 0
```
The output files should be **output_name.part1, ..., output_name.partN**.

To get WIG format output, run:
```
bash path_to_BIRD/BIRD_bash.sh path_to_BIRD input_file.txt output_name partition_size 1
```
####Example
Run BIRD:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i path_to_BIRD/example/Exon_K562_lab.txt -o K562_DNase.txt
```
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i path_to_BIRD/example/Exon_K562_lab.txt -o K562_DNase -w
```
You should get three output files: K562_DNase.txt, K562_DNase.Duke.wig, K562_DNase.UW.wig

Run BIRD_bash.sh:
```
bash path_to_BIRD/BIRD_bash.sh path_to_BIRD path_to_BIRD/example/Exon_K562_lab.txt K562_DNase.txt 1 0
```
```
bash path_to_BIRD/BIRD_bash.sh path_to_BIRD path_to_BIRD/example/Exon_K562_lab.txt K562_DNase 1 1
```
You should get four output files: K562_DNase.txt.part1, K562_DNase.txt.part2, K562_DNase.Duke.wig, K562_DNase.UW.wig

####How to use (for RNA-seq)
First, use **Tophat** and **Cufflinks** to obtain the gene expression (i.e. FPKM) for the input sample.

To download and install Tophat and Cufflinks, see https://ccb.jhu.edu/software/tophat/tutorial.shtml and http://cole-trapnell-lab.github.io/cufflinks/getting_started/

How to use Tophat and Cufflinks, see https://github.com/WeiqiangZhou/BIRD/blob/master/RNAseq_instruction.md

After obtaining the gene expression data **genes.fpkm_tracking**, used the R script **match_input.r** to prepare the input file for BIRD.
```
Rscript path_to_BIRD/R_script/match_input.r genes.fpkm_tracking path_to_BIRD/model/ref_gene.txt genes.fpkm_tracking.match
```
Then run the **BIRD_predict** program for prediction (use model file **RNAseq_model_file.bin**):

To get data matrix format output:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/RNAseq_model_file.bin -i genes.fpkm_tracking.match -o output_file.txt
```
To get WIG format output, run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/RNAseq_model_file.bin -i genes.fpkm_tracking.match -o output_name -w
```
####How to build the prediction model
The BIRD software package contains the pre-built prediction model for both exon array and RNA-seq data. If you would like to know how to build a prediction model, see https://github.com/WeiqiangZhou/BIRD/blob/master/build_prediction_model.md for details.
####Note:
Change **path_to_BIRD** to the path where you install BIRD.

Genomic information in the current version of BIRD is based on **human genome hg19**.

####Contact
wzhou14@jhu.edu

