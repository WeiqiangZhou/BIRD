## BIRD: Big data Regression for predicting DNase I hypersensitivity

### Overview
BIRD is a software to predict DNase I hypersensitivity (DNase-seq signal) based on gene expression data (support both human exon array and RNA-seq data). Using a pre-built model and input gene expression data, BIRD is capable to predict DNase-seq signal genome-wide (~1M genomic loci). BIRD provided two types of outputs: (1) data matrix format or (2) WIG format. Users can easily visualize the predicted DNase-seq signals in UCSC genome browser. 

### News
**The training data for BIRD using DNase-seq and exon array is now available in:
https://github.com/WeiqiangZhou/BIRD-data**

**BIRD is updated to v1.1 which provides more options in the prediction function**
```
Options:
-b   Specify library file. If not sepecified,the program will search for model_file.bin in the current directory.
-i   Specify input file (gene expression obtained from GeneBASE).
-o   Specify output file.
-u   Set upper bound for predicted values (default:14).
-w   Output WIG file for each sample.
-l   Use locus-level model for prediction.
```
Specifically, the predicted values are now bounded from 0 to 14. Users can use the -u option to change the upper bound when using their own prediction model. Users can also use -l option to perform prediction using the locus-level model rather than using the full model. This is useful when you build your own prediction model but you are not sure if the cluster-level model works or not.

To be consistent, **ref_gene.txt** file is renamed as **gene_name.txt**.

### Installation
Currently, BIRD supports Linux/Unix system. Download the latest version of BIRD from: 

https://github.com/WeiqiangZhou/BIRD/releases/download/v1.1/BIRD_v1.1.zip
```
unzip BIRD_v1.1.zip
cd BIRD_v1.1
make
```
### How to use (for exon array)
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
-b   Specify library file. If not sepecified,the program will search for model_file.bin in the current directory.
-i   Specify input file (gene expression obtained from GeneBASE).
-o   Specify output file.
-u   Set upper bound for predicted values (default:14).
-w   Output WIG file for each sample.
-l   Use locus-level model for prediction.                                                                 
```
### For large dataset
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
### Example
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

### How to use (for RNA-seq)
First, use **Tophat** and **Cufflinks** to obtain the gene expression (i.e. FPKM) for the input sample.

To download and install Tophat and Cufflinks, see https://ccb.jhu.edu/software/tophat/tutorial.shtml and http://cole-trapnell-lab.github.io/cufflinks/getting_started/

How to use Tophat and Cufflinks, see https://github.com/WeiqiangZhou/BIRD/blob/master/RNAseq_instruction.md

After obtaining the gene expression data **genes.fpkm_tracking**, used the bash script **match_gene.sh** to prepare the input file for BIRD.
```
bash path_to_BIRD/match_gene.sh genes.fpkm_tracking path_to_BIRD/model/gene_name.txt genes.fpkm_tracking.match
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

### For RNA-seq with FPKM data matrix
If you have a data matrix containing the gene names and FPKM for multiple samples (see **FPKM_data_matrix.txt** in the **example** folder as an example), use the R script **match_input_matrix.r** to prepare the input data before running BIRD prediction.

For example:
```
Rscript --vanilla path_to_BIRD/R_script/match_input_matrix.r path_to_BIRD/model/gene_name.txt FPKM_data_matrix.txt FPKM_data_matrix_match.txt

path_to_BIRD/BIRD_predict -b path_to_BIRD/model/RNAseq_model_file.bin -i FPKM_data_matrix_match.txt -o output_file.txt
```

### How to build the prediction model
The BIRD software package contains the pre-built prediction model for both exon array and RNA-seq data.

The training data for BIRD using DNase-seq and exon array is now available in:
https://github.com/WeiqiangZhou/BIRD-data

If you would like to know how to build a prediction model, see https://github.com/WeiqiangZhou/BIRD/blob/master/build_prediction_model.md for details.
### Note:
Change **path_to_BIRD** to the path where you install BIRD.

Genomic information in the current version of BIRD is based on **human genome hg19**.

### PDDB: Predicted DNase I hypersensitivity database
A database of the predicted DNase I hypersensitivity for 2,000 GEO exon array samples is available at http://jilab.biostat.jhsph.edu/~bsherwo2/bird/index.php

### Interpretation of the prediction models
The predictor genes for each genomic locus and each DHS cluster, motif enrichment analysis results, GO analysis results, and the active cell types for each DHS cluster in the prediction models are provided as an online resource which is available at https://zhiji.shinyapps.io/CABS/.

Users can input a list of interested genomic regions and explore the overlapped DHSs from the prediction models. 

### Contact
Weiqiang Zhou: wzhou14@jhu.edu

### References
Zhou W, Sherwood B, Ji Z, Xue Y, Du F, Bai J, Ying M, Ji H. Genome-wide Prediction of DNase I Hypersensitivity Using Gene Expression. _Nature Communications_ 8, 1038 (2017). ([open access](https://www.nature.com/articles/s41467-017-01188-x))

Zhou W, Ji Z, Ji H. Global Prediction of Chromatin Accessibility Using RNA-seq from Small Number of Cells. bioRxiv doi: 10.1101/035816 (2016). ([preprint](https://www.biorxiv.org/content/early/2016/01/03/035816))
