###BIRD: Big data Regression for predicting DNase I hypersensitivity

####Overview
BIRD is a software to predict DNase I hypersensitivity (DNase-seq signal) based on gene expression data (currently support exon array data). Using a pre-built model and input gene expression data, BIRD is capable to predict DNase-seq signal genome-wide (~1M genomic loci). BIRD provided two types of outputs: (1) data matrix format or (2) WIG format. Users can easily visualize the predicted DNase-seq signals in UCSC genome browser. 

####Installation
Currently, BIRD supports Linux/Unix system. Download the release version of BIRD from: 

https://github.com/WeiqiangZhou/BIRD/releases/download/v1.0/BIRD_v1.0.zip
```
unzip BIRD_v1.0.zip
cd BIRD_v1.0
make
```
####How to use
BIRD accepts gene expression output file from geneBASE.
If you have the raw exon array data (CEL file), use geneBASE to get the gene expression. For detail, see http://web.stanford.edu/group/wonglab/GeneBASE/

After running geneBASE, you will get the gene expression data file (e.g. input_file.txt).

To get data matrix format output, run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i input_file.txt -o output_file.txt
```
To get WIG format output, run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i sample_expr.txt -o output_name -w
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
-i   Specify input file (gene expression obtained from geneBASE).                                           
-o   Specify output file.                                                                                   
-w   Output WIG file for each sample.                                                                       
```
Note: 

Change path_to_BIRD to the path where you install BIRD.

####Example
Run:
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i path_to_BIRD/example/Exon_K562_lab.txt -o K562_DNase.txt
```
```
path_to_BIRD/BIRD_predict -b path_to_BIRD/model/model_file.bin -i path_to_BIRD/example/Exon_K562_lab.txt -o K562_DNase -w
```
You should get three output files: K562_DNase.txt, K562_DNase.Duke.wig, K562_DNase.UW.wig

####Contact
wzhou14@jhu.edu

