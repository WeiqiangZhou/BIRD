### How to build the prediction model
To build the prediction model, you should have normalized and log2 transformed DNase-seq data and the log2 transformed and quantile normalized gene expression data in hand.

Example of the DNase-seq data (the first three columns are the genomic location of the DHSs):
```
chromosome  start    end  AG04449  AG04450   AG10803      AoAF        BJ
      chr1  10200  10399 2.153723 2.485426 0.9619968 2.2021475 0.9770899
      chr1  10400  10599 1.510176 1.910179 1.5340658 1.3807791 2.0641816
      chr1  11200  11399 1.473960 1.173254 1.5340658 0.6274026 1.5197826
      chr1 115600 115799 1.042391 0.000000 0.3623239 0.0000000 0.0000000
      chr1 235600 235799 1.305501 1.845236 0.6930962 2.1039988 1.3927213
```

Example of the gene expression data (the first column is the ids for each gene, e.g., transcript cluster ids from exon array or gene names, ensembl gene ids from RNA-seq):
```
Gene_id   AG04449   AG04450   AG10803      AoAF        BJ       Hac
2359453  6.289403  6.908867  6.438640  7.766761  5.423756  7.028196
2359470  5.683281  4.825203  5.211613  4.692368  5.240313  5.237280
3145953 12.405847 12.393643 12.415091 12.603896 12.590820 12.393643
3539201  5.400856  5.750283  5.157050  3.324667  5.982363  4.270757
2359483  8.043600  7.526427  7.528059  7.525314  7.805880  7.714811
```

You should first install **BigKmeans** and **BIRD**. See https://github.com/fangdu64/BDT for details of how to install BigKmeans. See https://github.com/WeiqiangZhou/BIRD for details of how to install BIRD.

Then, apply the following steps:
#### Step 1:
Use the R code in https://github.com/WeiqiangZhou/BIRD/blob/master/R_script/get_param.r to obtain the optimal number of gene clusters K and number of predictors N.

#### Step 2:
Use the R code in https://github.com/WeiqiangZhou/BIRD/blob/master/R_script/get_DHS_cluster.r to obtain the clustering results of the DNase-seq data.

#### Step 3:
Use the R code in https://github.com/WeiqiangZhou/BIRD/blob/master/R_script/get_model_data.r to obtain the required files to build the BIRD model.

#### Step 4:
Run the following command to compile the BIRD model file:
```
BIRD_build_library -i par_file.txt -o model_file.bin
```
Now, you should get the prediction model file **model_file.bin**.
