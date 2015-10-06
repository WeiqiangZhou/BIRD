####How to build the prediction model
To build the prediction model, you should have normalized and log2 transformed DNase-seq data and the quantile normalized gene expression data in hand.
Example of the DNase-seq data:
```
chromosome  start    end  AG04449  AG04450   AG10803      AoAF        BJ
      chr1  10200  10399 2.153723 2.485426 0.9619968 2.2021475 0.9770899
      chr1  10400  10599 1.510176 1.910179 1.5340658 1.3807791 2.0641816
      chr1  11200  11399 1.473960 1.173254 1.5340658 0.6274026 1.5197826
      chr1 115600 115799 1.042391 0.000000 0.3623239 0.0000000 0.0000000
      chr1 235600 235799 1.305501 1.845236 0.6930962 2.1039988 1.3927213
```
