####Tophat
Use human genome hg19 as the reference
```
tophat -p 10 -o /tophat_output_directory /genome_index_directory/hg19 RNA_seq_mate1.fastq RNA_seq_mate2.fastq
```
####Cufflinks
We use gene annotation (GTF) file obtained from GENCODE (Release 19 (GRCh37.p13)) ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

To run Cufflinks:
```
cufflinks -o /cufflinks_output_directory -G gencode.v19.annotation.gtf /tophat_output_directory/accepted_hits.bam
```
