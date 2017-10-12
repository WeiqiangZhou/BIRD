##match the gene order in a data matrix
##Usage: Rscript --vanilla match_input_matrix.r PATH_TO_BIRD/model/ref_gene.txt input_file.txt output_file.txt

args = commandArgs(trailingOnly=TRUE)

data_in <- as.matrix(read.table(file=args[2], header=TRUE, row.names=1, check.names=FALSE))
ref_gene <- read.table(file=args[1], stringsAsFactors=FALSE)

ref_gene_re <- do.call(rbind,strsplit(ref_gene[,1],"\\."))[,1]
input_gene_re <- do.call(rbind,strsplit(row.names(data_in),"\\."))[,1]

match_idx <- match(input_gene_re,ref_gene_re)
data_out <- matrix(data=0, ncol=ncol(data_in), nrow=length(ref_gene_re))
data_out[match_idx[which(match_idx != "NA")],] <- data_in[which(match_idx != "NA"),]
colnames(data_out) <- colnames(data_in)

write.table(data.frame(gene_id=ref_gene[,1], data_out), file=args[3], sep="\t", row.names=FALSE, quote=FALSE)
