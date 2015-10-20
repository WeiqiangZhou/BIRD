##Usage Rscript match_input.r input_file ref_gene_file output_file
input_param <- commandArgs(TRUE)

data_in <- read.table(file=input_param[1],header=TRUE)
input_gene_id <- data_in[,"gene_id"]

ref_gene_id <- as.vector(unlist(read.table(file=input_param[2])))

data_out <- rep(0,nrow=length(ref_gene_id))

for(i in 1:length(ref_gene_id)){
if(length(which(input_gene_id==ref_gene_id[i]))==1){
data_out[i] <- data_in[which(input_gene_id==ref_gene_id[i]),"FPKM"]
}
}

data_out <- log2(data_out+1)

test_data <- data.frame(Gene_name=ref_gene_id,value=data_out)
write.table(test_data,file=input_param[3],row.names=FALSE,sep="\t",quote=FALSE)
