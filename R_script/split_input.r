##Usage: Rscript split_input.r input_file partition_size##
input_param <- commandArgs(TRUE)

if(length(input_param)!=2){
print("Please input arguments.",quote=FALSE)
print("example: Rscript split_input.r input_file.txt 100",quote=FALSE)
q(save="no")
}

sep_size <- as.numeric(input_param[2])
data_in <- read.table(file=input_param[1],comment.char="")
sample_size <- dim(data_in)[2] - 1

if(sample_size<=sep_size){
print("Partition size is larger than sample size.",quote=FALSE)
q(save="no")
}

cut_size <- floor(sample_size/sep_size)
sub_idx <- cut(seq(1,sample_size),breaks=cut_size,label=FALSE)

for(i in 1:cut_size){
write.table(data.frame(data_in[,1],data_in[,which(sub_idx==i)+1]),file=paste0("./BIRD_tmp/BIRD_input_part",i),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

q(save="no")
