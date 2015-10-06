#######################################################################################################
#                                        README
# This file contains the R code for performing clustering to the DNase-seq data using bigKmeans.
# The bigKmeans program should be first installed. See https://github.com/fangdu64/BDT for details.
#######################################################################################################

####functions####
##standardize data
standardize_row_train <- function (Data_train)
{
    train_mean_sd <- matrix(data = NA, nrow = dim(Data_train)[1],
        ncol = 2)
    Data_train_sd <- matrix(data = NA, nrow = dim(Data_train)[1],
        ncol = dim(Data_train)[2])
    for (i in 1:dim(Data_train)[1]) {
        train_mean_sd[i, 1] <- mean(Data_train[i, ])
        train_mean_sd[i, 2] <- sd(Data_train[i, ])
        if (train_mean_sd[i, 2] == 0) {
            Data_train_sd[i, ] <- rep(0, dim(Data_train)[2])
        }
        else {
            Data_train_sd[i, ] <- (Data_train[i, ] - train_mean_sd[i,
                1])/train_mean_sd[i, 2]
        }
    }
    result <- list(mean_sd = train_mean_sd, train = Data_train_sd)
}

##convert clustering results from binary to text##
readBinary <- function(filename,vec_length){
   to.read = file(filename, "rb")
   input <- readBin(to.read, integer(), n=vec_length)
}

##generate script to run BigKmeans
##see https://github.com/fangdu64/BDT/blob/master/doc/bigKmeans.md for details of how to use bigKmeans
kmppScript <- function(data,k,thread_num=10,dist_type="Euclidean",max_iter=100,solveThresh=0.0001,
                       tmp_path="/tmp/tmp_data.txt", out_path="/tmp"){
     write.table(data,tmp_path,row.names=FALSE,col.names=FALSE)
     nCol <- ncol(data)
     nRow <- nrow(data)
     sink("kmpp.sh")
     cat("/bdt-0.1.1/bigKmeans \\") #path to the executable bigKmeans program
     cat("\n")
     cat(paste("--data-input text-mat@", tmp_path, " \\", sep="")) #input data path
     cat("\n")
     cat(paste("--data-nrow", nRow, "\\")) #number of rows in the data
     cat("\n")
     cat(paste("--data-ncol", nCol, "\\")) #number of columns in the data
     cat("\n")
     cat(paste("--k", k, "\\")) #number of clusters
     cat("\n")
     cat(paste("--out", out_path, "\\")) #output path
     cat("\n")
     cat(paste("--thread-num", thread_num, "\\")) #number of threads used for parallell process
     cat("\n")
     cat(paste("--dist-type", dist_type, "\\")) #type of distance measurement 
     cat("\n")
     cat(paste("--max-iter", max_iter, "\\")) #maximum number of iterations allowed
     cat("\n")
     cat(paste("--min-expchg", solveThresh)) #the minimum change of fraction of varaince explained
     sink()
     system(paste("chmod +x kmpp.sh"))
     system("./kmpp.sh")
  
     cluster_results <- readBinary(paste(out_path,"/run/3-run-kmeans/cluster_assignments.bfv",sep=""),nRow) #convert clustering results from binary file to text file
     cluster_results + 1 #start with group 1 instead of 0 
}

##main##
load("DNase_data.rda")  ##first three columns contain the genomic locus information
DNase_processed <- standardize_row_train(as.matrix(DNase_data[,-c(1:3)])) #standardize DNase-seq data
DNase_train_sd <- DNase_processed$train

#perform clustering to the DNase-seq data for K=1000, 2000 and 5000
clusters <- kmppScript(DNase_train_sd,1000)
write.table(clusters,file="DH_cluster_1000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

clusters <- kmppScript(DNase_train_sd,2000)
write.table(clusters,file="DH_cluster_2000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

clusters <- kmppScript(DNase_train_sd,5000)
write.table(clusters,file="DH_cluster_5000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

q(save="no")

