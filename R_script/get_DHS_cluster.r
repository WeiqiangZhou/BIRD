#######################################################################################################
#
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


##cluster function##
readBinary <- function(filename,vec_length){
   to.read = file(filename, "rb")
   input <- readBin(to.read, integer(), n=vec_length)
}

kmppScript <- function(data,k,thread_num=10,dist_type="Euclidean",max_iter=100,solveThresh=0.0001,
                       tmp_path="/dcs01/gcode/wzhou14/TF_project/paper_prepare/new_result/40cell_train/replicate/r4/tmp/tmp_data.txt",
                       out_path="/dcs01/gcode/wzhou14/TF_project/paper_prepare/new_result/40cell_train/replicate/r4/tmp"){
     write.table(data,tmp_path,row.names=FALSE,col.names=FALSE)
     nCol <- ncol(data)
     nRow <- nrow(data)
     sink("kmpp.sh")
     #cat("/dcs01/gcode/fdu1/BDT/pipeline-kmeans++.py \\")
     cat("/dcs01/gcode/fdu1/install/gcc-4.4.7/bdt-0.1.1/bigKmeans \\")
     cat("\n")
     cat(paste("--data-input text-mat@", tmp_path, " \\", sep=""))
     cat("\n")
     cat(paste("--data-nrow", nRow, "\\"))
     cat("\n")
     cat(paste("--data-ncol", nCol, "\\"))
     cat("\n")
     cat(paste("--k", k, "\\"))
     cat("\n")
     cat(paste("--out", out_path, "\\"))
     cat("\n")
     cat(paste("--thread-num", thread_num, "\\"))
     cat("\n")
     cat(paste("--dist-type", dist_type, "\\"))
     cat("\n")
     cat(paste("--max-iter", max_iter, "\\"))
     cat("\n")
     cat(paste("--min-expchg", solveThresh))
     sink()
     system(paste("chmod +x kmpp.sh"))
     system("./kmpp.sh")
  
     cluster_results <- readBinary(paste(out_path,"/run/3-run-kmeans/cluster_assignments.bfv",sep=""),nRow)
     #print(cluster_results)
     cluster_results + 1 #start with group 1 instead of 0 
}

##main##
load("DNase_data_r.rda")  ##first three columns contain the genomic locus information
DNase_processed <- standardize_row_train(as.matrix(DNase_data[,-c(1:3)]))
DNase_train_sd <- DNase_processed$train

clusters <- kmppScript(DNase_train_sd,1000)
write.table(clusters,file="DH_cluster_1000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

clusters <- kmppScript(DNase_train_sd,2000)
write.table(clusters,file="DH_cluster_2000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

clusters <- kmppScript(DNase_train_sd,5000)
write.table(clusters,file="DH_cluster_5000.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

q(save="no")

