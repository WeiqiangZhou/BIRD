#############################################################################################################################
#                                                    README
# This file contains the R code for generating the BIRD prediction model.
# Results for clustering DHSs into 1000, 2000 and 5000 clusters should be provided. 
# See https://github.com/fangdu64/BDT for details of how to install and use BigKmeans.
# After generating the required files, use [BIRD_build_library -i par_file.txt -o model_file.bin] to compile the model file 
#############################################################################################################################

####functions####
##get regression coefficient and predictors
get_regression_info <- function(DNase_train,Exon_train_mean,top_n){
y_cor <- cor(DNase_train,t(Exon_train_mean))
y_cor[is.na(y_cor)] <- 0

if(top_n > 1){
max_idx <- sort(y_cor, decreasing=TRUE, index.return=TRUE)$ix[1:top_n]
data_train <- data.frame(y=DNase_train, t(Exon_train_mean[max_idx,]))
fit <- lm(y~.,data_train)
}else{
max_idx <- which(y_cor == max(y_cor))[1]
data_train <- data.frame(y=DNase_train, x=Exon_train_mean[max_idx,])
fit <- lm(y~x,data_train)
}

return(list(lm_fit=fit,predictor=max_idx))
}

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

##get cluster mean##
cluster_data_mean <- function(Exon_data,cluster){
  
  cluster_unique <- unique(cluster)
  Exon_mean <- matrix(data=NA, nrow=length(cluster_unique), ncol=dim(Exon_data)[2])
  
  
  for(j in 1:length(cluster_unique)){
    cluster_idx <- which(cluster==j)
    Exon_cluster <- Exon_data[cluster_idx,]
    Exon_cluster <- as.matrix(Exon_cluster)
    
    for (i in 1:dim(Exon_data)[2]){
      if(length(Exon_cluster)==dim(Exon_data)[2]){
        Exon_mean[j,i] <- Exon_cluster[i]
      }
      else{
        Exon_mean[j,i] <- mean(Exon_cluster[,i])
      }
    }
  }
  
  return(Exon_mean)
}

##get data quantile##
compute_quantile <- function(train_exon){
  data_sorted <- train_exon
  for(i in 1:dim(train_exon)[2]){
    data_sorted[,i] <- sort(train_exon[,i],decreasing=TRUE)
  }
  data_quantile <- apply(data_sorted,1,mean)
  return(data_quantile)
}

##main script##
Num_cluster <- 1500   ##different number of cluster can be set here
Num_predictor <- 7  ##different number of cluster can be set here
Bin_size <- 200  ##bin size depends on the preprocess of DNase data, we used 200 in BIRD

##input data are Exon_data_sample.txt and DNase_data_sample.rda
Exon_train <- as.matrix(readRDS("Exon_data_57_cells.rds"))  ##read quantile normalized gene expression data
Exon_quantile <- compute_quantile(Exon_train) #calculate the quantiles of the gene expression data
Exon_processed <- standardize_row_train(Exon_train) #standardize gene expression data
Exon_mean_sd <- Exon_processed$mean_sd
Exon_train_sd <- Exon_processed$train
set.seed(2015)
Exon_cluster <- kmeans(Exon_train_sd, centers=Num_cluster, nstart=10, iter.max=50)    ##cluster gene expression

write.table(Exon_mean_sd[,1],file="gene_mean.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Exon_mean_sd[,2],file="gene_sd.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Exon_quantile,file="gene_quantile.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(row.names(Exon_train),file="gene_name.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Exon_cluster$cluster,file="cluster_idx.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

Exon_train_mean <- cluster_data_mean(Exon_train_sd,Exon_cluster$cluster) #calculated average gene expression within each cluster

##Locus-level prediction##
DNase_data <- readRDS("DNase_data_57_cells.rds")  ##read DNase-seq data, first three columns contain the genomic locus information
DNase_processed <- standardize_row_train(as.matrix(DNase_data[,-c(1:3)])) #standardize DNase-seq data
DNase_mean_sd <- DNase_processed$mean_sd
DNase_train_sd <- DNase_processed$train

write.table(DNase_data[,1:3],file="genomic_loci.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
write.table(DNase_mean_sd[,1],file="DNase_mean.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
write.table(DNase_mean_sd[,2],file="DNase_sd.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

coef_all <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=Num_predictor+1)
predictor_idx <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=Num_predictor)

for (i in 1:dim(DNase_train_sd)[1]){
regress_result <- get_regression_info(DNase_train_sd[i,],Exon_train_mean,Num_predictor) #build regression model for each DHS
coef_all[i,] <- coefficients(regress_result$lm_fit) #get regression coefficients
predictor_idx[i,] <- regress_result$predictor #get predictor index
}

write.table(coef_all[,-1],file="regress_coef.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE) ##b0 is zero and excluded
write.table(predictor_idx,file="regress_predictor.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

##DH-cluster-level prediction##
DH_cluster <- list()
DH_cluster_data <- list()
DH_cluster_num <- c(1000,2000,5000)
##Input DH cluster files DH_cluster_1000.txt, DH_cluster_2000.txt and DH_cluster_5000.txt should be included in the same folder.##
for(j in 1:length(DH_cluster_num)){
DH_cluster[[j]] <- as.vector(unlist(read.table(file=paste0("DH_cluster_",DH_cluster_num[j],".txt"))))
DH_cluster_data[[j]] <- cluster_data_mean(DNase_train_sd,DH_cluster[[j]])
}

Dis_matrix <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=length(DH_cluster_num))
for(j in 1:length(DH_cluster_num)){
for(i in 1:dim(DNase_train_sd)[1]){
Dis_matrix[i,j] <- cor(DNase_train_sd[i,],DH_cluster_data[[j]][DH_cluster[[j]][i],]) #calculate the weights for model averaging
}
}
Dis_matrix[which(Dis_matrix<0)] <- 0

coef_cluster <- list()
predictor_idx_cluster <- list()

for(j in 1:length(DH_cluster_num)){
coef_cluster[[j]] <- matrix(data=NA,nrow=dim(DH_cluster_data[[j]])[1],ncol=Num_predictor+1)
predictor_idx_cluster[[j]] <- matrix(data=NA,nrow=dim(DH_cluster_data[[j]])[1],ncol=Num_predictor)
for (i in 1:dim(DH_cluster_data[[j]])[1]){
regress_result <- get_regression_info(DH_cluster_data[[j]][i,],Exon_train_mean,Num_predictor)
coef_cluster[[j]][i,] <- coefficients(regress_result$lm_fit)
predictor_idx_cluster[[j]][i,] <- regress_result$predictor
}
}

DH_cluster_matrix <- DH_cluster[[1]]
for(j in 2:length(DH_cluster_num)){
DH_cluster_matrix <- cbind(DH_cluster_matrix,DH_cluster[[j]])
}

write.table(DH_cluster_matrix,file="DH_cluster_matrix.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Dis_matrix,file="distance_matrix.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
DH_coef_filenames <- rep(NA,length(DH_cluster_num))
DH_predictor_filenames <- rep(NA,length(DH_cluster_num))
for(j in 1:length(DH_cluster_num)){
write.table(coef_cluster[[j]][,-1],file=paste0("DH_cluster",DH_cluster_num[j],"_coef.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE) ##b0 is zero and excluded
write.table(predictor_idx_cluster[[j]],file=paste0("DH_cluster",DH_cluster_num[j],"_predictor.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
DH_coef_filenames[j] <- paste0("DH_cluster",DH_cluster_num[j],"_coef.txt")
DH_predictor_filenames[j] <- paste0("DH_cluster",DH_cluster_num[j],"_predictor.txt")
}

filenames <- paste(getwd(),c("gene_quantile.txt","gene_mean.txt","gene_sd.txt","cluster_idx.txt","DNase_mean.txt","DNase_sd.txt","regress_coef.txt","regress_predictor.txt","genomic_loci.txt","gene_name.txt","distance_matrix.txt","DH_cluster_matrix.txt",DH_coef_filenames,DH_predictor_filenames),sep="/")
param_list <- c(dim(DNase_data)[1],dim(Exon_train)[1],Num_cluster,Bin_size,Num_predictor,DH_cluster_num,filenames)
names(param_list) <- c("NumLoci","NumGene","NumCluster","NumBin","NumVar",paste0("DHCluterNum",c(1:length(DH_cluster_num))),"GeneQuantile","GeneMean","GeneSD","ClusterIndex","DNaseMean","DNaseSD","RegressionCoef","RegressionPredictor","GenomicLoci","GeneName","DistanceMatrix","DHCluster",paste0("DHClusterCoef",c(1:length(DH_cluster_num))),paste0("DHClusterPredictor",c(1:length(DH_cluster_num))))
write.table(param_list,file="par_file.txt",sep="\t",col.names=FALSE,quote=FALSE)

q(save="no")
