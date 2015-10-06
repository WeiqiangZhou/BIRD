########################################################################################################################
#                                                       README
#This file contains the R code for selecting the optimal parameters for BIRD using 5-fold cross-validation on 1% of DHSs.
#The number of gene clusters K is set to 100,200,500,1000,1500 and 2000.
#The number of predictors N is set to 1,2,3,4,5,6,7 and 8.
########################################################################################################################

####functions####
##regression##
cluster_regression_topn <- function(DNase_train,Exon_train_mean,Exon_test_mean,top_n){
y_cor <- cor(DNase_train,t(Exon_train_mean))
y_cor[is.na(y_cor)] <- 0

if(top_n > 1){
max_idx <- sort(y_cor, decreasing=TRUE, index.return=TRUE)$ix[1:top_n]      #if the N>1, select the N most correlated predictiors
data_train <- data.frame(y=DNase_train, t(Exon_train_mean[max_idx,]))
data_test <- data.frame(t(Exon_test_mean[max_idx,]))
fit <- lm(y~.,data_train)
y_pre <- predict(fit,data_test)
}else{
max_idx <- which(y_cor == max(y_cor))[1]      #if N=1, select the most correlated predictor
data_train <- data.frame(y=DNase_train, x=Exon_train_mean[max_idx,])
data_test <- data.frame(x=Exon_test_mean[max_idx,])
fit <- lm(y~x, data_train)
y_pre <- predict(fit,data_test)
}

return(y_pre)
}

##standardize data
standardize_rows <- function(Data_train, Data_test){

  train_mean_sd <- matrix(data=NA, nrow=dim(Data_train)[1], ncol=2)
  Data_train_sd <- matrix(data=NA, nrow=dim(Data_train)[1], ncol=dim(Data_train)[2])
  Data_test_sd <- matrix(data=NA, nrow=dim(Data_test)[1], ncol=dim(Data_test)[2])
  for(i in 1:dim(Data_train)[1]){
   train_mean_sd[i,1] <- mean(Data_train[i,])
   train_mean_sd[i,2] <- sd(Data_train[i,])
   if(train_mean_sd[i,2]==0){
     Data_train_sd[i,] <- rep(0,dim(Data_train)[2])
     Data_test_sd[i,] <- rep(0,dim(Data_test)[2])
   }
   else{
   Data_train_sd[i,] <- (Data_train[i,] - train_mean_sd[i,1])/train_mean_sd[i,2]
   Data_test_sd[i,] <- (Data_test[i,] - train_mean_sd[i,1])/train_mean_sd[i,2]
   }
  }

  result <- list("mean_sd"=train_mean_sd, "train"=Data_train_sd, "test"=Data_test_sd)

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

##main script##
Exon_data <- as.matrix(read.table(file="Exon_data_train.txt",header=TRUE,row.names=1))   #input gene expression data
load("DNase_data_train.rda")  #input DNase-seq data, first three columns contain the genomic locus information 

set.seed(2015)
DNase_sub <- as.matrix(DNase_data[,-c(1:3)])[sample(1:dim(DNase_data)[1],floor(0.01*dim(DNase_data)[1])),]  #randomly sample 1% of DHSs
group_idx <- sample(cut(seq(1,dim(Exon_data)[2]),breaks=5,label=FALSE))  #randomly partition the data into 5 groups

Num_cluster <- c(100,200,500,1000,1500,2000) #set range of K
Num_predictor <- c(1:8)  #set range of N

##5-fold cross-validation##
K_result <- list()
for(K in 1:length(Num_cluster)){
group_result <- list()
for(group in 1:5){

Exon_test <- Exon_data[,which(group_idx==group)]
Exon_train <- Exon_data[,setdiff(1:dim(Exon_data)[2],which(group_idx==group))]

Exon_processed <- standardize_rows(Exon_train,Exon_test) #standardize gene expression data
Exon_train_sd <- Exon_processed$train
Exon_test_sd <- Exon_processed$test

Exon_cluster <- kmeans(Exon_train_sd, centers=Num_cluster[K], nstart=10, iter.max=50)  #cluster gene expression
 
Exon_train_mean <- cluster_data_mean(Exon_train_sd,Exon_cluster$cluster)
Exon_test_mean <- cluster_data_mean(Exon_test_sd,Exon_cluster$cluster)

DNase_test <- DNase_sub[,which(group_idx==group)]
DNase_train <- DNase_sub[,setdiff(1:dim(DNase_sub)[2],which(group_idx==group))]

DNase_processed <- standardize_rows(DNase_train,DNase_test) #standardize DNase-seq data
DNase_mean_sd <- DNase_processed$mean_sd
DNase_train_sd <- DNase_processed$train

N_result <- list()
for(N in 1:length(Num_predictor)){

Y_pre <- matrix(data=NA,nrow=dim(DNase_train_sd)[1],ncol=dim(Exon_test_sd)[2])
for (i in 1:dim(DNase_train_sd)[1]){

if(DNase_mean_sd[i,2]!=0){
Y_pre[i,] <- cluster_regression_topn(DNase_train_sd[i,],Exon_train_mean,Exon_test_mean,Num_predictor[N]) #make prediction for each DHS
}
else{
Y_pre[i,] <- 0
}

Y_pre[i,] <- Y_pre[i,]*DNase_mean_sd[i,2] + DNase_mean_sd[i,1]
}
N_result[[N]] <- Y_pre
}

group_result[[group]] <- N_result

}
K_result[[K]] <- group_result
}

##evaluation##
##find the optimal K and N based on the average cross-cell-type correlation
DNase_true <- DNase_sub[,which(group_idx==1)]
for(i in 2:5){
DNase_true <- cbind(DNase_true,DNase_sub[,which(group_idx==i)])
}

evaluate_result <- matrix(data=NA,ncol=length(Num_cluster),nrow=length(Num_predictor))

for(K in 1:length(Num_cluster)){
for(N in 1:length(Num_predictor)){

DNase_pre <- K_result[[K]][[1]][[N]]
for(i in 2:5){
DNase_pre <- cbind(DNase_pre,K_result[[K]][[i]][[N]])
}

locus_cor <- rep(0,dim(DNase_true)[1])
for(i in 1:dim(DNase_true)[1]){
locus_cor[i] <- cor(DNase_true[i,],DNase_pre[i,])
}

locus_cor[is.na(locus_cor)] <- 0
evaluate_result[N,K] <- mean(locus_cor)
}
}

N_max <- Num_predictor[which(evaluate_result == max(evaluate_result), arr.ind = TRUE)[1,1]]
K_max <- Num_cluster[which(evaluate_result == max(evaluate_result), arr.ind = TRUE)[1,2]]

write.table(evaluate_result,file="cross_validation_result_final.txt",sep="\t") #record the prediction performance for all K,N pairs
write.table(c(N_max,K_max),file="best_parameter_final.txt") #record the optimal K,N
runtime <- proc.time()
save(runtime,file="parameter_select_time.rda") #record the computational time
q(save="no")
