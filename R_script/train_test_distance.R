##This function is used to calculate the distance between the input samples and the training samples in the model.
library(preprocessCore)
library(data.table)

get_dist <- function(train_expr,input_expr,cluster_num=1000){

	input_expr <- as.matrix(input_expr)
	match_idx <- match(row.names(input_expr),row.names(train_expr))
	input_expr_match <- matrix(data=0, ncol=ncol(input_expr), nrow=nrow(train_expr))
	input_expr_match[match_idx[which(match_idx != "NA")],] <- input_expr[which(match_idx != "NA"),]
	colnames(input_expr_match) <- colnames(input_expr)
	row.names(input_expr_match) <- row.names(train_expr)

	row_mean <- rowMeans(train_expr)
	row_sd <- apply(train_expr,1,sd)
	train_expr_center <- (train_expr - row_mean)/row_sd

	cluster_id <- kmeans(train_expr_center, centers=cluster_num, nstart=10, iter.max=50)$cluster

	train_quantile <- normalize.quantiles.determine.target(train_expr)
	input_expr_norm <- (normalize.quantiles.use.target(input_expr_match,train_quantile) - row_mean)/row_sd

	train_expr_cluster <- t(sapply(c(1:max(cluster_id)),function(x){
	  colMeans(as.matrix(train_expr_center[which(cluster_id==x),]))
	}))

	input_expr_cluster <- t(sapply(c(1:max(cluster_id)),function(x){
	  colMeans(as.matrix(input_expr_norm[which(cluster_id==x),]))
	}))

	output <- sapply(c(1:ncol(input_expr_cluster)),function(x){
    cor_tmp <- cor(input_expr_cluster[,x],train_expr_cluster)
    1 - mean(sort(cor_tmp,decreasing=TRUE)[1:round(ncol(train_expr_cluster)*0.1)])
  })

	return(output)
}

##example
##read in training RNA-seq data
train_expr <- readRDS("/BIRD_dir/data/RNA_data_167_cells.rds")
##read in input data, first column should be ensembl gene ids
input_expr <- data.frame(fread("input_expr.txt"),row.names=1)
##get training-test distance
train_test_dist <- get_dist(train_expr,input_expr)
