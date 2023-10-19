options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

#chr_ind <- as.numeric(args[1])
#cat(sprintf("%s", index), sep="\n")
DATA_PATH <- args[1]
RESULT_PATH <- args[2]
Rscript_path <- args[3]
Resolution <- as.numeric(args[4])


#DATA_PATH <- args[1]
#RESULT_PATH <- args[2]
#Rscript_path <- args[3]
#Resolution <- as.numeric(args[4])

library(msir, lib="/hpf/largeprojects/pmaass/Signature/pipeline/signature/scr/Rlib-4.2.1")
library(mclust, lib="/hpf/largeprojects/pmaass/Signature/pipeline/signature/scr/Rlib-4.2.1")
library(dplyr)
library(parallel)
library(MASS)

span_summary <- list()

##anchor_chromosomes = c(1:22, "X", "Y")
anchor_chromosomes = c(1:22, "X", "Y")

#anchor_chromosomes <- anchor_chromosomes[as.numeric(chr_ind) %% length(anchor_chromosomes) + 1]
#cat(sprintf("%s", anchor_chromosomes), sep="\n")

##for(k in anchor_chromosomes){
cis_cross <- function(k){

  
if (file.exists(sprintf("%s/chr%s.txt", RESULT_PATH, k)) == FALSE){
    
#target_chromosomes <- c(1:22, "X", "Y")
#for(i in target_chromosomes){
  
if (file.exists(sprintf("%s/chr%s.bed",DATA_PATH, k)) == TRUE){
interactions_table <- read.table(sprintf("%s/chr%s.bed",DATA_PATH, k))
#interactions_table <- interactions_table[1:14000,]
#interactions_table <-read.table(sprintf("%s",DATA_PATH))
#interactions_table<-read.table("chr1chr2.txt")
cat(sprintf("chr%s", k), sep="\n")

#source("LWLR_train_cis.R")
source(Rscript_path)

colnames(interactions_table)=c("chr1","st1","end1","chr2","st2","end2","dist","freq")
interactions_table <- interactions_table[order(interactions_table$dist),]
interactions_table$index <- 1:nrow(interactions_table)

######log2 transformation#################
interactions_table$freq <- scale(log2(interactions_table$freq))

#########smoothing parameter selection###########
data_length <- nrow(interactions_table) 

span_set <- c((1000000/Resolution)*5, (1000000/Resolution)*6, (1000000/Resolution)*7, (1000000/Resolution)*8)
#span_set <- c((1000000/Resolution)*5, (1000000/Resolution)*6, (1000000/Resolution)*7)

#span_set <- span_set <- c(5, 6, 7, 8, 9, 10, 11, 12)
fold_count <- 4


anchor_regions <- length(unique(interactions_table$st1))
if(anchor_regions <= min(span_set)){
  span_set = anchor_regions 
} else if(anchor_regions < max(span_set)){
  span_set <- span_set[1:(anchor_regions-(span_set[1]-1))]
} 

span_threshold = max(span_set)/anchor_regions 

###if (data_length < 200){
###  span_frac = 1
###} else if(span_threshold*data_length < 200){
###  span_frac = 200/data_length
###} else {

NRMSE_matrix <- matrix(NA, nrow = fold_count, ncol = length(span_set), dimnames = list(1:fold_count, sprintf("%g", span_set)))

set.seed(149)
#train_indices <- sample(interactions_table$index, ceiling(train_ratio * length(interactions_table$index)))
#allocation <- sample(rep(1:fold_count, ceiling(length(train_indices) / fold_count)), length(train_indices))
allocation <- sample(rep(1:fold_count, ceiling(length(interactions_table$index) / fold_count)), length(interactions_table$index))


for (fold in 1:fold_count) {
 # train_CV_indices <- c(train_indices[which(allocation != fold)], train_indices[which(allocation != fold)])
#  test_CV_indices <- c(train_indices[which(allocation == fold)], train_indices[which(allocation == fold)]) 
  train_CV_indices <- c(interactions_table$index[which(allocation != fold)])
  test_CV_indices <- c(interactions_table$index[which(allocation == fold)])
  
  
  X_train <- interactions_table[train_CV_indices,]
  X_test <- interactions_table[test_CV_indices,]
  distances <- X_train$dist
  
  for (span_size in span_set) {
    print(sprintf("running fold = %d, genomic regions in span = %g", fold, span_size))
    
    span_frac = span_size/length(unique(distances))
    sprintf("smoothing parameter is %s", span_frac)
    
    trained_set <- LWLR_train_cis(X_train , distances, span_frac)
    
    M <-aggregate(loess_y ~ dist, trained_set, mean )
    X_test <- merge(X_test, M, by="dist")
    
    NRMSE_matrix[fold, sprintf("%g", span_size)] <- mean(sqrt((X_test$freq - X_test$loess_y)^2) / sd(X_test$freq))
   # X_test <- select(X_test, -(loess_y))
     X_test <- X_test[,colnames(X_test)!="loess_y"]
  } 
  
}  
best_span <- span_set[max.col(t(colMeans(-NRMSE_matrix)), ties.method = "first")] 

cat(sprintf("optimal region number in span is %s", best_span), sep="\n")
dd <- c(sprintf("chr%s",k), best_span)
span_summary <- rbind(span_summary, dd)
 
span_frac = best_span/length(unique(interactions_table$dist))

#if (span_frac*data_length < 200) {
#  span_frac = 200/data_length
#}
###}


sprintf("smoothing parameter is %s", span_frac)

###########run the algorithm on full data#############

X_train <- interactions_table
distances <- interactions_table$dist
interactions_table <- LWLR_train_cis(X_train , distances, span_frac)

#####weighted mean and standard deviation for each start point of anchor chromosome####

M <-aggregate(loess_y ~ dist, interactions_table, mean )
S <-aggregate(loess.sd_y ~ dist, interactions_table, mean )

interactions_table <- merge(interactions_table, M, by="dist")
interactions_table <- merge(interactions_table, S, by="dist")

########filter data points with 0 SD########
interactions_table <- interactions_table[which(interactions_table$loess.sd_y.y != 0),]

#interactions_table <- select(interactions_table, -(loess_y.x))
#interactions_table <- select(interactions_table, -loess.sd_y.x)

#####calculate z-scores for each interaction#####
interactions_table$zscore<-(interactions_table$freq - interactions_table$loess_y.y)/interactions_table$loess.sd_y.y
interactions_table$zscore <- as.numeric(format(round(interactions_table$zscore, 3) , scientific=F))

###########converting zscore to pvalue using normal cumulative distribution##############
interactions_table$pvalue <- (1-pnorm(abs(-interactions_table$zscore)))*2
#interactions_table$pvalue <- as.numeric(format(round(interactions_table$pvalue, 7) , scientific=F))

interactions_table <- interactions_table[order(interactions_table$index),]
interactions_table$anchor_chromosome <- interactions_table$chr1

interactions_table <- interactions_table[,c("index", "chr1","st1","end1","chr2","st2","end2","dist","freq", "loess_y.y", "loess.sd_y.y", "zscore", "pvalue")]

write.table(interactions_table,file = sprintf("%s/chr%s.txt", RESULT_PATH, k), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}}
}
numCores <- detectCores()
sprintf("core number is %s", numCores)
mclapply(anchor_chromosomes, cis_cross, mc.cores = numCores)
##################save file##########################
#write.table(interactions_table,file = sprintf("%s/tmp.zscores.trans.LOESS1.txt", RESULT_PATH), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#write.table(interactions_table,gsub(" ", "", paste(DATA_PATH,"zscores.trans.LOESS1.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=F
#}
write.table(span_summary, file = sprintf("%s/span_summary_1MB.txt", RESULT_PATH), row.names = FALSE)
