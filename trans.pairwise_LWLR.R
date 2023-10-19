options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- args[1]
RESULT_PATH <- args[2]
Rscript_PATH <- args[3]
Resolution <- as.numeric(args[4])
print(DATA_PATH)
print(RESULT_PATH)



library(msir, lib="/hpf/largeprojects/pmaass/Signature/pipeline/signature/scr/Rlib-4.2.1")
library(mclust, lib="/hpf/largeprojects/pmaass/Signature/pipeline/signature/scr/Rlib-4.2.1")
library(dplyr)
library(parallel) 
library(MASS) 

source(Rscript_PATH)
#interactions_table <-read.table(sprintf("%s",DATA_PATH))

#anchor_chromosomes = c("Y")
#for(k in anchor_chromosomes){
#target_chromosomes <- c(1:22, "X", "Y")
#for(i in target_chromosomes){

#if (file.exists(sprintf("%s/chr%schr%s.txt",DATA_PATH, k, i)) == TRUE){
#interactions_table <- read.table(sprintf("%s/chr%schr%s.txt",DATA_PATH, k, i))
#trans_cross <- function(k){ 

###
#interactions_table <- interactions_table_anchor
###
trans_cross <- function(input_data){

interactions_table <- read.table(sprintf("%s",input_data))
  
  
#source(Rscript_PATH)

colnames(interactions_table)=c("id","chr1","st1","end1","chr2","st2","end2","dist","freq")
interactions_table$index <- 1:nrow(interactions_table)

######log2 transformation#################
interactions_table$freq <- scale(log2(interactions_table$freq))

#########smoothing parameter selection###########
data_length <- nrow(interactions_table)

span_set <- c((1000000/Resolution)*5, (1000000/Resolution)*6, (1000000/Resolution)*7, (1000000/Resolution)*8, (1000000/Resolution)*9, (1000000/Resolution)*10)
#span_set <- span_set <- c(5, 6, 7, 8, 9, 10, 11, 12)
fold_count <- 4


anchor_regions <- length(unique(interactions_table$st1))
if(anchor_regions <= min(span_set)){
  span_set = anchor_regions
} else if(anchor_regions < max(span_set)){
  span_set <- span_set[1:(anchor_regions-(span_set[1]-1))]
}

span_threshold = max(span_set)/anchor_regions

if (data_length < 200){
  span_frac = 1
} else if(span_threshold*data_length < 200){
  span_frac = 200/data_length
} else {

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
  anchor <- X_train$st1

  for (span_size in span_set) {
    print(sprintf("running fold = %d, genomic regions in span = %g", fold, span_size))

    span_frac = span_size/length(unique(anchor))
    sprintf("smoothing parameter is %s", span_frac)

    trained_set <- LWLR_train(X_train , anchor, span_frac)

    M <-aggregate(loess_y ~ st1, trained_set, mean )
    X_test <- merge(X_test, M, by="st1")

    NRMSE_matrix[fold, sprintf("%g", span_size)] <- mean(sqrt((X_test$freq - X_test$loess_y)^2) / sd(X_test$freq))
 #   X_test <- select(X_test, -(loess_y))  
     X_test <- X_test[,colnames(X_test)!="loess_y"]

  }

}
best_span <- span_set[max.col(t(colMeans(-NRMSE_matrix)), ties.method = "first")]

cat(sprintf("optimal region number in span is %s", best_span), sep="\n")
#dd <- c(sprintf("chr%s(anchor)chr%s",k,i), best_span)
#span_summary <- rbind(span_summary, dd)

span_frac = best_span/length(unique(interactions_table$st1))

#if (span_frac*data_length < 200) {
#  span_frac = 200/data_length
#}
}


sprintf("smoothing parameter is %s", span_frac)

###########run the algorithm on full data#############

X_train <- interactions_table
anchor <- interactions_table$st1
interactions_table <- LWLR_train(X_train , anchor, span_frac)

#####weighted mean and standard deviation for each start point of anchor chromosome####

M <-aggregate(loess_y ~ st1, interactions_table, mean )
S <-aggregate(loess.sd_y ~ st1, interactions_table, mean )

interactions_table <- merge(interactions_table, M, by="st1")
interactions_table <- merge(interactions_table, S, by="st1")

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

interactions_table <- interactions_table[,c("index", "id","chr1","st1","end1","chr2","st2","end2","dist","freq", "loess_y.y", "loess.sd_y.y", "zscore", "pvalue", "anchor_chromosome")]

anc <- interactions_table$chr1[1]
tar <- interactions_table$chr2[1]

##################save file##########################
#write.table(interactions_table,gsub(" ", "", paste(DATA_PATH,"zscores.trans.LOESS1.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
write.table(interactions_table,file = sprintf("%s/%s%stmp.zscores.trans.LOESS1.txt", RESULT_PATH, anc, tar), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#write.table(interactions_table,gsub(" ", "", paste(DATA_PATH,"zscores.trans.LOESS1.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
######
} ###### end of function
numCores <- detectCores()
sprintf("core number is %s", numCores)
mclapply(DATA_PATH, trans_cross, mc.cores = numCores)
#lapply(DATA_PATH, trans_cross)


#################Second chromosome as anchor####################
#interactions_table2 <- read.table(sprintf("%s",DATA_PATH))

trans_cross2 <- function(input_data2){

interactions_table2 <- read.table(sprintf("%s",input_data2))
  
  
colnames(interactions_table2)=c("id","chr1","st1","end1","chr2","st2","end2","dist","freq")
interactions_table2$index <- 1:nrow(interactions_table2)
######log2 transformation#################
interactions_table2$freq <- scale(log2(interactions_table2$freq))

##########################################################

data_length <- nrow(interactions_table2)

span_set <- c((1000000/Resolution)*5, (1000000/Resolution)*6, (1000000/Resolution)*7, (1000000/Resolution)*8, (1000000/Resolution)*9, (1000000/Resolution)*
10)

#span_set <- c(5, 6, 7, 8, 9, 10, 11, 12)
fold_count <- 4
#train_ratio <- 0.8

anchor_regions <- length(unique(interactions_table2$st2))

if(anchor_regions <= min(span_set)){
  span_set = anchor_regions
} else if(anchor_regions < max(span_set)){
  span_set <- span_set[1:(anchor_regions-(span_set[1]-1))]
}

span_threshold = max(span_set)/anchor_regions

if (data_length < 200){
  span_frac = 1
} else if(span_threshold*data_length < 200){
  span_frac = 200/data_length
} else {

NRMSE_matrix2 <- matrix(NA, nrow = fold_count, ncol = length(span_set), dimnames = list(1:fold_count, sprintf("%g", span_set)))

set.seed(623)
#train_indices <- sample(interactions_table$index, ceiling(train_ratio * length(interactions_table$index)))
#allocation <- sample(rep(1:fold_count, ceiling(length(train_indices) / fold_count)), length(train_indices))
allocation <- sample(rep(1:fold_count, ceiling(length(interactions_table2$index) / fold_count)), length(interactions_table2$index))


for (fold in 1:fold_count) {
  # train_CV_indices <- c(train_indices[which(allocation != fold)], train_indices[which(allocation != fold)])
  #  test_CV_indices <- c(train_indices[which(allocation == fold)], train_indices[which(allocation == fold)])
  train_CV_indices <- c(interactions_table2$index[which(allocation != fold)])
  test_CV_indices <- c(interactions_table2$index[which(allocation == fold)])


  X_train <- interactions_table2[train_CV_indices,]
  X_test <- interactions_table2[test_CV_indices,]

    anchor <- X_train$st2

    for (span_size in span_set) {
    print(sprintf("running fold = %d, genomic regions in span = %g", fold, span_size))

    span_frac = span_size/length(unique(anchor))
    sprintf("smoothing parameter is %s", span_frac)

    trained_set <- LWLR_train(X_train , anchor, span_frac)

    M <-aggregate(loess_y ~ st2, trained_set, mean )
    X_test <- merge(X_test, M, by="st2")

    NRMSE_matrix2[fold, sprintf("%g", span_size)] <- mean(sqrt((X_test$freq - X_test$loess_y)^2) / sd(X_test$freq))
#    X_test <- select(X_test, -(loess_y))
    X_test <- X_test[,colnames(X_test)!="loess_y"]

	}

}
best_span <- span_set[max.col(t(colMeans(-NRMSE_matrix2)), ties.method = "first")]

cat(sprintf("optimal region number in span is %s", best_span), sep="\n")
#ff <- c(sprintf("chr%schr%s(anchor)",k,i), best_span)
#span_summary <- rbind(span_summary, ff)

span_frac = best_span/length(unique(interactions_table2$st2))

#if (span_frac*data_length < 200) {
#  span_frac = 200/data_length
#}
}


sprintf("smoothing parameter is %s", span_frac)

###########run the algorithm on full data#############

X_train <- interactions_table2
anchor <- interactions_table2$st2
interactions_table2 <- LWLR_train(X_train , anchor, span_frac)

#####weighted mean and standard deviation for each start point of anchor chromosome####

M <-aggregate(loess_y ~ st2, interactions_table2, mean )
S <-aggregate(loess.sd_y ~ st2, interactions_table2, mean )

interactions_table2 <- merge(interactions_table2, M, by="st2")
interactions_table2 <- merge(interactions_table2, S, by="st2")

########filter data points with 0 SD########
interactions_table2 <- interactions_table2[which(interactions_table2$loess.sd_y.y != 0),]

#interactions_table2 <- select(interactions_table2, -(loess_y.x))
#interactions_table2 <- select(interactions_table2, -loess.sd_y.x)

#####calculate z-scores for each interaction#####
interactions_table2$zscore<-(interactions_table2$freq - interactions_table2$loess_y.y)/interactions_table2$loess.sd_y.y
interactions_table2$zscore <- as.numeric(format(round(interactions_table2$zscore, 3) , scientific=F))

###########converting zscore to pvalue using normal cumulative distribution##############
interactions_table2$pvalue <- (1-pnorm(abs(-interactions_table2$zscore)))*2
#interactions_table2$pvalue <- as.numeric(format(round(interactions_table2$pvalue, 7) , scientific=F))

interactions_table2 <- interactions_table2[order(interactions_table2$index),]
interactions_table2$anchor_chromosome <- interactions_table2$chr2

interactions_table2 <- interactions_table2[,c("index", "id","chr1","st1","end1","chr2","st2","end2","dist","freq", "loess_y.y", "loess.sd_y.y", "zscore", "pvalue", "anchor_chromosome")]

anc <- interactions_table2$chr1[1]
tar <- interactions_table2$chr2[1]

#write.table(interactions_table2,gsub(" ", "", paste(DATA_PATH,"zscores.trans.LOESS2.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
write.table(interactions_table2, file = sprintf("%s/%s%stmp.zscores.trans.LOESS2.txt", RESULT_PATH,anc,tar), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#write.table(interactions_table2,gsub(" ", "", paste(DATA_PATH,"zscores.trans.LOESS2.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#}}}
#######
} ###### end of function
numCores2 <- detectCores()
sprintf("core number is %s", numCores2)
mclapply(DATA_PATH, trans_cross2, mc.cores = numCores2)
#lapply(DATA_PATH, trans_cross2)



#warnings()
