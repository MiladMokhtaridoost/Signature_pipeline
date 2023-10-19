options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
tmp <- args[1]
#install.packages("dplyr", repos="http://cran.r-project.org", lib="/hpf/largeprojects/pmaass/3D-flow/scripts/src/R-lib")
#install.packages("FSA", repos="http://cran.r-project.org", lib="/hpf/largeprojects/pmaass/3D-flow/scripts/src/R-lib")


#########################################################################

#read in data
#tmp <- "pairwise_test_dat.txt"
#tmp <- "two_zscore_values_1vsAll_comparison/Cardiomyocites_primitive_rep1.1000000_zscore_1vAll.txt"
interactions<-read.table(tmp)
##pairwise
colnames(interactions)=c("index","id","chr1","st1","end1","chr2","st2","end2","dist","freq","mean","sd","zscore","pvalue")
#interactions <- interactions[,-1]
interactions$interID <- do.call(paste0, interactions[c("chr1","st1","end1", "chr2", "st2", "end2")])
#summary(interactions)
#head(interactions)
#interactions[which(interactions$interID == "chr1001000000chr1201000000"),]
sDat <- interactions[order(interactions$interID, abs(interactions$zscore) ), ]
#sDat[which(sDat$interID == "chr1001000000chr1201000000"),]
#remove absolute value highest zscore
uDat <- sDat[ !duplicated(sDat$interID), ]
#uDat[which(uDat$interID == "chr1001000000chr1201000000"),]

##################save file##########################
write.table(uDat,gsub(" ", "", paste(tmp,".tmp.zscores.trans.lowest.score.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
