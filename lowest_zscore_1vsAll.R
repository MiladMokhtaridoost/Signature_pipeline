options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
tmp <- args[1]


#read in data
interactions<-read.table(tmp)
#1vsAll
colnames(interactions)=c("id","chrA","stA","endA","chrB","stB","endB","freq", "mean", "sd", "zscore", "pvalue", "anchor_chromosome", "qvalue")
                          
print("order all zscores")
interactions$interID <- do.call(paste0, interactions[c("chrA", "stA", "endA", "chrB", "stB", "endB")])
sDat <- interactions[order(interactions$interID, abs(interactions$zscore) ), ]

print("remove absolute highest zscore")
uDat <- sDat[ !duplicated(sDat$interID), ]

##################save file##########################
print("write output file")
write.table(uDat,gsub(" ", "", paste(tmp,".tmp.zscores.trans.lowest.score.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
