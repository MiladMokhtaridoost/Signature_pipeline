library(dplyr)


#_____Set variables and create null dataframes__________________________________

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

pathway <- args[1]

unknown_list <- read.table(args[2])
unknown_list <- unknown_list[,1]

gonosome_unknown_data <- data.frame(matrix(ncol = 25, nrow = 1, "NA"))
colnames(gonosome_unknown_data) <- c('Name',
                                   'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chrX','chr8','chr9','chr11','chr10','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr20','chr19','chrY','chr22','chr21')




#_____Preparing input data______________________________________________________


unknown <- list(NULL)

for (i in 1:length(unknown_list)){
  sample=unknown_list[i]
  unknown[[i]]=list.files(pathway,pattern=sample)
}

# remove empty elements from list
unknown <- Filter(function(x) length(x) > 0, unknown)



# -- read -- #
print("reading in individual PMR files")
unknown_dat <- noquote(unlist(unknown))

read_function <- function(r){read.table(paste0(pathway,"/",unknown_dat[r]), header = T)}
unknown_dat_dfs <- lapply(1:length(unknown_dat), read_function)



# -- merge -- #
print("merge files into final dataframe")
for (i in 1:length(unknown_dat)){

  gonosome_unknown_data[i,1] <- gsub(".ProportionMappedReads.txt","",unknown_dat[i])
  gonosome_unknown_data[i,2:25] <- unknown_dat_dfs[[i]]$mapped_coverage

}

#remove any rows with empty runs
gonosome_unknown_data <- gonosome_unknown_data[apply(gonosome_unknown_data,1, function(x) all(x!=0)),]




#_____Modifying final datast____________________________________________________

# take average of autosomes
average <- lapply(1:nrow(gonosome_unknown_data), function(m){mean(as.numeric(gonosome_unknown_data[m,2:ncol(gonosome_unknown_data)]))})
gonosome_unknown_data$chrAll <- unlist(average)
gonosome_unknown_data <- gonosome_unknown_data %>% select(c("Name","chrAll","chrX","chrY"))

# re-formatting final 
gonosome_unknown_data <- data.frame(gonosome_unknown_data, row.names = 1)

# exporting
write.table(gonosome_unknown_data, file=paste0(pathway,"/gonosome_unknown_data.txt"), row.names = T, col.names = T, quote = F)
print(sprintf("Output: %s", paste0(pathway,"/gonosome_unknown_data.txt"))) 





