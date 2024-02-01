################################################################################
# Machine learning approach for gonosome determination: logistic regression
#
# Purpose: Gonosomal characterization is important for sex-specific analyses. 
#          This demographic is typically not recorded in papers which is why 
#          a machine learning approach was utilized to classify a sample/dataset
#          as either male (XY) or female (XX). 
#
# Input: The merged and processed Proportion Mapped Reads file (one male one 
#        female)
#
# Developed by: Jordan Chalmers
#               email: jordanchalmers9@gmail.com
#
################################################################################
library(dplyr)


#_____Load in data______________________________________________________________

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

output <- args[1]

data_f <- read.table(args[2])
data_m <- read.table(args[3])

unknown_data <- read.table(args[4])
unknown_list <- read.table(args[5])
unknown_list <- unknown_list[,1]


#_____Train the model on the full dataset_______________________________________

# full dataset
full_train_data <- rbind(data_f,data_m)

# set each binary label as either 1 or 0 for later specification in predication
full_train_data$Label[which(full_train_data$Label == "female")] <- 1
full_train_data$Label[which(full_train_data$Label == "male")] <- 0

# fit the model
full_train_data$Label <- as.factor(full_train_data$Label)
GD_model <- glm(Label~., data = full_train_data, family = binomial)



#_____Run the model on your unknown data________________________________________

# run predication using trained model
probabilities <- GD_model %>% predict(unknown_data, type = "response")
predicted <- ifelse(probabilities > 0.5, "female", "male")

# export predicted values
predicted <- as.data.frame(predicted)
predicted <- data.frame("sample"=rownames(predicted), predicted)
write.csv(predicted, paste0(output,"/gonosome_prediction_raw.csv"), row.names = F)



#_____Ensure all runs are concordant____________________________________________
library(stringr)

predicted_final <- data.frame(matrix(ncol = 2, nrow = length(unknown_list), "NA"))
colnames(predicted_final) <- c('sample','predicted.sex')


for (i in 1:length(unknown_list)){

  # subset based on sample name
  sample=unknown_list[i]
  predicted_final$sample[i] <- sample
  index <- grep(sample, predicted$sample)
  
  # calculate percent female
  predicted[index,]
  ct.f <- sum(str_count(predicted[index,]$predicted, "female"))
  total <- length(predicted[index,]$predicted)
  
  # annotate
  if (ct.f / total > 0.5){
    predicted_final$predicted.sex[i] <- "female"
  } else if (ct.f / total < 0.5) {
    predicted_final$predicted.sex[i] <- "male"
  } else {
    predicted_final$predicted.sex[i] <- "undetermind"
  }

}


# export final predicted values
write.csv(predicted_final, paste0(output,"/gonosome_prediction_final.csv"), row.names = F)
paste0("Final gonosome prediction output: ",output,"/gonosome_prediction_final.csv")
