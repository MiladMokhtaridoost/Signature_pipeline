LWLR_train <- function(X_train, anchor, span_frac) {
library(purrr)  

###############estimating interaction frequency using weighted regression################

possloes <- possibly(.f = loess.sd, otherwise = "unfit_span")
first.pass <- possloes(anchor, X_train$freq, span=span_frac)

if("unfit_span" %in% first.pass){
  X_train$loess_y <- 0  #fitted
  X_train$loess.sd_y <- 0.001
} else {
X_train$loess_y<-first.pass$y  #fitted
X_train$loess.sd_y<-first.pass$sd
}

return(X_train)
}
