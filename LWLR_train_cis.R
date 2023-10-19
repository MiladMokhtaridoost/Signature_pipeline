LWLR_train_cis <- function(X_train, distances, span_frac) {
  
  ###############estimating interaction frequency using weighted regression################
  
  first.pass<-loess.sd(distances, X_train$freq, span=span_frac)
  X_train$loess_y<-first.pass$y  #fitted
  X_train$loess.sd_y<-first.pass$sd
  
  return(X_train)
}
