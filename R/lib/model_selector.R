# model selector
modelSelector <- function(results, obj=1){
  modop <- list(results[[2]][obj], c(results[[3]][[1]][obj], results[[3]][[2]][obj]))
  return(modop)
}