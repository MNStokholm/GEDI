#' @importFrom stats binomial complete.cases quantile
#' @importFrom utils read.table
#' @importFrom Rdpack reprompt

.checkClass<-function(X,X.class){
  if (class(X.class) != "character") stop(gettextf("Object %s must be of class %s",
                                                   sQuote("X.class"), dQuote("character")))
  X.name = deparse(substitute(X))
  if(!all(class(X) %in% X.class)){
    X.class.msg = paste(paste(dQuote(X.class[-length(X.class)]),collapse=", "),
                        "or",dQuote(X.class[length(X.class)]),sep=" ")
    stop(gettextf("Object %s must be of class %s",
                  sQuote(X.name), X.class.msg))
  }else if(X.class[1] == "logical"){
    if(is.na(X)[1]) stop(gettextf("%s is NA when TRUE/FALSE is required",sQuote(X.name)))
  }
}

.onlyFirstElement <- function(X){
  X.name = deparse(substitute(X))
  if(length(X) > 1){
    warning(gettextf("%s has length > 1. Only the first element, %s, is used",
                     sQuote(X.name),dQuote(X[1])))
  }
  return(X[1])
}


.NAtoFALSE <- function(bool){
  if (is.na(bool)) return(FALSE)
  else if (bool %in% c(TRUE,FALSE)) return(bool)
  else return(NA)
}


