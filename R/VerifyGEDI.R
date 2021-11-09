#' @title Verify transcriptomic data integration
#'
#' @description
#' \code{VerifyGEDI} allows users to verify if the \code{\link{GEDI}} function
#' successfully integrated the transcriptomic datasets. Before using this function,
#' the \code{\link{BatchCorrection}} function must remove the integrated dataset's
#' batch effect. The transcriptomic data integration is verified using a supervised
#' machine learning model to predict the samples' status in one batch based on the
#' remaining batches' samples. The forward stepwise selection (FSS) algorithm picks
#' which genes to use as features in the chosen model.
#' @name VerifyGEDI
#' @usage
#' VerifyGEDI(X, y, batch, model = "logistic", stop.at.score = 1,
#'            verbose = TRUE)
#' @param X A \code{data.frame} or \code{matrix} of integrated gene expression
#' data with batch effect removed. Columns contain samples, and rows are genes.
#'
#' @param y A character vector containing a status variable for the samples.
#' For example, healthy/sick or treated/control. \code{y} does not have to
#' be binary. \code{y} must satisfy: \code{length(y) = ncol(X)}.
#'
#' @param batch Character vector that describes which sample is in what batch.
#' \code{batch} must satisfy: \code{length(batch)=ncol(Data)}.
#'
#' @param model (Optional) A character describing which model \code{VerifyGEDI}
#' uses to predict the status \code{y} in each sample in every batch.
#' This function currently supports two models; a logistic regression model
#' (\dQuote{logistic}) and a decision tree classifier (\dQuote{tree}).
#' \dQuote{logistic} only works if Y contains two classes.
#'
#' @param stop.at.score (Optional) The threshold is a numeric value between 0 and 1.
#' This value determines at what score to force the FSS algorithm to terminate.
#' The score is how many per cent of the samples are correctly classified.
#' The default value is 1, meaning that the FSS algorithm will end when the model
#' makes 0 errors.
#'
#' @param verbose (Optional) Logical value. If \code{TRUE}, the function prints
#' descriptive messages about the progress.
#'
#'
#' @details
#' There are more genes than samples in gene expression datasets in almost all
#' cases. Therefore, the \code{VerifyGEDI} function must select some of the genes
#' as features for the classifier. \code{VerifyGEDI} uses the Forward stepwise
#' selection (FSS) algorithm to do this.
#'
#' In each iteration of the FSS algorithm, the chosen model is cross-validated
#' across all the batches. The model uses one batch as test data and the
#' remaining batches as training data. This is done for all the batches. The
#' algorithm uses the lowest test score across all the batches to evaluate which
#' gene is the best new feature for the model. The FSS algorithm will continue
#' to do this and add features to the model until the model no longer improves
#' or that the test score exceeds the threshold determined by the
#' \code{stop.at.score} argument.
#'
#' @return
#' A list with the following elements:
#' \describe{
#' \item{\code{$Feature.idx}}{Indices of rows in \code{X} used as features in the model.}
#' \item{\code{$Feature.names}}{Names of rows in \code{X} used as features in the model.
#' These would be Ensembl gene IDs if \code{X} was created using the \code{\link{GEDI}}
#' function.}
#' \item{\code{$TrainScores}}{Numeric vector. How the model scores on the training data for each batch.}
#' \item{\code{$TestScores}}{Numeric vector. Describes how the model scores on each batch
#' when training on the remaining batches.}
#' }
#'
#' @importFrom progress progress_bar
#' @importFrom stats glm predict
#' @importFrom rpart rpart
#'
#' @seealso
#' \code{\link{ReadGE}}
#' \code{\link{GEDI}}
#' \code{\link{BatchCorrection}}
#' @examples
#' \dontrun{
#' ### From README.md example
#' ## Continuing example from BatchCorrection()'s example,
#' ## see help("BatchCorrection")
#'
#' # The gene expression data integration is verified:
#' res <- VerifyGEDI(X = cData, y = status, batch = batch, model = "logistic")
#' #> Number of features:  2
#' #> ENSBTAG00000003530, ENSBTAG00000000476
#' #> Worst score across batches:  1
#'
#' # With an accuracy of 100%, it can be concluded
#' # that the data integration was successful
#' }
#' @export

VerifyGEDI <- function(X,y,batch,model = "logistic",stop.at.score = 1,verbose = TRUE){
  ## Input control ##
  if(class(model)=="character") model <- .onlyFirstElement(model)
  .VerifyGEDI.checkInput(X,y,batch,model,stop.at.score,verbose)
  stop.at.score <- .onlyFirstElement(stop.at.score)

  X <- data.matrix(t(X))
  y <- as.factor(y)

  var.n <- ncol(X)
  var.av <- c(1:var.n)
  var.ch.best <- c()
  score.best <- 0
  y0 = levels(y)[which.max(summary(y))]
  score <- sum((y==y0))/length(y);score
  scores<-NULL

  while(score > score.best & score < stop.at.score){
    score.best <- score
    if(!is.null(scores)){
      var.ch.best <- c(var.ch.best,var.av[which.max(scores)])
      var.av <- var.av[-which.max(scores)]
    }

    pb <- progress_bar$new(total = length(var.av),
                           format = paste("FSS algoithm step",length(var.ch.best)+1,"[:bar] :percent eta: :eta"))
    scores = c()
    for(attr in var.av){
      pb$tick()
      var.ch <- c(var.ch.best,attr)
      scores <- c(scores, suppressWarnings(min(.batchCV(X[,var.ch],y,batch,model = model)$TestScores)))
    }
    score = max(scores)
  }

  if (score > score.best){
    score.best <- score
    var.ch.best <- c(var.ch.best,var.av[which.max(scores)])
  }
  if (verbose) {
    #cat("__________________________________________________\n")
    #cat("Final result\n")
    cat("Number of features: ",length(var.ch.best),"\n")
    if(length(var.ch.best)>0) {cat(paste(colnames(X)[var.ch.best],collapse=", "),"\n")}
    cat("Worst score across batches: ",score.best,"\n")
  }

  ## This plot is quite useless
  # if(length(var.ch.best)>=2){
  #   plot.data = data.frame(X[,var.ch.best[1:2]],y,batch)
  #   colnames(plot.data)[3:4] = c("Status","Batch")
  #   p<-ggplot(data=plot.data, aes_string(x=colnames(plot.data)[1],
  #                                        y=colnames(plot.data)[2],
  #                                        color="Status", shape="Batch")) +
  #     geom_point(size=3) +
  #     xlab(colnames(plot.data)[1]) +
  #     ylab(colnames(plot.data)[2]) +
  #     ggtitle("Scatterplot with the first two selected genes")
  #
  #   print(p)
  # }

  finalScores <- .batchCV(X[,var.ch.best],y,batch,model = model)
  out <- list(Feature.idx = var.ch.best, Feature.names = colnames(X)[var.ch.best],
              TrainScores = t(finalScores)[1,],
              TestScores= t(finalScores)[2,])

  return(out)
}

.VerifyGEDI.checkInput <- function(X,y,batch,model,stop.at.score,verbose){
  .checkClass(X,c("data.frame","matrix","array"))
  if(any(dim(X)==0)){
    stop(gettextf("'X' is empty"))
  }else if(any(apply(X,2,class)!="numeric")){
    colError <- match(TRUE,apply(X,2,class)!="numeric")
    stop(gettextf("All columns in 'X' must all be of class %s or %s. Column %d has class %s",
                  dQuote("numeric"),dQuote("integer"),colError,dQuote(class(X[,colError]))))
  }
  .checkClass(y,c("character","numeric","factor"))
  if(ncol(X) != length(y)){
    stop(gettextf("'ncol(X) == length(y)' is not TRUE.
                  'ncol(X)=%d' and 'length(y)=%d'",
                  ncol(X),length(y)))
  }
  .checkClass(batch,c("character","factor","numeric"))
  if(length(y) != length(batch)){
    stop(gettextf("'length(y) == length(batch)' is not TRUE.
                  'length(y)=%d' and 'length(batch)=%d'",
                  length(y),length(batch)))
  }
  .checkClass(model,"character")
  validModels=c("logistic","tree")
  if(!model %in% validModels){
    n.models = length(validModels)
    model.msg = paste(paste(dQuote(validModels[-n.models]),collapse=", "),
                      "and",dQuote(validModels[n.models]),sep=" ")
    stop(gettextf("Invalid 'model' chosen. Valid options are: %s",model.msg))
  }
  if(model=="logistic" & length(levels(as.factor(y))) != 2){
    stop(gettextf("Logistic regression model requires %s to have 2 levels, not %d",
                  sQuote("y"),length(levels(as.factor(y)))))
  }

  .checkClass(stop.at.score,"numeric")
  if (stop.at.score <= 0 | 1 < stop.at.score){
    stop("'stop.at.score' out of range. Must be a value between 0 and 1.")
  }
  .checkClass(verbose,"logical")
}

.batchCV <- function(X,y,batch,model){
  err.train <- c()
  err.test <- c()
  batch = as.factor(batch)
  n.batches <- length(levels(as.factor(batch)))
  for(i in 1:n.batches){
    if(class(X)[1]=="numeric") X=matrix(X,nrow = length(y))
    test.batch <- batch==levels(batch)[i]
    #X.test=data.frame(t(X[i,]));y.test=data.frame(t(y[test.batch]))

    X.test=X[test.batch,];y.test=y[test.batch]
    X.train=X[!test.batch,];y.train=y[!test.batch]
    trainData<-data.frame(y.train,X.train);colnames(trainData)[1]="y"
    testData<-data.frame(y.test,X.test);colnames(testData)=colnames(trainData)

    if(model == "logistic"){
      trainData$y = unclass(trainData$y)-1
      testData$y = unclass(testData$y)-1
      fit <- glm(y ~ ., data=trainData, family=binomial(link="logit"))
      errorRate <- function(fit,dat) return(sum((predict(fit,dat, type="response")>0.5)==(dat[,"y"]==1))/dim(dat)[1])
    }else if(model == "tree"){
      fit <- rpart(y~.,method="class", data=trainData)
      errorRate <- function(fit,dat) return(sum(predict(fit, newdata = dat, type = "class") == dat[,1])/dim(dat)[1])
    }

    # Train error
    err.train=c(err.train,errorRate(fit,trainData))
    # Test error
    err.test=c(err.test,errorRate(fit,testData))
  }
  # w.mean <- function(X,w) {sum(X*w)/sum(w)}
  # w<-summary(as.factor(batch))

  out <- data.frame(TrainScores = err.train, TestScores = err.test)
  # print(names(batch))
  # print(nrow(out))
  # print(length(names(batch)))
  # cat("\n",paste(names(batch),collapse = " "),"\n")
  rownames(out) = paste(levels(batch))

  # outlst = rep(list(NULL),2);names(outlst)=c("train","test")
  #
  #
  # outlst$train <- list(err.train,
  #                      w.mean(err.train,w),
  #                      min(err.train))
  # names(outlst$train)=c("errors","mean","min")
  # outlst$test <- list(err.test,
  #                     w.mean(err.test,w),
  #                     min(err.test))
  # names(outlst$test)=c("errors","mean","min")

  return(out)
}
