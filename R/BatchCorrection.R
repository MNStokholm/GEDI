#### Roxygen ####
#' @title Remove batch effect from integrated gene expression data
#'
#' @description The \code{BatchCorrection} function removes the batch effect
#' in an integrated gene expression dataset using functions inspired by
#' \code{sva::ComBat} and \code{sva::ComBat_seq} from the \code{sva}
#' package \insertCite{sva}{GEDI}. This function verifies itself by calculating
#' the mean and standard deviation for each batch before and after the correction.
#' If the means of the batches are similar, then the batch effect has been removed.
#' Optionally, the verification can be supported visually with PCA and RLE plots
#' before and after the batch correction. The RLE plots can be replaced with boxplots.
#'
#' @name BatchCorrection
#' @usage BatchCorrection(Data, batch, visualize = TRUE, status = NULL, boxplot.type = "rle")
#'
#' @param Data A \code{data.frame} or \code{matrix} of integrated gene expression
#' data with multiple batches. Columns contain samples, and rows are genes.
#'
#' @param batch Character vector that describes which sample is in what batch.
#' \code{length(batch) = ncol(Data)}.
#'
#' @param visualize (Optional) Logical value. If \code{TRUE}, the batch correction
#' is verified using PCA plots and boxplots that visualize the data before and
#' after the batch correction.
#'
#' @param status (Optional) Character vector that describes the samples' status,
#' e.g. treated or control. \code{BatchCorrection} only uses the \code{status}
#' argument to make plots, i.e. when \code{visualize = TRUE}.
#'
#' @param boxplot.type (optional) Character that describes which type of boxplot
#' to use if \code{visualize=TRUE}. Default is \code{"rle"} resulting in a RLE plot.
#' \code{"normal"} results in a regular boxplot.
#'
#'
#' @details
#' Bioconductor's \code{sva} package does not remove batch effect in genes with
#' zero variance or zero counts within any batch. Therefore, the \code{BatchCorrection}
#' function uses modified versions of the \code{sva::ComBat} and
#' \code{sva::ComBat_seq} functions to remove the batch effect in all genes.
#' The modified \code{ComBat_seq} function corrects data for the batch effect
#' like the modified \code{ComBat} function, but only for RNA-seq count data.
#'
#' The \code{BatchCorrection} function verifies the batch correction by calculating
#' and printing the mean and the standard deviation of the gene expressions in each batch.
#' If the means of the batches are similar, then the batch effect has been removed.
#' By setting \code{visualize = TRUE}, the \code{BatchCorrection} function
#' supports the verification visually. This produces a PCA plot and an RLE plot
#' for the data before and after the batch correction. If the dataset contains
#' count data, the \code{\link[DESeq2:varianceStabilizingTransformation]{DESeq2::varianceStabilizingTransformation}} function
#' transforms the dataset first. This transformation is only used for verification and
#' is not permanent.
#'
#' The PCA plot is a scatterplot where the x- and y-axis are the first two principal
#' components. Initially, the samples aggregate in batches. After batch correction,
#' the samples should no longer aggregate in batches, but instead, the batches
#' should lie on top of each other. Suppose it was possible to distinguish between
#' samples based on some other variable within each batch before the batch effect was
#' removed. In that case, it must still be possible to differentiate the samples
#' after the batch correction.
#'
#' The RLE plots are boxplots with log-transformed data subtracted the median of
#' each row. This plot shows the distribution of gene expression values for all samples.
#' After the batch correction, the distribution should be similar across all
#' batches. Normal boxplots shows the same thing, but unlike RLE, the data is
#' not transformed.
#'
#' @return
#' A \code{data.frame} with the same dimensions as the \code{Data} argument.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point xlab ylab ggtitle theme geom_boxplot
#' aes_string element_blank aes coord_cartesian
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation
#' @importFrom SummarizedExperiment assay
#' @importFrom grid viewport grid.newpage pushViewport grid.layout
#' @importFrom grDevices boxplot.stats
#' @importFrom stats median
#' @seealso
#' \code{\link{ReadGE}}
#' \code{\link{GEDI}}
#' \code{\link{VerifyGEDI}}
#' @examples
#' \dontrun{
#' ### From example in README.md
#' ## Continuing example from GEDI()'s example, see help("GEDI")
#' # Each dataset is a batch
#' # Each sample has a status of 'positive' or 'negative'
#' # Batch and Status variables is seen below:
#' summary(as.factor(batch))
#' #> B1 B2 B3
#' #>  8  6  6
#' summary(as.factor(status))
#' #> control treated
#' #>      10      10
#'
#' # The batch effect is removed
#' cData <- BatchCorrection(dat, batch, status = status)
#' #> Found3batches
#' #> Adjusting for0covariate(s) or covariate level(s)
#' #> Standardizing Data across genes
#' #> Fitting L/S model and finding priors
#' #> Finding parametric adjustments
#' #> Adjusting the Data
#'
#' # To see the resulting figure, check the README.md file in the
#' # package GitHub repository: github.com/s184257/GEDI
#' }
#' @export

#### Code ####
BatchCorrection <- function(Data, batch, visualize = TRUE, status = NULL,
                            boxplot.type = "rle"){
  ### Check inputs ###
  .BatchCorrection.checkInput(Data, batch,visualize, status, boxplot.type)
  visualize <- .onlyFirstElement(visualize)

  # Check if 'Data' contains count data
  counts = all(apply(Data,2,class)=="integer")

  if (counts){
    Data = data.matrix(Data)
    cData = .myComBat_seq(counts=Data, batch=batch, group = NULL)
  }else{
    modcombat = model.matrix(~1,data=data.frame(batch))#, data=pheno
    cData = .myComBat(dat=Data, batch=batch, mod=modcombat,
                               par.prior=TRUE, prior.plots=FALSE)
  }
  cData = as.data.frame(cData)

  ## Calculate mean (M) and standard deviation (SD) for each group
  ## before and after batch correction
  # Before
  df <- data.frame(GeneExp = apply(Data,2,mean), Group = batch)
  # After
  cdf <- data.frame(GeneExp = apply(cData,2,mean), Group = batch)
  # Create and print table
  res <- rbind(
    with(df, tapply(GeneExp, Group, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))})),
    with(cdf, tapply(GeneExp, Group, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))}))
  )
  rownames(res) <- c("Before", "After")
  print(t(res))

  ## Make visualisations
  if(visualize){
    .visualize(Data,cData,
               batch=batch,
               status=status,
               counts=counts,
               boxplot.type = boxplot.type)
  }

  return(cData)
}

.BatchCorrection.checkInput <- function(Data, batch,visualize, status, boxplot.type){
  ## Data ##
  .checkClass(Data,c("data.frame","matrix","array"))
  if(any(dim(Data)==0)){
    stop(gettextf("'Data' is empty"))
  }else if(any(!apply(Data,2,class) %in% c("numeric","integer"))){
    colError <- match(TRUE,!apply(Data,2,class) %in% c("numeric","integer"))
    stop(gettextf("All columns in 'Data' must all be of class %s or %s. Column %d has class %s",
                  dQuote("numeric"),dQuote("integer"),colError,dQuote(class(Data[,colError]))))
  }
  ## batch ##
  .checkClass(batch,c("numeric","character","factor"))
  if (ncol(Data) != length(batch)){
    stop(gettextf("'ncol(Data) == length(batch)' is not TRUE.
                  'ncol(Data)=%d' and 'length(batch)=%d'",
                  ncol(Data),length(batch)))
  }
  ## Validation plot ##
  .checkClass(visualize, "logical")
  if(!is.null(status)){
    .checkClass(status,c("character","factor","numeric"))
    if(length(status) != ncol(Data)){
      stop(gettextf("'ncol(Data) == length(status)' is not TRUE.
                  'ncol(Data)=%d' and 'length(status)=%d'",
                    ncol(Data),length(status)))
    }
  }
  .checkClass(boxplot.type, "character")
  validType=c("normal","rle")
  if(!tolower(boxplot.type) %in% validType){
    n.boxplot.type = length(validType)
    boxplot.type.msg = paste(paste(dQuote(validType[-n.boxplot.type]),collapse=", "),
                      "and",dQuote(validType[n.boxplot.type]),sep=" ")
    stop(gettextf("Invalid 'boxplot.type' chosen. Valid options are: %s",boxplot.type.msg))
  }
}

.visualize <- function(Data,cData,batch,status,counts,boxplot.type){
  Data = as.data.frame(Data)
  if(counts){
    coldata <- data.frame(rep(NA,ncol(Data)))
    # print(nrow(coldata));print(ncol(Data));print(ncol(cData))
    dds <- DESeqDataSetFromMatrix(Data, coldata, design = ~1)
    vsd <- varianceStabilizingTransformation(dds, blind=TRUE, fitType="local")
    Data <- as.data.frame(assay(vsd))

    cdds <- suppressMessages(DESeqDataSetFromMatrix(cData, coldata, design = ~1))
    cvsd <- varianceStabilizingTransformation(cdds, blind=TRUE, fitType="local")
    cData <- as.data.frame(assay(cvsd))
  }

  pca.before <- .PCAplot(Data,batch=batch,status = NULL,
                         plot.title = "Before batch correction\nPCA plot")
  pca.after <- .PCAplot(cData,batch=batch,status = status,legend=TRUE,
                        plot.title = "After batch correction\nPCA plot")
  if (boxplot.type == "rle"){
    boxplot.before <- .rle(Data,batch,plot.title = "RLE plot")
    boxplot.after <- .rle(cData,batch,plot.title = "RLE plot")
  }else if (boxplot.type == "normal"){
    boxplot.before <- .box(Data,batch,plot.title = "Boxplot")
    boxplot.after <- .box(cData,batch,plot.title = "Boxplot")
  }



  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))

  print(pca.before, vp = vplayout(1,1))
  print(pca.after, vp = vplayout(1,2))
  print(boxplot.before, vp = vplayout(2,1))
  print(boxplot.after, vp = vplayout(2,2))
}

.PCAplot <- function(dat,batch,status,legend=FALSE,plot.title=""){
  # calculate the variance for each gene
  ntop=500
  rv <- apply(dat,1,var)

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(dat[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Batch=batch, name=colnames(dat))
  if(is.null(status)){
    Aes <- aes_string(x="PC1", y="PC2",color="Batch")
  }else{
    d<-data.frame(d,Status=status)
    Aes <- aes_string(x="PC1", y="PC2",color="Batch", shape="Status")
  }

  legend.pos <- if(legend) {"right"} else {"none"}


  ggplot(data=d) +
    Aes +
    geom_point(size=3) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    ggtitle(plot.title) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = legend.pos)
}

.box <- function(dat,batch,plot.title){
  df <- melt(dat, measure = 1:dim(dat)[2])
  p0 <- ggplot(data = df, aes_string(x="variable", y="value")) +
    ggtitle(plot.title) +
    geom_boxplot(aes(fill=rep(batch,each=nrow(dat))), outlier.shape = NA) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())
  # compute lower and upper whiskers
  ylim1 = boxplot.stats(df$value)$stats[c(1, 5)]

  # scale y limits based on ylim1
  p0 + ggplot2::coord_cartesian(ylim = ylim1*1.05)
}

.rle <- function(dat,batch,plot.title){
  logDat <- log(dat+1)
  med <- apply(logDat,1,median)
  rle <- as.data.frame(apply(logDat, 2, function(x) x - med))

  df <- melt(rle, measure = 1:dim(rle)[2])
  p0 <- ggplot(data = df, aes_string(x="variable", y="value")) +
    ggtitle(plot.title) +
    geom_boxplot(aes(fill=rep(batch,each=nrow(dat))), outlier.shape = NA) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())
  # compute lower and upper whiskers
  ylim1 = boxplot.stats(df$value)$stats[c(1, 5)]

  # scale y limits based on ylim1
  p0 + ggplot2::coord_cartesian(ylim = ylim1*1.05)
}

