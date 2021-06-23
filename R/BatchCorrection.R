#' @title Remove batch effect from integrated gene expression data
#'
#' @description The \code{BatchCorrection} function removes the batch effect
#' in an integrated gene expression dataset using functions inspired by
#' \code{sva::ComBat} and \code{sva::ComBat_seq} from the \code{sva}
#' package \insertCite{sva}{GEDI}. This function also allows users to visualise the data before and after
#' the batch correction with principal component analysis (PCA) plots and boxplots
#' to verify the batch correction.
#'
#' @name BatchCorrection
#' @usage BatchCorrection(Data, batch, visualise = TRUE, status = NULL)
#'
#' @param Data A \code{data.frame} or \code{matrix} of integrated gene expression
#' data with multiple batches. Columns contain samples, and rows are genes.
#'
#' @param batch Character vector that describes which sample is in what batch.
#' \code{length(batch) = ncol(Data)}.
#'
#' @param visualise (Optional) Logical value. If \code{TRUE}, the batch correction
#' is verified using PCA plots and boxplots that visualise the data before and
#' after the batch correction.
#'
#' @param status (Optional) Character vector that describes the samples' status,
#' e.g. treated or control. \code{BatchCorrection} only uses the \code{status}
#' argument to make plots, i.e. when \code{visualise = TRUE}.
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
#' By setting \code{visualise = TRUE}, the \code{BatchCorrection} function
#' verifies the batch correction visually. This produces a PCA plot and a boxplot
#' for the data before and after the batch correction. If the dataset contains
#' count data, the \code{\link[DESeq2:varianceStabilizingTransformation]{DESeq2::varianceStabilizingTransformation}} function
#' transforms the dataset. This transformation is only used for verification and
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
#' The boxplots show the distribution of gene expression values for all samples.
#' After the batch correction, the distribution should be identical across all
#' batches.
#'
#' @return
#' A \code{data.frame} with the same dimensions as the \code{Data} argument.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point xlab ylab ggtitle theme geom_boxplot
#' aes_string element_blank aes
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation
#' @importFrom SummarizedExperiment assay
#' @importFrom grid viewport grid.newpage pushViewport grid.layout
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


BatchCorrection <- function(Data, batch, visualise = TRUE, status = NULL){
  ### Check inputs ###
  .BatchCorrection.checkInput(Data, batch,visualise, status)
  visualise <- .onlyFirstElement(visualise)

  ### Remove rows with zero variance in any batch ###
  # zero.rows.lst <- lapply(levels(as.factor(batch)), function(batch_level){
  #   if(sum(batch==batch_level)>1){
  #     return(which(apply(Data[, batch==batch_level], 1, function(x){var(x)==0})))
  #   }else{
  #     return(integer(0))
  #   }
  # })
  # zero.rows <- Reduce(union, zero.rows.lst)
  # keep.rows <- setdiff(1:nrow(Data), zero.rows)
  #
  # if (length(zero.rows) > 0) {
  #   cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these are removed from the corrected Data\n", length(zero.rows)))
  #   Data <- Data[keep.rows, ]
  # }

  #Data = data.matrix(Data)
  #counts=all(apply(Data,2,function(x) x%%1)==0)
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

  if(visualise){
    .Visualise(Data,cData,
               batch=batch,
               status=status,
               counts=counts)
  }

  return(cData)
}

.BatchCorrection.checkInput <- function(Data, batch,visualise, status){
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
  .checkClass(visualise, "logical")
  if(!is.null(status)){
    .checkClass(status,c("character","factor","numeric"))
    if(length(status) != ncol(Data)){
      stop(gettextf("'ncol(Data) == length(status)' is not TRUE.
                  'ncol(Data)=%d' and 'length(status)=%d'",
                    ncol(Data),length(status)))
    }
  }
}

.Visualise <- function(Data,cData,batch,status,counts){
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
  box.before <- .box(Data,batch,plot.title = "Boxplot")
  box.after <- .box(cData,batch,plot.title = "Boxplot")

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))

  print(pca.before, vp = vplayout(1,1))
  print(pca.after, vp = vplayout(1,2))
  print(box.before, vp = vplayout(2,1))
  print(box.after, vp = vplayout(2,2))
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
  ggplot(data = suppressMessages(melt(dat)), aes_string(x="variable", y="value")) +
    ggtitle(plot.title) +
    geom_boxplot(aes(fill=rep(batch,each=nrow(dat)))) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())
}


