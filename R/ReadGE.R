#' @title Read multiple gene expression datasets
#'
#' @description \code{ReadGE} allows users to read multiple transcriptomic
#' datasets and store them in a list of \code{data.frames}. A subdirectory
#' for each dataset contains the data files. The function can read Affymetrix,
#' Agilent and RNA sequencing data.
#' @name ReadGE
#' @usage ReadGE(dataFolders, sources, path = NULL, verbose = TRUE)
#'
#' @param dataFolders Character vector with the names of the subdirectories
#' containing gene expression datasets.
#' @param sources Character vector with the data source for each dataset. Valid
#' sources are \dQuote{affy} or \dQuote{affymetrix}, \dQuote{agilent} and
#' \dQuote{RNAseq}.
#' @param path (Optional) Logical value. If \code{TRUE}, the function prints
#' descriptive messages about the progress.
#' @param verbose (Optional) Logical value. If \code{TRUE}, the function prints
#' descriptive messages about the progress.
#'
#' @details
#' \code{ReadGE} reads datasets with different methods depending on the source.
#' Sample files for Agilent and Affymetrix data are read in alphabetical order.
#' Samples are presented in the same order on Gene Expression Omnibus (GEO);
#' \url{https://www.ncbi.nlm.nih.gov/geo/}.
#'
#' \dQuote{affy}/\dQuote{affymetrix} corresponds to Affymetrix data. \code{ReadGE}
#' reads Affymetrix data using the \cr\code{\link[affy:ReadAffy]{affy::ReadAffy}} function and
#' converts it to an expression data table using the \code{\link[gcrma:gcrma]{gcrma::gcrma}}
#' function. Affymetrix data files have the .CEL or .CEL.gz file extension.
#'
#' \dQuote{agilent} corresponds to Agilent data. \code{ReadGE} reads Agilent data
#' using the \code{\link[limma:read.maimages]{limma::read.maimages}} function. Agilent sample files
#' have the file extension .txt or .txt.gz. The filenames \dQuote{targets.txt} and
#' \dQuote{annot.txt} are ignored.
#'
#' \dQuote{RNAseq} corresponds to RNA sequencing data. \code{ReadGE} reads
#' RNA-seq data using the \code{DESeq2} package. \dQuote{RawCounts.txt} is a
#' tab-separated table containing RNA-seq count data. This file must contain raw
#' counts and not transformed or normalised data. \code{ReadGE} transforms the
#' raw count data using the \code{\link[DESeq2:varianceStabilizingTransformation]{DESeq2::varianceStabilizingTransformation}}
#' function so the RNA sequencing data can be integrated with microarray data.
#' \code{ReadGE} does not transform the counts if all datasets are \dQuote{RNAseq}.
#'
#' @return
#' A list of all transcriptomic datasets that were successfully read where
#' \code{names(datasets) = dataFolders}.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix counts varianceStabilizingTransformation
#' @importFrom limma read.maimages backgroundCorrect normalizeBetweenArrays avereps
#' @importFrom affy ReadAffy
#' @importFrom gcrma gcrma
#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assay
#' @seealso
#' \code{\link{GEDI}}
#' \code{\link{BatchCorrection}}
#' \code{\link{VerifyGEDI}}
#' @export
#' @examples
#' \dontrun{
#' ### From example in README.md
#' ## Three transcriptomic datasets are read
#' library(GEDI)
#' # Names of folders with data
#' dataFolders <- c("RNAseq_data", "affy_data", "agilent_data")
#' # Which type of data is each dataset
#' sources <- c("RNAseq", "affy", "agilent")
#'
#' # Read the data
#' datasets <- ReadGE(dataFolders, sources,
#'                    path = "PATH_TO_DATA_FOLDERS")
#' #> Reading RNAseq_data...
#' #> Reading affy_data...
#' #> Reading agilent_data...
#' }


ReadGE <- function(dataFolders, sources, path = NULL, verbose = TRUE){
  ## Check input ##
  if(is.null(path)) path="."
  .ReadGE.checkInputs(dataFolders,sources,path,verbose)
  verbose <- .onlyFirstElement(verbose)

  ## Define other variables ##
  # Do not transform counts if all sources are RNAseq
  transformCounts <- !all(sources == "RNAseq")
  n = length(dataFolders) # Number of datasets


  # Define empty data collection to store data sets in
  collection <- rep(list(NULL), n)
  names(collection) <- dataFolders

  for (i in 1:n){
    dataFolder <- dataFolders[i]
    dataFolder.path = paste(path,"/",dataFolder,sep="")
    source <- sources[i]
    if(verbose) cat(paste("Reading ",dataFolder,"...\n",sep=""))

    if (!file.exists(dataFolder.path)){
      message(gettextf("Data folder, %s, not found", dQuote(dataFolder)))
      message(gettextf("Skipping dataset: %s", dQuote(dataFolder)))
      next
    }

    ## Affymetrix data ##
    if (source %in% c("affy","affymetrix")){
      if(length(list.files(path = dataFolder.path, pattern = ".CEL"))==0){
        message(gettextf("No affymetrix data (.CEL files) for %s", dQuote(dataFolder)))
        message(gettextf("Skipping dataset: %s", dQuote(dataFolder)))
        next
      }
      Data <- ReadAffy(celfile.path = dataFolder.path)
      # Be careful with this
      eset <- suppressPackageStartupMessages(suppressWarnings(gcrma(Data, verbose = FALSE)))
      #Calculating Expression
      Data <- data.frame(exprs(eset))

      collection[dataFolder] <- list(Data)
    }
    ## Agilent ##
    else if (source == "agilent"){
      # Create list of targets
      targets <- list.files(path = dataFolder.path, pattern = ".txt")
      targets <- targets[!targets%in%c("annot.txt","targets.txt")]
      if(length(targets)==0){
        message(gettextf("No agilent data for %s", dQuote(dataFolder)))
        message(gettextf("Skipping dataset: %s", dQuote(dataFolder)))
        next
      }

      # Read data
      x <- read.maimages(targets, source="agilent", green.only=TRUE,
                         path = dataFolder.path, verbose = FALSE)
      y <- backgroundCorrect(x, method="minimum")
      y <- normalizeBetweenArrays(y, method="quantile")
      neg95 <- apply(y$E[y$genes$ControlType==-1,], 2, function(x) quantile(x,p=0.95))
      cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
      isexpr <- rowSums(y$E > cutoff) >= 4 ##This number should be equal than the number of replicates in each grou
      y0 <- y[y$genes$ControlType==0 & isexpr,]

      Data <- avereps(y0,ID=y0$genes[,"ProbeName"])
      rowID <- Data$genes$ProbeName
      Data <- data.frame(Data$E)
      rownames(Data) <- rowID

      collection[dataFolder] <- list(Data)
    }
    ## RNA-seq ##
    else if (source == "RNAseq"){
      if(!file.exists(paste(dataFolder.path,"/RawCounts.txt",sep=""))){
        message(gettextf("The RNA-sequencing data file, %s was not found in %s", dQuote("RawCounts.txt"), dQuote(dataFolder)))
        message(gettextf("Skipping dataset: %s", dQuote(dataFolder)))
        next
      }
      countTable = read.table(paste(dataFolder.path,"/RawCounts.txt",sep=""), header=TRUE, row.names=1)

      # Make arbitrary design file
      coldesign = data.frame(rep(dataFolder,dim(countTable)[2]))
      rownames(coldesign) <- colnames(countTable)
      colnames(coldesign) <- c("Folder.name")

      dds <- DESeqDataSetFromMatrix(countData = countTable, colData = coldesign, design= ~1)

      if (transformCounts){
        # There are other sources than 'RNAseq'
        # Therefore the data is transformed
        dds <- varianceStabilizingTransformation(dds, blind=T, fitType="local")
      }

      Data <-as.data.frame(assay(dds))
      collection[dataFolder] <- list(Data)
    }
    else{
      message(gettextf("Unknown source for %s: %s", dQuote(dataFolder), dQuote(source)))
      message(gettextf("Valid sources are %s or %s, %s and %s",
                       dQuote("affy"), dQuote("affymetrix"),
                       dQuote("agilent"), dQuote("RNAseq")))
      message(gettextf("Skipping dataset: %s", dQuote(dataFolder)))
    }

    if(is.null(collection[[dataFolder]])) collection[[dataFolder]] = NULL
  }
  return(collection)
}

.ReadGE.checkInputs <- function(dataFolders, sources, path, verbose){
  .checkClass(dataFolders,"character")
  if (!all(is.na(sources))){
    .checkClass(sources,"character")
  }
  if(length(dataFolders) != length(sources)){
    stop(gettextf("'length(dataFolders) == length(sources)' is not TRUE.
                  'length(dataFolders)=%d' and 'length(sources)=%d'",
                  length(dataFolders),length(sources)))
  }

  .checkClass(verbose,"logical")
  if(class(path) != "character"){
    stop(gettextf("Invalid class for %s argument: %s. 'class(path)' must return %s",
                  sQuote("path"), dQuote(class(path)), dQuote("character")))
  }
}

