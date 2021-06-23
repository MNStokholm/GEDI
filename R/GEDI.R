#### Info ####
#' @title Re-annotate and integrate gene expression datasets
#'
#' @description \code{GEDI} allows users to integrate transcriptomic datasets
#' into one data table containing samples from all the datasets. The first step
#' is to map the reporter IDs to Ensembl gene IDs. The second step is to collapse
#' rows mapped to the same Ensembl gene ID. The third step is to merge the
#' datasets. Optionally, \code{GEDI} can fill in missing data using a linear
#' regression model.
#'
#' @name GEDI
#' @usage
#' GEDI(datasets, attributes, species = "hsapiens", BioMart = TRUE,
#'      collapseMethod = "MaxMean", predict.missing = FALSE,
#'      predict.threshold = 0.8, path=NULL)
#'
#' @param datasets A list of transcriptomic datasets. The \code{\link{ReadGE}}
#' function can create this list. The \code{GEDI} function integrates all datasets
#' in the list.
#'
#' @param attributes Character vector of BioMart attributes that correspond to
#' each datasets reporter IDs. The \code{\link{BM_attributes}} function can find
#' all available BioMart attributes for a given \code{species}.
#'
#' @param species (Optional) Character that refers to a species covered by BioMart
#' Ensembl. The \code{species} must be written with the first letter of the genus
#' name followed by the specific name, e.g. \emph{Homo sapiens} is \dQuote{hsapiens},
#' and \emph{Mus musculus} is \dQuote{mmusculus}. The \code{\link{listSpecies}} function
#' can find all valid species inputs.
#'
#' @param BioMart (Optional) Logical value. When \code{TRUE} (default), the
#' reporter IDs are mapped to Ensembl gene IDs using the \code{biomaRt} package.
#' If \code{FALSE}, the \code{GEDI} function maps reporter IDs to Ensembl IDs
#' using a tab-separated annotation table named \dQuote{annot.txt} for each dataset.
#' The reporter IDs must be in the first column of the table. The second column holds
#' Ensembl gene IDs. The second columns name must correspond to the BioMart attribute,
#' e.g. \dQuote{ensembl_gene_id}.
#'
#' @param collapseMethod (Optional) Character that determines which method the
#' \code{\link[WGCNA:collapseRows]{WGCNA::collapseRows}} function uses to collapse two rows
#' mapped to the same Ensembl gene ID. Valid options are: \code{MaxMean} (default),
#' \code{MinMean}, \code{absMaxMean}, \code{absMinMean}, \code{Average},
#' \code{maxRowVariance} and \code{ME}. See \code{\link[WGCNA:collapseRows]{WGCNA::collapseRows}}
#' documentation for more info.
#'
#' @param predict.missing (Optional) Logical value. When \code{TRUE}, \code{GEDI}
#' predicts missing expression values for a gene if the percentage of known values
#' is higher than the \code{predict.threshold}. This does not currently work for
#' count data.
#'
#' @param predict.threshold (Optional) Numeric argument between 0 and 1. This
#' argument is ignored if \code{predict.missing = FALSE}.
#'
#' @param path (Optional) (Optional) Path to the directory where all data folders
#' are located. The data folders store annotation files if there are any.
#' The current working directory is the default. The names of the data folders
#' are \code{names(datasets)}.
#'
#'
#' @details
#' The \code{GEDI} function maps reporter IDs of the datasets to Ensembl gene IDs
#' using the \code{biomaRt} package \insertCite{biomart1,biomart2}{GEDI}.
#' However, if the BioMart attribute for a dataset is invalid or missing, an
#' annotation table called \dQuote{annot.txt} is needed;
#' see the \code{BioMart} argument. If the second column does not contain Ensembl
#' gene IDs, then it must contain another type of reporter ID with a known BioMart
#' attribute.
#' The annotation file must be in the same subdirectory as the corresponding
#' dataset. The subdirectory's name must correspond to the dataset's name; see
#' \code{names(datasets)}. The annotation file is also used if the user chooses
#' not to use BioMart.
#'
#' Multiple reporter IDs can map to the same Ensembl gene ID. Therefore, these
#' reporter IDs' expression values must be collapsed using the
#' \code{\link[WGCNA:collapseRows]{WGCNA::collapseRows}} function. The method used to collapse
#' the rows can be specified with the \code{collapseMethod} argument. After this,
#' the datasets are merged.
#'
#' If \code{predict.missing = TRUE}, then the missing expression values are
#' predicted for all Ensembl gene IDs where the percentage of known expression
#' values exceeds the \code{predict.threshold} argument. \code{GEDI} uses a
#' linear regression model to predict the missing expression values. The linear
#' regression model uses only one feature; the gene most correlated with the gene
#' with missing values. Only one feature is chosen to get a quick estimate.
#' However, it can still be time-consuming if many rows have missing values.
#'
#' @return
#' A \code{Data.frame} containing columns from all datasets. All rows are mapped
#' to Ensembl gene IDs. Rows with missing data have been removed from the output.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom progress progress_bar
#' @importFrom WGCNA collapseRows
#' @export
#' @seealso
#' \code{\link{ReadGE}}
#' \code{\link{BatchCorrection}}
#' \code{\link{VerifyGEDI}}
#' @examples
#' \dontrun{
#' ### From example in README.md
#' ## Continuing example from ReadGE()'s example, see help("ReadGE")
#' # Biomart attributes defined for the three datasets
#' # The species is Bos taurus
#' attr <- c("ensembl_gene_id", affy_bovine, NA)
#'
#' # There are no attribute for the third dataset
#' # An annotation file (annot.txt) are used
#' # It is located in the datafolder for the dataset
#'
#' dat <- GEDI(datasets, attributes = attr, BioMart = TRUE,
#'            species = "btaurus", path = PATH_TO_DATA_FOLDERS)
#'            #> Connecting to BioMart...
#'            #> Invalid BioMart attribute, NA in agilent_data
#'            #> Reading annotation file ('annot.txt') for agilent_data
#' dim(dat)
#' #> [1] 8328   20
#'
#' # The integrated dataset has 20 samples with expression values
#' # for 8328 genes
#' }


#### GEDI ####
GEDI <- function(datasets, attributes,
                 species = "hsapiens",
                 BioMart = TRUE,
                 collapseMethod = "MaxMean",
                 predict.missing = FALSE, predict.threshold = 0.8,path=NULL){
  ### Check input ###
  if(is.null(path)) path="."
  if(class(collapseMethod)=="character") collapseMethod <- .onlyFirstElement(collapseMethod)
  .GEDI.checkInput(datasets, attributes,
                   species = species,
                   BioMart = BioMart,
                   collapseMethod = collapseMethod,
                   predict.missing = predict.missing,
                   predict.threshold = predict.threshold, path=path)

  species <- .onlyFirstElement(species)
  BioMart <- .onlyFirstElement(BioMart)
  predict.missing <- .onlyFirstElement(predict.missing)
  predict.threshold <- .onlyFirstElement(predict.threshold)
  path <- .onlyFirstElement(path)



  ## BioMart ##
  if (BioMart==TRUE & any(attributes != "ensembl_gene_id")){
    cat("Connecting to BioMart...\n")
    if(curl::has_internet()){
      if (species %in% listSpecies()$species){
        mart <- useDataset(paste(species,"_gene_ensembl",sep=""), useMart("ensembl"))
      }else{
        stop(gettextf("Invalid argument for %s: %s. See all species covered by BioMart and ensembl using the listSpecies() function.",
                      sQuote("species"),dQuote(species)))
      }
    }else{
      stop("Could not connect to BioMart. Check your internet connection")
    }

  }else mart=NULL

  ### Initialization ###
  ## Progress messages ##
  n.steps = length(datasets)*3+1
  pb<-progress_bar$new(total = n.steps,
                       format = ":percent\t:msg",clear=TRUE)
  pb$tick(0);Sys.sleep(0.2)
  ## Variables ##
  attributes[attributes==""] <- NA
  total <- NULL
  Data.names <- names(datasets)
  for (i in 1:length(datasets)){
    working.on = paste("Working on ",Data.names[i],":",sep="")
    pb$tick(0,tokens=list(msg=working.on))

    # Extract data set
    Data <- datasets[[i]]
    # Check if Data is correct
    if (class(Data) != "data.frame"){
      msg=gettextf("Invalid class for the %s dataset in 'datasets'. Must be of class 'data.frame' and not %s",
                   dQuote(Data.names[i]),dQuote(class(Data)))
      pb$message(msg);warning(msg)
      pb$message(paste("Skipping dataset:",Data.names[i]))
      next
    }else if (any(dim(Data)==0)){
      msg=gettextf("Dataset %s in 'datasets' is empty", dQuote(Data.names[i]))
      pb$message(msg);warning(msg)
      pb$message(paste("Skipping dataset:",Data.names[i]))
      next
    }else if(any(!apply(Data,2,class)%in% c("numeric","integer"))){
      colError <- match(TRUE,!apply(Data,2,class) %in% c("numeric","integer"))
      pb$message(paste("All columns in '",Data.names[i],"' must be of class 'numeric' or 'integer'."))
      pb$message(paste("Column ",colError," has class '",class(Data[,colError]),"'.",sep=""))
      pb$message(paste("Skipping dataset:",Data.names[i]))
      next
    }

    ID <- attributes[i]
    if (.NAtoFALSE(ID != "ensembl_gene_id") | is.na(ID)){
      pb$update(((i-1)*3+1)/n.steps,
                tokens=list(msg=paste(working.on,"Re-annotating data. ")))
      Data <- .Map2ensembl(Data, ID, Data.names[i],
                          species = species,
                          BioMart = BioMart,
                          mart = mart, path = path,
                          progress.bar = pb, progress.msg = working.on)
      if(is.null(Data)) {
        pb$message(paste("Skipping dataset:",Data.names[i]))
        next
      }
      pb$update(((i-1)*3+2)/n.steps,
                tokens=list(msg=paste(working.on,"Collapsing rows. ")))
      Data <- .CollapseRows(Data, collapseMethod)
    }
    pb$update(i*3/n.steps,
              tokens=list(msg=paste(working.on,"Merging. ")))
    if (is.null(total)){
      total <- Data
    }else{
      total <- merge(total,Data,by = "row.names",all.x=TRUE)
      rowID <- total$Row.names
      total <- total[,-c(1)]
      rownames(total) <- rowID
    }
  }
  pb$update(1,tokens = list(msg=""))

  if (predict.missing){
    cat("___________________________________\n")
    cat("Predicting missing values using linear regression...\n")
    cat(paste("Can only predict values for genes where ",
              predict.threshold*100,"% of the values are known.\n",sep=""))
    total <- .EstMissing(total,threshold = predict.threshold)
  }

  if(!is.null(total)) total <- total[complete.cases(total),]
  return(total)
}

#### .GEDI.checkInput ####
.GEDI.checkInput <- function(datasets, attributes, species, BioMart, collapseMethod,
                             predict.missing, predict.threshold, path){
  ## datasets ##
  .checkClass(datasets,"list")
  if (is.null(names(datasets))) {
    names(datasets)=c(1:length(datasets))
  }
  ## attributes ##
  .checkClass(attributes,"character")
  if (length(datasets) != length(attributes)){
    stop(gettextf("'length(datasets) == length(attributes)' is not TRUE.
                  'length(datasets)=%d' and 'length(attributes)=%d'",
                  length(datasets),length(attributes)))
  }
  ## species ##
  .checkClass(species,"character")

  ## collapseMethod ##
  validCollapseMethods = c("maxRowVariance", "MaxMean",
                           "MinMean", "absMaxMean", "absMinMean", "ME","Average")
  if (!collapseMethod %in% validCollapseMethods){
    n.methods = length(validCollapseMethods)
    collapseMethod.msg = paste(paste(dQuote(validCollapseMethods[-n.methods]),collapse=", "),
                               "and",dQuote(validCollapseMethods[n.methods]),sep=" ")
    stop(gettextf("Invalid 'collapseMethod' chosen. Valid options are: %s",collapseMethod.msg))
  }

  .checkClass(path,"character")
  .checkClass(BioMart,"logical")

  ## predict.missing ##
  .checkClass(predict.missing,"logical")
  ## predict.threshold ##
  .checkClass(predict.threshold,"numeric")
  if (predict.threshold <= 0 | 1 <= predict.threshold){
    stop("'predict.threshold' out of range. Must be a value between 0 and 1.")
  }
}
#### .Map2ensembl ####
.Map2ensembl <- function(Data, ID, Data.name,
                         species = "hsapiens",
                         BioMart = TRUE,
                         mart = NULL,path=".", progress.bar, progress.msg){
  annot <- .getAnnot(Data, ID, Data.name,
            species = species,
            BioMart = BioMart,
            mart = mart, path = path,
            progress.bar = progress.bar, progress.msg = progress.msg)

  if(is.null(annot)){
    # Message that explains why annot is NULL is provided by .getAnnot
    return(NULL)
  }

  Data <- merge(annot, Data, by.x=ID, by.y="row.names")
  rowID <- Data[,1]
  Data <- Data[,-c(1)]
  rownames(Data) <- rowID
  Data[Data==""] <- NA
  Data <- Data[complete.cases(Data),]
  return(Data)
}

#### .getAnnot ####
.getAnnot <- function(Data, ID, Data.name,
                      species = "hsapiens",
                      BioMart = TRUE,
                      mart = NULL,path=".",progress.bar,progress.msg){
  if (BioMart){
    if (.NAtoFALSE(ID %in% BM_attributes(species)$name)){
      # BioMart is used to make annotation table
      genes <- rownames(Data)
      annot <- getBM(filters= ID, attributes= c("ensembl_gene_id", ID),values=genes,mart= mart)
      annot = annot[!duplicated(annot[,ID]),]
      return(annot)
    }else{
      progress.bar$message(paste("Invalid BioMart attribute,",
                                 ID,"in",Data.name))
      progress.bar$message(paste("Reading annotation file ('annot.txt') for",Data.name))
    }
  }

  progress.bar$tick(0, tokens=list(msg=paste(progress.msg,
                                             "Re-annotating data...",
                                             "Reading annotation file. ")))
  annot.path = paste(path,"/",Data.name,"/annot.txt",sep="")
  if (file.exists(annot.path)) {
    annot <- read.table(annot.path,header = TRUE)
  } else {
    progress.bar$message(paste("Annotation file not found at location:",annot.path))
    return(NULL)
  }
  colnames(annot)[1] <- ID
  ID2 <- colnames(annot)[2]

  if (.NAtoFALSE(ID2 != "ensembl_gene_id") | is.na(ID2)){
    if (BioMart){
      if (.NAtoFALSE(ID2 %in% BM_attributes(species)$name)){
        annot[annot==""] <- NA
        annot <- annot[complete.cases(annot),]

        # Map ID2 to 'ensembl_gene_id'
        genes <- annot[,ID2]
        G_list <- getBM(filters= ID2, attributes= c("ensembl_gene_id", ID2),values=genes,mart= mart)
        G_list = G_list[!duplicated(G_list[,ID2]),]
        annot <- merge(annot,G_list,by=ID2)
        annot <- annot[,-c(1)]
      }else{
        progress.bar$message(paste("Invalid BioMart attribute,'",ID2,"',in second column of the annotation file:",annot.path))
        progress.bar$message(paste("Unable to re-annotate dataset,",Data.name))
        return(NULL)
      }
    }else{
      progress.bar$message("Second column of annotation file is not 'ensembl_gene_id'. Change annotation file or set BioMart=TRUE")
      return(NULL)
    }
  }
  return(annot)
}

#### .CollapseRows ####
.CollapseRows <- function(Data,method){
  datET <- Data[,-c(1)]
  rowID <- rownames(Data)
  rowGroup <- Data[,1]
  rownames(datET) <- rowID

  collapse = collapseRows(datET, rowGroup, rowID,
                          method=method, connectivityBasedCollapsing=FALSE,
                          methodFunction=NULL, connectivityPower=1,
                          selectFewestMissing=TRUE, thresholdCombine=NA)

  Data <- collapse$datETcollapsed
  return(Data)
}
#### .EstMissing ####
.EstMissing <- function(Data,threshold){
  canFill.idx = c()
  for(i in 1:dim(Data)[1]){
    row <- Data[i,]
    if(.canFill(row,threshold)) {canFill.idx=c(canFill.idx, i)}
  }
  cat(length(canFill.idx),"genes fulfil the criteria.\n")
  # Only genes with no missing values
  complete.data <- Data[complete.cases(Data),]

  pb <- progress_bar$new(total = length(canFill.idx),
                         format = "Predicting missing values [:bar] :percent eta: :eta")

  pb$tick(0)
  on.exit(pb$terminate())

  for(i in canFill.idx){
    pb$tick()
    row <- Data[i,]
    misCol <- is.na(row) %in% TRUE
    # Train
    trainData <- complete.data[,!misCol];trainData=data.frame(t(trainData))
    Ytrain <- t(row[!misCol]);colnames(Ytrain)=c("Y")
    Xidx = order(apply(trainData,2,cor,Ytrain)**2, decreasing = TRUE)[1]
    Xtrain <- trainData[,Xidx]
    trainData = data.frame(Ytrain,Xtrain);

    Xtest <- complete.data[Xidx,misCol]
    Xtest=data.frame(t(Xtest));
    colnames(Xtest)=colnames(trainData)[-1]
    #Xtest <- Xtest[,Xidx]
    if(class(Xtest)=="numeric"){
      Xtest <- data.frame(Xtest)
      colnames(Xtest) <- colnames(trainData)[-1]
    }
    fit <- lm(Y~.,data = trainData)
    Data[i,misCol] <- predict(fit,Xtest)
  }
  return(Data)
}



