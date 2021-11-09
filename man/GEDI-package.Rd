\name{GEDI-package}
\alias{GEDI-package}
\docType{package}
\title{
\packageTitle{GEDI}
}
\description{
\packageDescription{GEDI}
}
\details{
The main functions are:
\itemize{
\item \code{\link{ReadGE}} - reads all the gene expression datasets and stores them in a list
\item \code{\link{GEDI}} - maps probe or read IDs to Ensembl or Entrez gene IDs using the \code{biomaRt} package and integrates the datasets into one data table
\item \code{\link{BatchCorrection}} - removes the batch effect from the integrated dataset
\item \code{\link{VerifyGEDI}} - verifies the transcriptomic data integration using one of two supervised machine learning models to predict the samples' status in one batch based on the remaining batches' samples
}

The package source files as well as a README file containing an example pipeline can be found in the github repository:
\url{www.github.com/MNStokholm/GEDI}.

}
\author{
\packageAuthor{GEDI}

Maintainer: \packageMaintainer{GEDI}
}

\keyword{ package }
