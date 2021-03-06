% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BioMart_input_help.R
\name{BM_attributes}
\alias{BM_attributes}
\title{List valid BioMart attributes}
\usage{
BM_attributes(species)
}
\arguments{
\item{species}{A species covered by BioMart Ensembl. The \code{species} must
be written with the first letter of the genus name followed by the specific
name, e.g. \emph{Homo sapiens} is \dQuote{hsapiens}, and \emph{Mus musculus} is
\dQuote{mmusculus}. Use \code{\link{listSpecies}} to find valid species inputs.}
}
\value{
A table with names and descriptions for all BioMart attributes for a
given \code{species}
}
\description{
\code{BM_attributes} allows users to overview all valid BioMart attributes
for a given species using the \code{\link[biomaRt:listAttributes]{biomaRt::listAttributes}} function
from the \code{biomaRt} package \insertCite{biomart1,biomart2}{GEDI}.
By searching through the available BioMart attributes, the user can find the
BioMart attribute corresponding to the type of reporter IDs used in a dataset.
}
\examples{
## For species: Bos taurus (cow)
# Find all BioMart Attributes for Bos taurus
# if Bos taurus is a valid species in BioMart
species <- "btaurus"
if (species \%in\% listSpecies()){
    BM_attributes(species = species)
}
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{GEDI}} \code{\link{listSpecies}} \code{\link[biomaRt:listDatasets]{biomaRt::listDatasets}} \code{\link[biomaRt]{listAttributes}}
}
