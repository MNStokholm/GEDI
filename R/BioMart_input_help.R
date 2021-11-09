#' @title List valid species for the biomaRt package
#'
#' @description
#' \code{listSpecies} allows the user to overview all available species for
#' BioMart Ensembl using the \code{\link[biomaRt:listDatasets]{biomaRt::listDatasets}} function from
#' the \code{biomaRt} package \insertCite{biomart1,biomart2}{GEDI}. \code{listSpecies} makes it easy to find the
#' correct species input for the \code{\link{GEDI}} function.
#'
#' @name listSpecies
#' @usage
#' listSpecies()
#'
#' @return A table with names and descriptions of all the valid species
#'
#' @seealso \code{\link{GEDI}} \code{\link{BM_attributes}} \code{\link[biomaRt:listDatasets]{biomaRt::listDatasets}} \code{\link[biomaRt:listAttributes]{biomaRt::listAttributes}}
#' @examples
#' ## Get an overview of all valid species in BioMart
#' tab <- listSpecies()
#'
#' ## Check if Mus musculus is coverd by BioMart
#' "mmusculus" %in% listSpecies()$species
#'
#' @references
#' \insertAllCited{}
#'
#' @export
listSpecies <- function(){
  if(!curl::has_internet()){
    stop("Could not connect to BioMart. Check your internet connection")
  }
  tab <- biomaRt::listDatasets(biomaRt::useMart("ensembl"))
  tab[,1] <- gsub("_gene_ensembl","",tab[,1]); colnames(tab)[1]="species"
  return(tab)
}

#' @title List valid BioMart attributes
#'
#' @description
#' \code{BM_attributes} allows users to overview all valid BioMart attributes
#' for a given species using the \code{\link[biomaRt:listAttributes]{biomaRt::listAttributes}} function
#' from the \code{biomaRt} package \insertCite{biomart1,biomart2}{GEDI}.
#' By searching through the available BioMart attributes, the user can find the
#' BioMart attribute corresponding to the type of reporter IDs used in a dataset.
#' @name BM_attributes
#' @usage
#' BM_attributes(species)
#' @param species A species covered by BioMart Ensembl. The \code{species} must
#' be written with the first letter of the genus name followed by the specific
#' name, e.g. \emph{Homo sapiens} is \dQuote{hsapiens}, and \emph{Mus musculus} is
#' \dQuote{mmusculus}. Use \code{\link{listSpecies}} to find valid species inputs.
#'
#' @return A table with names and descriptions for all BioMart attributes for a
#' given \code{species}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{GEDI}} \code{\link{listSpecies}} \code{\link[biomaRt:listDatasets]{biomaRt::listDatasets}} \code{\link[biomaRt]{listAttributes}}
#' @examples
#' ## For species: Bos taurus (cow)
#' # Find all BioMart Attributes for Bos taurus
#' # if Bos taurus is a valid species in BioMart
#' species <- "btaurus"
#' if (species %in% listSpecies()){
#'     BM_attributes(species = species)
#' }


#' @export
BM_attributes <- function(species){
  if(!curl::has_internet()){
    stop("Could not connect to BioMart. Check your internet connection")
  } else if(!species %in% listSpecies()$species){
    stop(gettextf("Invalid argument for %s: %s. See all species covered by BioMart and ensembl using the listSpecies() function.",
                  sQuote("species"),dQuote(species)))
  }
  mart <- biomaRt::useDataset(paste(species,"_gene_ensembl",sep=""),biomaRt::useMart("ensembl"))
  tab <- biomaRt::listAttributes(mart)
  colnames(tab)[1] = "name"
  return(tab)
}
