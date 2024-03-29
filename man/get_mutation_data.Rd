% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcgaretriever_v_1_9e.R
\name{get_mutation_data}
\alias{get_mutation_data}
\title{Retrieve Mutation Data corresponding to a Genetic Profile of Interest.}
\usage{
get_mutation_data(
  case_list_id,
  gprofile_id,
  glist = c("TP53", "E2F1"),
  dryrun = FALSE
)
}
\arguments{
\item{case_list_id}{String corresponding to the Identifier of a list of cases.}

\item{gprofile_id}{String corresponding to the Identifier of a genetic Profile of interest.}

\item{glist}{Vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMOL). ENTREZID
gene identifiers should be passed as numeric.}

\item{dryrun}{Logical. If TRUE, all other arguments (if any) are ignored and 
a representative example is returned as output. No Internet connection is 
required for executing the operation when `dryrun` is TRUE.}
}
\value{
data Frame inluding one row per mutation
}
\description{
Retrieve DNA Sequence Variations (Mutations) identified by exome sequencing projects.
This function is the workhorse of the TCGAretriever package for mutation data and can be used to fetch data 
concerning several genes at once. For retrieving non-mutation data, please use the `get_molecular_data()` function. 
For large queries (more than 500 genes), please use the `fetch_all_tcgadata()` 
function.
}
\examples{
# Set `dryrun = FALSE` (default option) in production!
x <- get_mutation_data(case_list_id = 'blca_tcga_sequenced',
                       gprofile_id = 'blca_tcga_mutations',
                       glist = c('TP53', 'PTEN'), dryrun = TRUE)
utils::head(x[, c(4, 7, 23, 15, 16, 17, 24, 18, 21)]) 
 

}
\references{
\url{https://www.data-pulse.com/dev_site/TCGAretriever/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
