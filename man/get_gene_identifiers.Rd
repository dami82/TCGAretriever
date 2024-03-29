% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcgaretriever_v_1_9e.R
\name{get_gene_identifiers}
\alias{get_gene_identifiers}
\title{Retrieve All Gene Identifiers}
\usage{
get_gene_identifiers(dryrun = FALSE)
}
\arguments{
\item{dryrun}{Logical. If TRUE, all other arguments (if any) are ignored and 
a representative example is returned as output. No Internet connection is 
required for executing the operation when `dryrun` is TRUE.}
}
\value{
Data Frame including gene identifiers.
}
\description{
Obtain all valid gene identifiers, including ENTREZ gene identifiers and 
HUGO gene symbols. Genes are classified according to the gene type 
(*e.g.*, 'protein-coding', 'pseudogene', 'miRNA', ...). Note that 
miRNA and phosphoprotein genes are associated with a negative entrezGeneId.
}
\examples{
# Set `dryrun = FALSE` (default option) in production!
x <- get_gene_identifiers(dryrun = TRUE)


}
\references{
\url{https://www.data-pulse.com/dev_site/TCGAretriever/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
