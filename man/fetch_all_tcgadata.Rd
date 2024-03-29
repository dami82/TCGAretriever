% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcgaretriever_v_1_9e.R
\name{fetch_all_tcgadata}
\alias{fetch_all_tcgadata}
\title{Fetch All Molecular Data for a Cancer Profile of Interest.}
\usage{
fetch_all_tcgadata(case_list_id, gprofile_id, mutations = FALSE)
}
\arguments{
\item{case_list_id}{string corresponding to the identifier of the TCGA Case List of interest}

\item{gprofile_id}{string corresponding to the identifier of the TCGA Profile of interest}

\item{mutations}{logical. If TRUE, extended mutation data are fetched instead of the standard TCGA data}
}
\value{
A data.frame is returned, including the desired TCGA data. 
Typically, rows are genes and columns are cases. 
If "extended mutation" data are retrieved (mutations = TRUE), 
rows correspond to individual mutations 
while columns are populated with mutation features
}
\description{
Recursively query cbioportal to retrieve data corresponding to all 
available genes. Data are returned as a `data.frame` 
that can be easily manipulated for downstream analyses.
}
\examples{
# The examples below require an active Internet connection.
# Note: execution may take several minutes.
\dontrun{
# Download all brca_pub mutation data (complete samples)
all_brca_MUT <- fetch_all_tcgadata(case_list_id = "brca_tcga_pub_complete", 
                                   gprofile_id = "brca_tcga_pub_mutations", 
                                   mutations = TRUE)

# Download all brca_pub RNA expression data (complete samples)
all_brca_RNA <- fetch_all_tcgadata(case_list_id = "brca_tcga_pub_complete", 
                                   gprofile_id = "brca_tcga_pub_mrna", 
                                   mutations = FALSE)
}



}
\references{
\url{https://www.data-pulse.com/dev_site/TCGAretriever/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
