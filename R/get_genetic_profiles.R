get_genetic_profiles <-
function(csid = NULL){
  if(!is.null(csid)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getGeneticProfiles"
    my_url <- paste(my_url, "&cancer_study_id=", tolower(as.character(csid)), sep = "")
    return(basic_tcga_query(my_url))
  } else {
    message("Missing cancer study ID")
  }
}
