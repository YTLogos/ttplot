#' search.names
#'
#' @keywords internal
#'
#' @return Nothing; internal function
#'
search.names <- function(term,dfnames){
  for(i in 1:length(term)){
    res <- grep(paste0("\\b",term[i],"\\b"),dfnames, ignore.case = TRUE)
    if(length(res)>0){
      if(length(res)==1){
        return(dfnames[res])
      }
    }
  }
}
