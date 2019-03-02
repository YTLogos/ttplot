#' Convert the value of genotype into 0,1,2
#'
#' This function converts the value of genotype into 0,1,2 form 0/0,0/1,1/1.
#' @param snp_data The snp_data matrix.
#' @return The genodata with the value of 0,1,2.
#' @export
#' @examples
#' convertTonumeric(snp_data)

convertToNumeric <- function(snp_data){
  genodata <- matrix(NA, nrow = nrow(snp_data), ncol = ncol(snp_data))
  for (a in 1:nrow(snp_data)){
    for (b in 1:ncol(snp_data)){
      m <- as.numeric(unlist(strsplit(snp_data[a,b], "/"))[1])
      n <- as.numeric(unlist(strsplit(snp_data[a,b], "/"))[2])
      genodata[a,b] <- m+n
    }
  }
  rownames(genodata) <- rownames(snp_data)
  colnames(genodata) <- colnames(snp_data)
}
