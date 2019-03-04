#' Obtain the snp info from the vcf file
#'
#' This function defined to get the snp info (ID, POS) from the vcf file.
#' @param file The vcf file name.
#' @return The snp info (ID, POS).
#' @export
#' @examples
#' snpInfo(sample.vcf)

getsnpInfo <- function(file){
  if(!require(vcfR)) BiocManager::install("vcfR")
  suppressMessages(require(vcfR))
  snp <- read.vcfR(file,verbose = FALSE)
  snp_info <- snp@fix
  snp_info <- as.data.frame(snp_info)
  snp_info$ID <- paste(snp_info$CHROM, snp_info$POS, sep = "_")
  snp_info <- snp_info%>%dplyr::select(ID, POS)
  snp_info$POS <- as.integer(as.character(snp_info$POS))
  return(snp_info)
}
