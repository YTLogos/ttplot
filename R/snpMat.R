#' Get the snpMatrix from the vcf file (plink format)
#'
#' This function defined to obtain the snpMatrix(with snp ID) from the vcf file
#' @param vcffile The vcf file (plink format).
#' @return The snpMAtrix for the visulization of LDheatmap.
#' @export
#' @examples
#' snpMat(vcffile)

snpMat <- function(vcffile){
  if(!require(vcfR)) BiocManager::install("vcfR")
  if(!require(snpStats)) BiocManager::install("snpStats")
  require(vcfR)
  require(snpStats)
  snp <- read.vcfR(vcffile, verbose = FALSE)
  snp_gt <- snp@gt
  snp_gt <- snp_gt[,-1]
  snpMat <- t(snp_gt)
  gdat_snp <- convertToNumeric(snpMat)
  info <- snpInfo(vcffile)
  snpNames <- info$ID
  colnames(gdat_snp) <- snpNames
  gdat_snp <- as(gdat_snp, "SnpMatrix")
}
