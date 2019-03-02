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
  if(!require(genetics)) BiocManager::install("genetics")
  require(vcfR)
  require(snpStats)
  snp <- read.vcfR(vcffile, verbose = FALSE)
  snp_gt <- snp@gt
  snp_gt <- as.data.frame(snp_gt[,-1])
  info <- as.data.frame(snp@fix)
  info$ID <- paste(info$CHROM,info$POS,sep = "_")
  snp_data <- cbind(info[,c(4,5)],snp_gt)
  snp_data[,1:ncol(snp_data)] <- lapply(snp_data[,1:ncol(snp_data)], as.character)
  for (i in 1:nrow(snp_data)){
    for (j in 3:ncol(snp_data)){
      if (is.na(snp_data[i,j])){
        snp_data[i,j] <- NA
      }else if (snp_data[i,j]=="0/0"){
        snp_data[i,j] <- paste(snp_data[i,1],snp_data[i,1], sep = "/")
      }else if (snp_data[i,j]=="0/1"){
        snp_data[i,j] <- paste(snp_data[i,1],snp_data[i,2], sep = "/")
      }else if (snp_data[i,j]=="1/1"){
        snp_data[i,j] <- paste(snp_data[i,2],snp_data[i,2], sep = "/")
      }
    }
  }
  snpMat <- t(snp_data)
  info <- snpInfo(vcffile)
  snpNames <- info$ID
  colnames(snpMat) <- snpNames
  for (i in 1:ncol(snpMat)){
    snpMat[, i] <- genetics::as.genotype(snpMat[, i])
  }
  return(snpMat)
}
