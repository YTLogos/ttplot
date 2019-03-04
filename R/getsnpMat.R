#' Get the snpMatrix from the vcf file (plink format)
#'
#' This function defined to obtain the snpMatrix(with snp ID) from the vcf file
#' @param vcffile The vcf file (plink format).
#' @return The snpMAtrix for the visulization of LDheatmap.
#' @export
#' @examples
#' getsnpMat(vcffile)

getsnpMat <- function(vcffile){
  if(!require(vcfR)) BiocManager::install("vcfR")
  if(!require(genetics)) BiocManager::install("genetics")
  suppressMessages(require(vcfR))
  suppressMessages(require(genetics))
  snp_data <- read.vcfR(vcffile, verbose = TRUE)
  snp_gt <- snp_data@gt
  snp_gt <- as.data.frame(snp_gt[,-1])
  info <- as.data.frame(snp_data@fix)
  info$ID <- paste(info$CHROM,info$POS,sep = "_")
  snp_geno <- cbind(info[,c(4,5)],snp_gt)
  snp_geno[,1:ncol(snp_geno)] <- lapply(snp_geno[,1:ncol(snp_geno)], as.character)
  for (i in 1:nrow(snp_geno)){
    for (j in 3:ncol(snp_geno)){
      if (is.na(snp_geno[i,j])){
        snp_geno[i,j] <- NA
      }else if (snp_geno[i,j]=="0/0"){
        snp_geno[i,j] <- paste(snp_geno[i,1],snp_geno[i,1], sep = "/")
      }else if (snp_geno[i,j]=="0/1"){
        snp_geno[i,j] <- paste(snp_geno[i,1],snp_geno[i,2], sep = "/")
      }else if (snp_geno[i,j]=="1/1"){
        snp_geno[i,j] <- paste(snp_geno[i,2],snp_geno[i,2], sep = "/")
      }
    }
  }
  snpMat <- t(snp_geno[,-c(1,2)])
  snpNames <- info$ID
  colnames(snpMat) <- snpNames
  rownames(snpMat) <- NULL
  snpMat <- as.data.frame(snpMat)
  for (i in 1:ncol(snpMat)){
    snpMat[, i] <- genetics::as.genotype(snpMat[, i])
  }
  return(snpMat)
}
