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
  suppressWarnings(suppressMessages(library(vcfR,quietly = T)))
  suppressWarnings(suppressMessages(library(genetics, quietly = T)))
  snp_data <- read.vcfR(vcffile, verbose = FALSE)
  snp_gt <- snp_data@gt
  snp_gt <- as.data.frame(snp_gt[,-1])
  info <- as.data.frame(snp_data@fix)
  info$ID <- paste(info$CHROM,info$POS,sep = "_")
  snp_geno <- cbind(info[,c(4,5)],snp_gt)
  snp_geno[,1:ncol(snp_geno)] <- lapply(snp_geno[,1:ncol(snp_geno)], as.character)
  for (i in 1:nrow(snp_geno)){
    snp_geno[i,] <- ifelse(snp_geno[i,]=="0/0",paste(snp_geno[i,1],snp_geno[i,1], sep = "/"),ifelse(snp_geno[i,]=="0/1",paste(snp_geno[i,1],snp_geno[i,2], sep = "/"),ifelse(snp_geno[i,]=="1/1",paste(snp_geno[i,2],snp_geno[i,2], sep = "/"),snp_geno[i,])))
  }
  snpMat <- t(snp_geno[,-c(1,2)])
  snpNames <- info$ID
  colnames(snpMat) <- snpNames
  rownames(snpMat) <- NULL
  snpMat <- as.data.frame(snpMat)
  snpMat[,1:ncol(snpMat)] <- lapply(snpMat[,1:ncol(snpMat)], genetics::as.genotype)
  return(snpMat)
}
