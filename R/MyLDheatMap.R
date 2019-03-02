#' Get the LDhaetmap from the vcf file (plink format)
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file.
#' @return the LDheatmap.
#' @export
#' @examples
#' MyLDheatMap(vcffile)

MyLDheatMap <- function(vcffile){
  title <- sub(".R","",basename(vcffile))
  gdat_snp <- snpMat(vcffile)
  info <- snpInfo(vcffile)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  myLDheatmap <- LDheatmap(gdat_snp, info$POS, color = rgb.palette(100), flip = TRUE, title = title)
  print(myLDheatmap)
}
