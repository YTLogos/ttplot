#' Get the LDhaetmap from the vcf file (plink format)
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file.
#' @return the LDheatmap.
#' @export
#' @examples
#' MyLDheatMap(vcffile)

MyLDheatMap <- function(vcffile){
  name <- basename(vcffile)
  title <- sub(".vcf","",name)
  gdat_snp <- ttplot::snpMat(vcffile)
  info <- ttplot::snpInfo(vcffile)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  myLDheatmap <- LDheatmap(gdat_snp, info$POS,
                           color = rgb.palette(100),
                           flip = TRUE,
                           title = paste0("The LDheatmap of ",title))
#  print(myLDheatmap)
}
