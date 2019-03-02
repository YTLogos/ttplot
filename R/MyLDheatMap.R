#' Get the LDhaetmap from the vcf file (plink format)
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file.
#' @return the LDheatmap.
#' @export
#' @examples
#' MyLDheatMap(vcffile)

MyLDheatMap <- function(vcffile){
  if(!require(LDheatmap)) BiocManager::install("LDheatmap")
  require(LDheatmap)
  name <- basename(vcffile)
  title <- sub(".vcf","",name)
  gdat_snp <- ttplot::snpMat(vcffile)
  info <- ttplot::snpInfo(vcffile)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  MyHeatmap <- LDheatmap(gdat_snp,
                         genetic.distances = info$POS,
                         color = grey.colors(20))
  myLDheatmap <- LDheatmap(MyHeatmap,
                           color = rgb.palette(100),
                           flip = TRUE,
                           title = paste0("The LDheatmap of ",title))
}
