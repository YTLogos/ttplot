#' Get the LDhaetmap from the vcf file (plink format)
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file.
#' @return the LDheatmap.
#' @export
#' @examples
#' MyLDheatMap(vcffile)

MyLDheatMap <- function(vcffile,title=""){
  if(!require(LDheatmap)) BiocManager::install("LDheatmap")
  require(LDheatmap)
  gdat_snp <- ttplot::getsnpMat(vcffile)
  info <- ttplot::getsnpInfo(vcffile)
  snp_dist <- as.numeric(info$POS)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  MyHeatmap <- LDheatmap(gdat_snp,
                         genetic.distances = info$POS,
                         color = grey.colors(20))
  myLDheatmap <- LDheatmap(MyHeatmap,
                           color = rgb.palette(100),
                           flip = TRUE,
                           title = paste0("The LDheatmap of ",title))
}
