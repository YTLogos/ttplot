#' Get the LDhaetmap from the vcf file (plink format)
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file.
#' @return the LDheatmap.
#' @export
#' @examples
#' MyLDheatMap(test_vcf, title="your title")

MyLDheatMap <- function(
  vcffile,
  file.output=TRUE,
  file="png",
  output="region",
  title="",
  verbose=TRUE,
  dpi=300
  )
  {
  if(!require(LDheatmap)) BiocManager::install("LDheatmap")
  suppressMessages(require(LDheatmap))
  gdat_snp <- ttplot::getsnpMat(vcffile)
  info <- ttplot::getsnpInfo(vcffile)
  snp_dist <- as.numeric(info$POS)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  if (!file.output){
    if (verbose) print("The Ldheatmap Plotting...")
    LDheatmap(gdat_snp,
            genetic.distances = snp_dist,
            color = rgb.palette(100),
            flip = TRUE,title = paste0("The LDheatmap of ",title))
  }
  if (file.output){
    if (verbose) print("The Ldheatmap Plotting...")
    if(file=="jpg")	jpeg(paste("LDheatmap of ", output, ".jpg", sep=""), width = 9*dpi, height=7*dpi, res=dpi, quality = 100)
    if(file=="pdf")	pdf(paste("LDheatmap of ", output, ".pdf", sep=""), width = 9, height=7)
    if(file=="tiff")	tiff(paste("LDheatmap of ", output, ".tiff", sep=""), width = 9*dpi, height=7*dpi, res=dpi)
    if(file=="png")	png(paste("LDheatmap of ", output, ".png", sep=""), width = 9*dpi, height=7*dpi, res=dpi)
    par(xpd=TRUE)
  }else{
    if(is.null(dev.list())) dev.new(width =9, height=7)
    par(xpd=TRUE)
  }
  LDheatmap(gdat_snp,
            genetic.distances = snp_dist,
            color = rgb.palette(100),
            flip = TRUE,
            title = paste0("The LDheatmap of ",title))
  if (file.output) dev.off()
  if(file.output & verbose)	print(paste("Plots are stored in: ", getwd(), sep=""))
}
