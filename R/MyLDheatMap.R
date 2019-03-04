#' Get the LDhaetmap from the vcf file (plink format) directly
#'
#' This function defined to obtain the LDheatmap from the vcf file directly.
#' @param vcffile The plink format vcf file. More detail can see View(test_vcf)
#' @param file.output a logical, if file.output=TRUE, the result will be saved.
#' if file.output=FALSE, the result will be printed. The default is TRUE
#' @param file a character, users can choose the different output formats of plot, so far, "jpeg", "pdf", "png", "tiff" can be selected by users. The default is "png".
#' @param title a character, the title of the LDheatmap will be "The LDheatmap of title".
#' the default is "region:". I suggest users use your own title.
#' @param dpi a number, the picture element for .jpeg, .png and .tiff files. The default is 300.
#' @param verbose whether to print the reminder.
#'
#' @author Tao Yan <\email{tyan@zju.edu.cn}> |
#' <\href{https://taoyan.netlify.com/}{https://taoyan.netlify.com/}>
#'
#' @return the LDheatmap.
#'
#' @export MyLDheatMap
#'
#' @examples
#' MyLDheatMap(system.file("extdata","test.vcf", package = "ttplot"), title="your title")


MyLDheatMap <- function(
  vcffile,
  file.output=TRUE,
  file="png",
  output="region",
  title="region:",
  verbose=TRUE,
  dpi=300
  )
  {
  if(!require(LDheatmap)) BiocManager::install("LDheatmap")
  suppressWarnings(suppressMessages(library(LDheatmap, quietly = T)))
  gdat_snp <- ttplot::getsnpMat(vcffile)
  info <- ttplot::getsnpInfo(vcffile)
  snp_dist <- as.numeric(info$POS)
  rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
  if (!file.output){
    if (verbose) print("The Ldheatmap is Plotting! Please wait for a moment...")
    LDheatmap(gdat_snp,
            genetic.distances = snp_dist,
            color = rgb.palette(100),
            flip = TRUE,title = paste0("The LDheatmap of ",title))
  }
  if (file.output){
    if(verbose) print("The Ldheatmap is Plotting! Please wait for a moment...")
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
  if(file.output) dev.off()
  if(file.output & verbose)	print(paste("Plot is stored in: ", getwd(), sep=""))
}
