# global variables to escape r cmd check
utils::globalVariables(c("index","marker","chrom_alt","xbreaks"))
#' Make manhattan plot with full ggplot customizability
#'
#' This function is provided to make manhattan plot with full ggplot customizability. So next
#' we can customize the manhattan plot with kinds of functions of ggplot2 and add additional layers.
#' @param gwasres a data frame of gwas results.
#' @param snp Name of the column containing SNP identifiers; default is NA.
#' @param bp Name of the column containing the SNP positions; default is NA.
#' @param chrom Name of the column containing the chromosome identifers; default is NA.
#' @param pvalue Name of the column containing the p values; default is NA.
#' @param vlinetype the type of vline (geom_vline()). The default is "solid".
#' @param vlinesize the size of the vline. The default is 0.75.
#' @param title the title of manhattan plot. The default is "Manhattan Plot".
#' @param color the colors of alternate chromosome. The default is "#FF8C00" and "#556B2F"
#' @param pointsize the size of point. The default is 1.25.
#' @param file.output a logical, if file.output=TRUE, the result will be saved.
#' if file.output=FALSE, the result will be printed. The default is TRUE.
#' @param file a character, users can choose the different output formats of plot,
#' so far, "jpeg", "pdf", "png", "tiff" can be selected by users. The default is "png".
#' @param dpi a number, the picture element for .jpeg, .png and .tiff files. The default is 300.
#' @param output a character, the name of your trait. The default is "Trait"
#'
#' @author Tao Yan <\email{tyan@zju.edu.cn}> |
#' <\href{https://taoyan.netlify.com/}{https://taoyan.netlify.com/}>
#'
#' @import ggplot2
#' @import gtools
#' @return a manhattan plot based on ggplot2.
#'
#' @export
#'
#' @examples
#' ggmanhattan(gwas_test)

ggmanhattan <- function(
              gwasres,
              snp=NA,
              bp=NA,
              chrom=NA,
              pvalue=NA,
              index=NA,
              file.output=FALSE,
              file="png",
              output="Trait",
              dpi=300,
              vlinetype="solid",
              vlinesize=0.75,
              title="Manhattan Plot",
              color= c("#FF8C00", "#556B2F"),
              pointsize= 1.25,
              verbose=TRUE,...)
  {
  dfnames <- names(gwasres)
  if(is.na(chrom)){
    chrom <- search.names(c("chr","chrom","chromosome"), dfnames)
    if(is.null(chrom)){
      stop("Couldn't find the chromosome column. Please specify the name of the column with chromosome ids(chr,chrom or chromosome)")
    }
  }
  if(is.na(snp)){
    snp <- search.names(c("snp","snpid","rs","rsid"), dfnames)
    if(is.null(snp)){
      stop("Couldn't find the snp column. Please specify the name of the column with snp ids(snp, snpid, rs or rsid)")
    }
  }
  if(is.na(bp)){
    bp <- search.names(c("bp","pos","position"), dfnames)
    if(is.null(bp)){
      stop("Couldn't find the bp column. Please specify the name of the column with bp ids(bp, pos or position)")
    }
  }
  if(is.na(pvalue)){
    pvalue <- search.names(c("p", "p-value", "pvalue", "pval"), dfnames)
    if(is.null(chrom)){
      stop("Couldn't find the pvalue column. Please specify the name of the column with pvalue ids(p, pvalue, pval or p-value)")
    }
  }
  df <- as.data.frame(gwasres)
  df$chrom <- df[ ,chrom]
  df$chrom <- as.character(df$chrom)
  df$bp <- as.numeric(as.character(df[ ,bp]))
  df$pvalue <- as.numeric(as.character(df[ ,pvalue]))
  df$snp <- df[ ,snp]

  suppressMessages(require(gtools, quietly = TRUE))
  if(is.na(index)){
    df <- df[order(df$bp), ]
    df <- df[mixedorder(df$chrom), ]
    df$index <- 1:nrow(df)
  }
  else {
    df$index <- df[ ,index]
    df <- df[order(df$index), ]
    df <- df[mixedorder(df$chrom), ]
  }

  #calculate the numbers of chromosome
  chrnum <- data.frame(table(df$chrom))
  chrnum$Var1 <- as.character(chrnum$Var1)
  chrnum <- chrnum[mixedorder(chrnum$Var1), ]

  #marker the odd and even chromosome
  chrom_odd <- as.character(chrnum$Var1[seq(1,nrow(chrnum),2)])
  df$chrom_alt <- replace(df$chrom, df$chrom %in% chrom_odd, 0)
  df$chrom_alt <- replace(df$chrom_alt, df$chrom_alt !=0, 1)

  #-Log10 transform the pvalue

  df$marker <- -log10(df$pvalue)

  #specify the y limit
  ymax <- max(df$marker)+1.5

  #specify x axis tick points

  df_split <- split(df, df$chrom)
  xbreaks <- sapply(df_split, function(x){
    midpoint <- length(x$index)/2
    if(midpoint < 1) midpoint <- 1
    return(x$index[midpoint])
  })

  #calculate the number of SNPs

  snp_num <- cumsum(chrnum$Freq)

  #make the manhattan plot

  p1 <- ggplot(df, aes(x=index, y=marker, colour=as.factor(chrom_alt)))+
    geom_point(size=pointsize)+
    scale_x_continuous(breaks = xbreaks, labels = names(xbreaks), expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), limits = c(0, ymax))+
    guides(colour=FALSE)+
    labs(x="Chromosomes",y=expression(bold(-log[10]~(pvalue))), title = "Manhattan Plot")+
    geom_vline(xintercept = snp_num, colour="#C0C0C0", size=vlinesize, linetype=vlinetype)+
    scale_color_manual(values = color)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size=1),
          axis.title = element_text(face = "bold",size = 20),
          axis.text = element_text(face = "bold",size=16, colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  class(p1) <- append(class(p1),"ggman")

  #whether to save the plot or not

  if(!file.output){
    if(verbose) cat("The Manhattan Plot is Plotting! Please wait for a moment...", sep = "\n")
    return(p1)
  }
  if (file.output){
    if(verbose) cat("The Manhattan Plot is Plotting! Please wait for a moment...", sep = "\n")
    if(file=="jpg")	jpeg(paste("The_Manhattan_Plot_of_", output, ".jpg", sep=""), width = 15*dpi, height=7*dpi, res=dpi, quality = 100)
    if(file=="pdf")	pdf(paste("The_Manhattan_Plot_of_", output, ".pdf", sep=""), width = 15, height=7)
    if(file=="tiff")	tiff(paste("The_Manhattan_Plot_of_", output, ".tiff", sep=""), width = 15*dpi, height=7*dpi, res=dpi)
    if(file=="png")	png(paste("The_Manhattan_Plot_of_", output, ".png", sep=""), width = 15*dpi, height=7*dpi, res=dpi)
    par(xpd=TRUE)
  }else{
    if(is.null(dev.list())) dev.new(width =9, height=7)
    par(xpd=TRUE)
  }
  print(p1)
  if(file.output) dev.off()
  if(file.output & verbose) cat(paste("The Manhattan Plot is stored in: ", getwd(), sep = ""), sep = "\n")
  }














