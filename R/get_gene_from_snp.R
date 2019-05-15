#' Get the candidate genes in regions based on significant SNPs of GWAS results.
#'
#' This function is developed to get the candidate genes in regions based on significant SNPs of GWAS results.
#' @param gff a data frame of all the gene (transcript), must have column names.
#' @param sig.snp a data frame of significant SNPs.
#' @param distance numeric (bp), it is to define the region. The default is 50000, you need to choose it based on the LD distance in your study.
#' @param geneid Name of the column containing the geneid in gff file; default is NA.
#'
#' @param gff.chrom Name of the column containing the chromosome identifers in the gff file; default is NA.
#' @param snp.chrom Name of the column containing the chromosome identifers in the snp.sig file; default is NA.
#' @param pvalue Name of the column containing the p values in snp.sig file; default is NA.
#' @param snp_location Name of the column containing the snp position in snp.sig file; default is NA.
#' @param file.save a logical, if file.output=TRUE, the result will be saved.
#' if file.output=FALSE, the result will be printed. The default is TRUE.
#' @param file.type a character, users can choose the different output formats,
#' so far, "csv", "txt", "xlsx" can be selected by users. The default is "csv".
#'
#' @author Tao Yan <\email{tyan@zju.edu.cn}> |
#' <\href{https://taoyan.netlify.com/}{https://taoyan.netlify.com/}>
#'
#' @import writexl
#' @return a data.frame contain the candidate genes with start,end,genid etc
#'
#' @export
#'
#' @examples
#' get_gene_from_snp(gff,sig.snp)

get_gene_from_snp <- function(
                    gff,
                    sig.snp,
                    distance=75000,
                    file.save=TRUE,
                    file.type="csv",
                    gff.chrom=NA,
                    snp.chrom=NA,
                    geneid=NA,
                    pvalue=NA,
                    gene_start=NA,
                    gene_end=NA,
                    snp_location=NA,
                    verbose=TRUE,...)
  {
  gff_names <- names(gff)
  snp_names <- names(sig.snp)
  if(is.na(gff.chrom)){
    gff.chrom <- search.names(c("chr", "chrom", "chromosome"), gff_names)
    if(is.null(gff.chrom)){
      stop("Couldn't find the chromosome column. Please specify the name of the column with chromosome ids(chr, chrom or chromosome)")
    }
  }

  if(is.na(gene_start)){
    gene_start <- search.names(c("start", "gene_start", "genestart", "begin"), gff_names)
    if(is.null(gene_start)){
      stop("Couldn't find the gene_start column. Please specify the name of the column with gene_start ids(start, gene_start, genestart or begin)")
    }
  }
  if(is.na(gene_end)){
    gene_end <- search.names(c("end", "gene_end", "geneend"), gff_names)
    if(is.null(gene_end)){
      stop("Couldn't find the gene_end column. Please specify the name of the column with gene_end ids(end, gene_end or geneend)")
    }
  }
  if(is.na(geneid)){
    geneid <- search.names(c("gene", "geneid", "gene_id", "gene_name","genename","transcript","transcript_id","transcriptid"), gff_names)
    if(is.null(geneid)){
      stop("Couldn't find the geneid column. Please specify the name of the column with geneid ids(gene, geneid, gene_id, gene_name, genename, transcript,transcript_id or transcriptid)")
    }
  }
  if(is.na(pvalue)){
    pvalue <- search.names(c("p", "p-value", "pval", "pvalue"), snp_names)
    if(is.null(pvalue)){
      stop("Couldn't find the pvalue column. Please specify the name of the column with pvalue ids(p, pvalue, pval or p-value)")
    }
  }
  if(is.na(snp_location)){
    snp_location <- search.names(c("bp","position","pos","snp_location","location"), snp_names)
    if(is.null(snp_location)){
      stop("Couldn't find the bp column. Please specify the name of the column with bp ids(bp, pos, snp-location, location or position)")
    }
  }
  if(is.na(snp.chrom)){
    snp.chrom <- search.names(c("chr", "chrom", "chromosome"), snp_names)
    if(is.null(snp.chrom)){
      stop("Couldn't find the chromosome column. Please specify the name of the column with chromosome ids(chr, chrom or chromosome)")
    }
  }
  gff$gff.chrom <- gff[ ,gff.chrom]
  gff$gene_start <- gff[ ,gene_start]
  gff$gene_end <- gff[ ,gene_end]
  gff$geneid <- gff[ ,geneid]
  sig.snp$pvalue <- sig.snp[ ,pvalue]
  sig.snp$snp_location <- sig.snp[ ,snp_location]
  sig.snp$snp.chrom <- sig.snp[ ,snp.chrom]

  r1 <- vector(mode = "numeric")
  r2 <- vector(mode = "numeric")

  for (i in 1:nrow(sig.snp)){
    if (sig.snp[i,][["snp_location"]]>=distance){
      r <- c(sig.snp[i,][["snp_location"]]-distance, sig.snp[i,][["snp_location"]]+distance)
    }else{
      r <- c(sig.snp[i,][["snp_location"]], sig.snp[i,][["snp_location"]]+distance)
    }
    r1 <- append(r1,r[1])
    r2 <- append(r2,r[2])
  }
  sig.snp$r1 <- r1
  sig.snp$r2 <- r2

  ## search significant genes
  findgene <- function(i){
    snp_gene <-  gff[gff$gene_start >= sig.snp$r1[i] & gff$gene_end <= sig.snp$r2[i] & gff$gff.chrom==sig.snp$snp.chrom[i], ]
    if (nrow(snp_gene)==0) {
      return(NULL)}else{
        snp_gene$snp_location <-  sig.snp$snp_location[i]
        snp_gene$p_value <- sig.snp$pvalue[i]
        return(snp_gene)
      }
  }

  genes_data <- future.apply::future_lapply(seq(1,nrow(sig.snp)),findgene)
  gene_snp <- do.call(rbind, genes_data)



  if(!file.save) {
    if(verbose) {
      cat(paste("The distance you choose is ", distance, "bp!", sep = ""), sep = "\n")
      cat(paste("You have ", nrow(sig.snp), " significant SNPs and ", nrow(gff), " genes!", sep = ""),sep = "\n")
      cat("Now we will extract the genes in the significant regions! This will need some time, please wait for severals minutes! ...", sep = "\n")
    }
  }
  if(file.save){
    if(verbose) {
      cat(paste("The distance you choose is ", distance, "bp!", sep = ""),sep = "\n")
      cat(paste("You have ", nrow(sig.snp), " significant SNPs and ", nrow(gff), " genes!", sep = ""),sep = "\n")
      cat("Now we will extract the genes in the significant regions! This will need some time, please wait! ...", sep = "\n")
    }
  }



  if(file.save){
    if(file.type=="csv"){
      write.csv(gene_snp,file = paste("genes_in_sig_regions","(",distance,"bp)",".csv"), row.names = FALSE, quote = FALSE)
    }
    if(file.type=="txt"){
      write.table(gene_snp,file = paste("genes_in_sig_regions","(",distance,"bp)",".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    if(file.type=="xlsx"){
      writexl::write_xlsx(gene_snp,path = paste("genes_in_sig_regions","(",distance,"bp)",".xlsx"))
    }
  }
  if(file.save & verbose) {
    cat("         \n")
    cat(paste("The output is stored in: ", getwd(), sep = ""))
  }
  if(!file.save & verbose) {
    cat("                                 ","The result is : ", "-------------------------------------------------------------",sep = "\n")
    return(gene_snp)
    }
}
