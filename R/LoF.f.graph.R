#'Create frequency spectra
#'
#'This function generates a site-frequnecy spectra based on a dataframe of data created
#'from the lof_funct function.
#'
#'@param data df generated from either lof_funct.R or lof_pos.R
#'@param vcf vcf file from which lof_funct.R was ran on
#'@param outfile file path for the output of the graph
#'@return null
#'@export

## -----------------------------------------------------------------------------------------------------------------------------------
funct_graph<-function(data, vcf, outfile){
  #number of columns (individuals)
  maxc<-ncol(vcf)
  #Converting df to double
  L<-(as.numeric(unlist(data[2])))
  print(L)
  print(maxc)
  print(typeof(L))
  freq <- L / (maxc - 8)
  freq = freq[! freq %in% c(0)]
  #bins for data set = b
  b = sqrt(length(freq))
  #rounding up to nearest integer
  b = ceiling(b)
  allele.df<- data.frame(freq)
  pdf(file=outfile, width=6.5, height = 6.5)
  print(
    funct.graph<-ggplot(allele.df, aes(x=freq)) +
      geom_histogram(bins=b, binwidth = (max(freq)/b), fill="#69b3a2", color="#e9ecef", alpha=0.9, aes(y=stat(count/sum(count)))) +
      ggtitle("Functional-Class LoF Frequency") +
      theme(panel.background=element_rect(fill='transparent', colour='black'), panel.grid = element_line(colour='grey75'), axis.text=element_text(colour='black'), axis.ticks=element_line(colour='black')) + scale_y_continuous(name="Proportion") + xlab("Allele Frequency")
  )
  dev.off()
  return(funct.graph)
}

## -----------------------------------------------------------------------------------------------------------------------------------
pos_graph<-function(data, vcf, outfile){
  maxc<-ncol(vcf)
  pos.lof.list <- (as.numeric(unlist(data[2])))
  freq <- pos.lof.list / (maxc - 8)
  freq = freq[! freq %in% c(0)]
  #bins for data set = b
  b = sqrt(length(freq))
  #rounding up to nearest integer
  b = ceiling(b)
  allele.df<- data.frame(freq)
  pdf(file=outfile, width=6.5, height = 6.5)
  print(
    pos.graph<-ggplot(allele.df, aes(x=freq)) +
      geom_histogram(bins=b, binwidth = (max(freq)/b), fill="#69b3a2", color="#e9ecef", alpha=0.9, aes(y=stat(count/sum(count)))) +
      ggtitle("Position Based LoF Frequency") +
      theme(panel.background=element_rect(fill='transparent', colour='black'), panel.grid = element_line(colour='grey75'), axis.text=element_text(colour='black'), axis.ticks=element_line(colour='black')) + scale_y_continuous(name="Proportion") + xlab("Allele Frequency")
  )
  dev.off()
  return(pos_graph)
}
