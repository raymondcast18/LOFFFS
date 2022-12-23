#' @title Load Vcf
#' This functions loads a vcf file and formats it to be readable in R.
#' @param infile path to infile
#' @return formatted vcf file
#' @export

load_vcf<-function(infile){
  tmp_vcf<-readLines(infile)
  tmp_vcf_data<-read.table(infile, stringsAsFactors = FALSE)
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  names(tmp_vcf_data)<-vcf_names
  save(tmp_vcf_data, file = "tmp_vcf_data")
  return(tmp_vcf_data)
}
