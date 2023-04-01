#'Count loss of function mutations based on functional status.
#'
#'This function counts LoF mutations based on functional status and
#'outputs a df and csv file with the data. This assumes the function load_vcf was run on
#'a vcf file and if you use the package to generate a SFS, you assign the output of the function to a variable
#'in your environment
#'
#'@param x formatted vcf file
#'@param y File path for .csv file
#'@param z optional file path for potential mutations that may not be LoF and need further investigation
#'@return data frame with gene name and the LoF mutations within the gene
#'@export
lof_funct<-function(x, y, z){
  ## -----------------------------------------------------------------------------------------------------------------------------------
  #Isolates the gene id in the vcf file
  maxc<-ncol(x)
  maxr<-nrow(x)
  i<-1
  j<-8
  #mutlimutant keeps track of mutations for the data set.
  multimutant<-c()
  maxc<-ncol(x)
  maxr<-nrow(x)
  i<-1
  j<-8
  bug1<-1
  while(i<=maxr){
    gene<-toString(x[bug1, j], width=NULL)
    k1<-strsplit(gene, ';', fixed = TRUE)
    k<-strsplit(gene, '|', fixed = TRUE)
    klength<-lengths(k)
    LoF.check<-grepl('LOF', k1[[1]][3])
    if(LoF.check==TRUE){
      x[bug1, j] = k[[1]][5]
      multimutant<-append(multimutant, k[[1]][2])
      bug1=bug1+1
    }
    else{
      x<-x[-bug1,]
    }
    i=i+1
  }


  ## ---------------------------------------------------------------------------------


  #summing LoF mutations in file
  #i=row in file
  #j=column in file

  #mdata and L is where final data will be store using datai and dataj to navigate the matrix
  #element<-1
  #m keeps track of gene number
  #lof1 and lof2 is what the code searches for in the vcf data file
  #sum used to check for lof1 and lof2 as a binary
  #gene2 checks if code is on a new gene
  #L is a list that will store our counts of LoF mutations
  L <- c()
  i<-1
  hom<-"1\\|1"
  het<-"0\\|1"
  lof.count<-0
  gene2<-0
  m<-1
  genenum<-1
  mutant.check<- c(0, 0)
  non.lof.list <- c()
  #This is the algorithm that iterates through the rows checking for variables lof1 and lof2 in the data and keeping track of LoF mutations per
  #individual will check if lof mut is on same gene for same individual
  #list of genes
  gene.list = c()
  individual = c(rep(0, maxc))
  for(i in 1:maxr){
    j<-10
    gene1<-x[i, 8]

    #grep checks if we are on a different gene or on the current one with a binary value
    genecheck<-grepl(gene1, gene2)
    if(genecheck==FALSE){
      gene.list=append(gene.list, gene1)
    }
    if(genecheck==FALSE & i!=1){
      mutant.check[1]=lof.count%%3
      if(mutant.check[1]==0 && mutant.check[2]!=1){
        non.lof.list<-append(non.lof.list, gene1)
      }
      L = append(L, lof.count)
      lof.count=0
      individual = c(rep(0, maxc))
      mutant.check<- c(0, 0)
      m=m+1
    }
    while(j<=maxc){
      value<-x[i, j]
      value_string<-toString(value, width = NULL)
      comp.hom<-grepl(hom, value_string)
      comp.het<-grepl(het, value_string)

      if(individual[j]!=1 && comp.hom==TRUE||comp.het==TRUE){
        lof.count = lof.count + 1
        individual[j] = 1
      }
      j=j+1
    }
    gene2<-x[i, 8]
    codon.check <- grepl(multimutant[i], "stop_gained")
    if(codon.check==TRUE){
      mutant.check[2] = 1
    }
  }
  #L = a list of final sums of LoF mutations
  L = append(L, lof.count)

  ##-------------------------------------------------------------------------------------------

  LoF.Funct.df<-data.frame(gene.list, L)
  colnames(LoF.Funct.df)[2]="LoF Count"
  write.csv(LoF.Funct.df, file=y, row.names=FALSE)
  if(missing(z)==FALSE){
    write.csv(non.lof.list, file = z, row.names = FALSE)
  }
  return(LoF.Funct.df)
}
