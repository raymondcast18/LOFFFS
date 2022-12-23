#'Count loss of function mutations based on functional status.
#'
#'This function counts LoF mutations based on functional status and
#'outputs a df and csv file with the data. This assumes the function load_vcf was run on
#'a vcf file and if you use the package to generate a SFS, you assign the output of the function to a variable
#'in your environment
#'
#'@param x formatted vcf file
#'@param y File path for .csv file
#'@return data frame with gene name and the LoF mutations within the gene
#'@export
lof_funct<-function(x, y){
  ## -----------------------------------------------------------------------------------------------------------------------------------
  #Isolates the gene id in the vcf file
  maxc<-ncol(x)
  maxr<-nrow(x)
  i<-1
  j<-8
  while(i<=maxr){
    gene<-toString(x[i, j], width=NULL)
    k<-strsplit(gene, '|', fixed = TRUE)
    x[i, j] = k[[1]][5]
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
  #pos.lof.list<-c()
  i<-1
  hom<-"1\\|1"
  het<-"0\\|1"
  lof.count<-0
  #pos.count<-0
  gene2<-0
  m<-1
  pos2<-0
  #This is the algorithm that iterates through the rows checking for variables lof1 and lof2 in the data and keeping track of LoF mutations per
  #individual will check if lof mut is on same gene for same individual
  gene.list = c()
  #pos.list<-c()
  individual = c(rep(0, maxc))
  for(i in 1:maxr){
    j<-10
    gene1<-x[i, 8]
    pos1<-x[i, 2]
    #grep checks if we are on a different gene or on the current one with a binary value
    genecheck<-grepl(gene1, gene2)
    if(genecheck==FALSE){
      gene.list=append(gene.list, gene1)
    }
    if(genecheck==FALSE & i!=1){
      L = append(L, lof.count)
      lof.count=0
      individual = c(rep(0, maxc))
      m=m+1
    }
    while(j<=maxc){
      value<-x[i, j]
      value_string<-toString(value, width = NULL)
      comp.hom<-grepl(hom, value_string)
      comp.het<-grepl(het, value_string)

      if(individual[j]!=1 && comp.hom==TRUE||comp.het==TRUE){
        lof.count = lof.count +1
        individual[j] = 1
      }
      #if(comp.hom==TRUEcomp.het==TRUE){
      # pos.count = pos.count + 1
      #}
      j=j+1
    }
    gene2<-x[i, 8]
    # pos2 = x[i+1, 2]
    #pos.check = grepl(pos1, pos2)
    #if(pos.check==FALSE){
    # pos.lof.list=append(pos.lof.list, pos.count)
    #pos.list=append(pos.list, pos1)
    # pos.count=0
    # pos1=pos2
    # }
  }
  #L = a list of final sums of LoF mutations
  L = append(L, lof.count)


  ##-------------------------------------------------------------------------------------------
  #creation of csv files

  LoF.Funct.df<-data.frame(gene.list, L)
  colnames(LoF.Funct.df)[2]="LoF Count"
  write.csv(LoF.Funct.df, file=y, row.names=FALSE)
  print(LoF.Funct.df)
  # LoF.pos.df<-data.frame(pos.list, pos.lof.list)
  # colnames(LoF.pos.df)[1]="Position"
  # colnames(LoF.pos.df)[2]="LoF Count"
  # write.csv(LoF.pos.df, file="LoF.pos.csv", row.names=FALSE)
  return(LoF.Funct.df)
}
