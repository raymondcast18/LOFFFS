% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lof_pos.R
\name{lof_pos}
\alias{lof_pos}
\title{Count loss of function mutations based on position in the genome.}
\usage{
lof_pos(x, y)
}
\arguments{
\item{x}{formatted vcf file}

\item{y}{File path for .csv file}
}
\value{
data frame with gene name and the LoF mutations within the gene
}
\description{
This function counts LoF mutations based on postion.
This function outputs a df and csv file with the data. This assumes the function load_vcf was run on
a vcf file. If you use the package to generate a SFS, you need assign the output of the function to a variable
in your environment
}
