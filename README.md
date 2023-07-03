# LoFFFS
Loss-of-Function Functional Frequency Spectra  R Package

This is an R package that can generate site grequency spectra of LoF mutations based on the the functional status of the gene they lie on
and/or the position of the LoF mutation. This package suggests the installation of ggplot2 and scales for the graph generation but is not required
if you require only the raw data. Additionally, this package assumes that the data taken in is a vcf file in the format of the 1001 Arabidopsis 
TAIR 10 Genome.

# DEPENDENCIES
The package dplyr is required for use of this package.
The packages scales and ggplot2 are recommended for SFS generation however, it is not required.

# INSTALLATION
The packages devtools or githubinstall can be used for installation.
To install the following commands can be used:

Devtools Installation

library(devtools)

install_github("raymondcast18/LOFFFS")

Githubinstall Installation

library(githubinstall)

githubinstall("LOFFFS")

# USAGE
To use this package, a vcf file annotated with the variant annotation software like SnpEff is required to be in the same format
as the 1001 Arabidopsis TAIR 10 Genome.

To begin load and assign the vcf to a variable the R enviroment as shown below:

vcf<-load_vcf("filepath")

Using the assigned variable, you can then use any of the core subfunctions of LoFFFS.

To obtain a functional frequency spectra of a dataset:

lof_funct(vcf, outfilepath, optionalfilepath).

To obtain a position-based spectra of a dataset:

pos_funct(vcf, outfilepath)

If you have a particular geneset that are of interest you can add a filepath with a textfile containing a list of genes of interest:

lof_genes_funct(vcf, outfilepath, genetextfilepath)

lof_genes_pos(vcf, outfilepath, genetextfilepath)


# EXAMPLE CODE

Functional Frequency Spectra Workflow:

vcf<-load_vcf(filepathtovcf)

functional.spectra<-lof_funct(vcf, outfilecsv, optionaloutcsv)

LoF.funct.
