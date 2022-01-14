# Laura Tibbs-Cortes
# Sept 7, 2021

# Calculating simpleM threshold for given genotype data
# For more information about simpleM, see doi 10.1002/gepi.20430

# You will need:
# genotype data in numeric (0,1,2) format, where the number at each locus represents the number of copies of the reference allele for that taxon at that locus.
# Genotype file should be formatted with taxa names in column names and SNP names in the row names.
# The name of each chromosome should be in a row labelled "chrom".

# Load packages:
# NOTE: this assumes that you have the tidyverse and data.table packages installed.
# If you do not, first run: install.packages(c("tidyverse", "data.table"))
library(data.table)
library(tidyverse)

print(sessionInfo())

# Settings: See doi 10.1002/gepi.20430 for more information, or use defaults below
pca_cutoff_perc <- 0.995 # set desired cutoff point for eigenvalues
sig.level <- 0.05 # set desired significance level

# make the simpleM function:
Meff.simpleM=function(genotype, cor_r, eigen_path, pca_cutoff_perc) {
  
  num_of_snps <- ncol(genotype) # how many SNPs are there?
  
  eigen_values <- eigen(cor_r, only.values = TRUE)$values # calculate and pull eigenvalues

  # save these eigenvalues! 
  # Can use them for altering threshold based on pca_cutoff_perc without re-running whole script
  fwrite(as.tibble(eigen_values), eigen_path)
  
  sum_eigen_values <- sum(eigen_values)
  eigen_values <- sort(eigen_values, decreasing = TRUE)
  
  M_eff_G <- 1
  for(k in 1:num_of_snps){
    temp <- sum(eigen_values[1:k])/sum_eigen_values # determine the proportion of variance explained by this # of eigenvalues
    if(temp >= pca_cutoff_perc){ # determine if the proportion of variance explained is above the provided cutoff
      M_eff_G <- k # save the effective marker value
      break
    }
  }
  return(M_eff_G)
}


# run simpleM -------------------------------------------------------------

# read in genotype data 
# for required genotype format: see notes at beginning of script
a=fread("genotype_file.csv")

# format genotype data:
a2=a[,-c(1:11)] # remove extra information columns: chromosome number, SNP position, etc. 
# BE CAREFUL-- hmp genotype files typically have 11 informational columns, 
# but different files will have different numbers of informational columns, 
# so you need to check your own data to find the appropriate number of columns to remove.

a2=as.matrix(a2) # format as matrix

stopifnot(!any(is.na(a2))) # make sure there are no missing genotypes; SimpleM does not work with missing genotypes

# format genotype matrix
class(a2) <- "numeric"
t.genotype=t(a2) # transpose

# use the simpleM function on my data, one chromosome at a time:
Meff.list <- vector("list", length=0)
for(i in unique(a$chrom)) { # loop through each chromosome present in the genotype data
  
  # pull genotype for chromosome:
  genotype <- t.genotype[,a$chrom==i]
  
  # SimpleM requires me to remove monomorphic SNPs:
  num.SNPs <- vector(length = ncol(genotype))
  for (j in 1:ncol(genotype)) {
    num.SNPs[j] <- length(unique(genotype[,j]))
  }
  genotype <- genotype[,num.SNPs>1]
  
  # use modified version of COMBAT's ld.Rsquare function:
  cor_G <- cor(as.matrix(genotype)) # need SNPs in row names, samples in column names for correlation
  stopifnot(ncol(cor_G)==ncol(genotype))
  
  # run simpleM function
  M_eff_G <- Meff.simpleM(genotype=genotype, cor_r = cor_G, 
                          eigen_path=paste0("simpleM.eigenvalues.chr", i,".csv"), pca_cutoff_perc = pca_cutoff_perc)
  # save results from this chromosome
  Meff.list[[i]] <- M_eff_G
  print(paste("finshed with chromosome", i, "Meff_G is", M_eff_G))
  rm(M_eff_G)
  
}

fwrite(Meff.list, paste0("simpleM.MeffG.cutoff", pca_cutoff_perc, ".csv")) # output results

Meff.G <- sum(do.call(cbind, (Meff.list))[1,]) # find total number of independent SNPs across genome

threshold <- sig.level/Meff.G # use Bonferroni correction on number of independent SNPs
print(threshold) #output SimpleM threshold

# Calculate other thresholds: ---------------------------------------------

# If desired, can use the eigenvalues calculated above to calculate new thresholds based on new pca_cutoff_perc and sig.level

pca_cutoff_perc=0.99 # set new cutoff point for eigenvalues
sig.level <- 0.05 # set new significance level

# read genotype data 
a=fread("genotype_file.csv")

Meff.list <- vector("list", length=0)
eigen.list <- vector("list", length=0)
for(i in unique(a$chrom)) { 
  
  # read in previous results
  eigen_path <- paste0("simpleM.eigenvalues.chr", i,".csv")
  eigen.list[[i]] <- fread(eigen_path)
  
  # calculate simpleM's MeffG:
  eigen_values <- sort(eigen.list[[i]]$value, decreasing = TRUE)
  sum_eigen_values <- sum(eigen_values)
  
  M_eff_G <- 1
  for(k in 1:nrow(eigen.list[[i]])){
    temp <- sum(eigen_values[1:k])/sum_eigen_values # determine the proportion of variance explained by this # of eigenvalues
    if(temp >= pca_cutoff_perc){
      M_eff_G <- k
      break
    }
  }
  Meff.list[[i]] <- M_eff_G
  print(paste("finshed with chromosome", i, "Meff_G is", M_eff_G))
  rm(M_eff_G)
}
fwrite(Meff.list, paste0("simpleM.MeffG.cutoff", pca_cutoff_perc, ".csv")) # output M_eff_G for each chromosome

# pull MeffG from file to calculate threshold
Meff.csomewise <- fread(paste0("simpleM.MeffG.cutoff", pca_cutoff_perc, ".csv"))
Meff.G <- sum(Meff.csomewise[1,])

threshold <- sig.level/Meff.G # use Bonferroni correction on number of independent SNPs
print(threshold) # output SimpleM threshold
