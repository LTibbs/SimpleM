# Laura Tibbs-Cortes
# Jan 30, 2024

# Calculating SimpleM threshold for high-density SNP data (>10,000 markers per chromosome)
# by separating each chromosome into multiple chunks

# Also contains other edits to make code more efficient for larger data.

# For more information about simpleM, see doi 10.1002/gepi.20430

# You will need: 
# NOTE: this is different from the regular simpleM_caculation.R requirements, for increased speed per run.
# Genotype data in "012" format, which can be produced by vcftools (https://vcftools.sourceforge.net/man_latest.html#:~:text=OUTPUT%20OTHER%20FORMATS-,%2D%2D012,-This%20option%20outputs) 
# This consists of three files, with suffixes .012, .012.indv, and .012.pos.
# There should be separate files for each chromosome (e.g., genotype_file_Chr1.012, genotype_file_Chr2.012, etc.)
# Missing data is NOT allowed, so use imputation or filtering before simpleM to address this. 


# Load packages:
# NOTE: this assumes that you have the tidyverse and data.table packages installed.
# If you do not, first run: install.packages(c("tidyverse", "data.table"))
library(data.table)
library(tidyverse)

print(sessionInfo())


# Settings: See doi 10.1002/gepi.20430 for more information, or use defaults below
pca_cutoff_perc <- 0.995 # set cutoff point for eigenvalues
sig.level <- 0.05 # set significance level

# Other settings:
marker.num <- 10000 # how many markers per "chunk" of the chromosome. Can increase if more computing resources available.
chr.arg <- "Chr1" # set chromosome name for this round; will need to re-run for each chromosome (alternatively, could make a loop to run all, but this way can easily run each chromosome separately in parallel)

# make the simpleM function:
Meff.simpleM=function(genotype, cor_r, eigen_path, pca_cutoff_perc) {
  
  num_of_snps <- ncol(genotype) # how many SNPs are there?
  
  eigen_values <- eigen(cor_r, only.values = TRUE)$values # calculate and pull eigenvalues
  
  # In this case, not writing eigenvalues to save time per run, but can choose to do so by un-commenting
  # fwrite(as.tibble(eigen_values), eigen_path)
  
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
  indv <- read_tsv(paste0("genotype_file_", chr.arg,".012.indv"),
                   col_names="Taxa")
  
  pos <- fread(paste0("genotype_file_", chr.arg,".012.pos"),
               col.names=c("Chromosome", "Position"), data.table=F)
  
  # Here, use fread to get the number of columns in the genotype file without reading the whole thing in
 num.cols <- fread(paste0("genotype_file_", chr.arg,".012"), sep="\t", data.table = F,
                  stringsAsFactors = F, header=F, verbose=F,
                  nrows=0)
 ncol(num.cols)
 
  # use the simpleM function on my data, one chromosome at a time:
  Meff.list <- vector("list", length=ceiling(ncol(num.cols)/marker.num) )
  # now start iterating here, read in every marker.num (default 10,000) markers
  for(k in 1:ceiling(nrow(pos)/ marker.num)) {# get how many 10,000s of SNPs there are
    
    start.snp <- (((k-1)*marker.num)+1) # start at 1 or 10,001 or...
    end.snp <- min(nrow(pos), k*marker.num) # end snp number for this chunk

    # output info for current round
    print(paste("chr", chr.arg, start.snp, end.snp, nrow(pos)))
    
    # read in current chunk of genotype file:
    num <- fread(paste0("genotype_file_", chr.arg,".012"), sep="\t", data.table = F,
                      stringsAsFactors = F, header=F, verbose=F,
                      na.strings="-1",
                      select=c((start.snp+1):(end.snp+1)) # need to do +1 here because num has a non-meaningful first column that would need to be dropped
    )

    # make into expected file format (based on GAPIT GM/GD format):
    GM <- pos[start.snp:end.snp,] %>%
      mutate(Name=paste(Chromosome, "_", Position, sep="")) %>%
      dplyr::select(Name, Chromosome, Position)
    
    colnames(num) <- GM$Name
    GD <- cbind(indv, num)
    colnames(GD) <- c("Taxa", GM$Name)
    
    # sanity checks
    stopifnot(all.equal(colnames(GD),c("Taxa",GM$Name)))
    stopifnot(nrow(indv)==nrow(num))
    
    # clear up workspace
    rm(num)
    gc()
    
    # remove Taxa so can make into numeric matrix:
    GD.matrix <- as.matrix(GD[,-1])
    stopifnot(!any(is.na(GD.matrix))) # confirm there are no missing genotypes; SimpleM does not work with missing genotypes
    # format genotype matrix
    class(GD.matrix) <- "numeric"# this is now the same format as t.genotype in simpleM_calculation.R
    
    rm(GD, GM) # tidy workspace
    
    # SimpleM requires me to remove monomorphic SNPs:
    num.SNPs <- vector(length = ncol(GD.matrix))
    for (j in 1:ncol(GD.matrix)) {
      num.SNPs[j] <- length(unique(GD.matrix[,j]))
    }
    GD.matrix <- GD.matrix[,num.SNPs>1]
    
    # use modified version of COMBAT's ld.Rsquare function:
    cor_G <- cor(as.matrix(GD.matrix))
    stopifnot(ncol(cor_G)==ncol(GD.matrix))
    
    # run simpleM function
    M_eff_G <- Meff.simpleM(genotype=GD.matrix, cor_r = cor_G,
                            eigen_path=paste0("simpleM.eigenvalues.chr", chr.arg,"_pt",k, ".csv"), pca_cutoff_perc = pca_cutoff_perc)
    # save results from this chromosome
    Meff.list[[k]] <- M_eff_G
    print(paste("finshed with ", file.name," chromosome", chr.arg, "part", k, "Meff_G is", M_eff_G))
    rm(M_eff_G, cor_G)
  }
  
  Meff.G <- sum(do.call(cbind, (Meff.list))[1,]) # find total number of independent SNPs across chromosome
  
  threshold <- sig.level/Meff.G # use Bonferroni correction on number of independent SNPs
  print(paste0(file.name, " chr ", chr.arg," Meff.G is: ", Meff.G)) #output Meff.G (total number of independent SNPs) for one chromosome
  

  fwrite(Meff.list, paste0("SimpleM/",file.name,".chr", chr.arg,".simpleM.MeffG.cutoff", pca_cutoff_perc, ".csv")) # output results
  
# To calculate final, genome-wide threshold:
  # sum up the total number of independent SNPs (Meff.G) from each chromosome to get the genome-wide number of independent SNPs (denoted genome.wide.Meff.G)
  # Then, threshold = sig.level/genome.wide.Meff.G