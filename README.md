# SimpleM
Calculate SimpleM threshold for GWAS significance as described in https://doi.org/10.1002/gepi.20430 and recommended in https://doi.org/10.1002/tpg2.20077.

# You will need:
Genotype data in numeric (0,1,2) format, where the number at each locus represents the number of copies of the reference allele for that taxon at that locus. The genotype file should be formatted with taxa names in column names and SNP names in the row names, and the name of each chromosome should be in a column labelled "chrom".

# Note for large data
Large matrices slow down SimpleM substantially. To increase computational efficiency, split each chromosome into smaller pieces and run SimpleM on each. Then, sum Meff across all pieces. See `simpleM_efficient.R` for an example using 10,000 markers at a time.
