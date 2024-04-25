# SimpleM
Calculate SimpleM threshold for GWAS significance as described in https://doi.org/10.1002/gepi.20430 and recommended in https://doi.org/10.1002/tpg2.20077.

# You will need:
See notes at beginning of each script for required input files. 

# Which script to use?
For relatively small data (~10,000 markers), `simpleM_calculation.R` can easily run the whole genome at once. However, large matrices slow down SimpleM substantially. Therefore, to increase computational efficiency, I added `simpleM_efficient.R` as an option to run each chromosome separately, in chunks, as well as including other edits to improve efficiency over the base `simpleM_calculation.R`. 
