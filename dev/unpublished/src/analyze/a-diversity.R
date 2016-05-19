# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# DESC Calculates the Shannon and Simpson measures for each sample
# PRODUCES xls

# Available variables:
# DATA
# DATA.WITH.MINUS1
# function.sourced
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load.function('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
load.function('vegan')

# Calculate Shannon and Simpson
for (measure in c("simpson", "shannon", "invsimpson")) {
  
  #### PROGRAM SPECIFIC ####
  file.result <- diversity(floor(DATA), index=measure)                  # Calculate, also floor the data
  save <- paste('a-diversity.', measure, sep='')
  #### PROGRAM SPECIFIC ####
  
  #### GENERIC FOR ALL FUNCTIONS ####
  file.result <- t(t(as.matrix(file.result)))                       # Transpose
  colnames(file.result) <- TYPE                                  # Add rownames
  RESULTS[[save]] <- cbind(RESULTS[[save]], file.result)   # Add to results list
  #### GENERIC FOR ALL FUNCTIONS ####
  
}

# END #
