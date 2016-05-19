# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# DESC Calculates the relative ratio of two specified features
# REQUIRE featureA
# REQUIRE featureB
# PRODUCES xls
# PRODUCES csv

# Available variables:
# IMPORTANT NOTE! featureA and featureB are actually stored in condA and condB
# DATA
# DATA.WITH.MINUS1
# FUNCTION.SOURCED
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
# Nothing to load

# Calculate ratio
  
  #### PROGRAM SPECIFIC ####
  condA <- SETTINGS['condA','setting']
  condB <- SETTINGS['condB','setting']
  file.result <- DATA[,condA] / DATA[,condB]
  #### PROGRAM SPECIFIC ####
  
  #### GENERIC FOR ALL FUNCTIONS ####
  file.result <- t(as.matrix(file.result))                       # Transpose
  rownames(file.result) <- TYPE                                  # Add rownames
  save <- paste('ratio.', condA, '.over.', condB, sep='')        # Store filename
  save <- gsub(" ", "_", save)                                   # Fix spaces in filename
  RESULTS[[save]] <- rbind(RESULTS[[save]], file.result)         # Add to results list
  #### GENERIC FOR ALL FUNCTIONS ####
  


# END #

