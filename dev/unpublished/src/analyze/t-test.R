# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC Calculates the t-test for taxo level over the specified group
# REQUIRE metadata
# REQUIRE groups|1

# Available variables:
# GROUPS = a vector with the selected metadata fields.
#          For the t-test script this is only one value.
#          (we know this, as we wanted only one input specified in the REQUIRE field)
# DATA
# DATA.WITH.MINUS1
# function.sourced
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages

#### PROGRAM SPECIFIC ####
cat("NOTE THIS PROGRAM IS ONLY A DUMMY PROGRAM TO TEST GROUPING CURRENTLY. SORRY!")
feature <- GROUPS[1]
file.result <- 1;
save <- paste('t-test.', feature, sep='')
#### PROGRAM SPECIFIC ####

#### GENERIC FOR ALL FUNCTIONS ####
file.result <- t(as.matrix(file.result))                       # Transpose
rownames(file.result) <- TYPE                                  # Add rownames
results[[save]] <- rbind(results[[save]], file.result)   # Add to results list
#### GENERIC FOR ALL FUNCTIONS ####



# END #
