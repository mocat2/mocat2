# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC generates iTOL tree with abundances of all taxa in the samples
## PRODUCES xls
## PRODUCES pdf
## REQUIRE metadata
## REQUIRE groups|1
## REQUIRE condA
## REQUIRE condB
# REQUIRE iTOL

# Available variables:
# SAMPLE.GROUPS = a list of the different sample groups defined by a selected 'groups|1'
# GROUPS = a vector with the selected metadata fields.
#          For the t-test script this is only one value.
#          (we know this, as we wanted only one input specified in the REQUIRE field)
#          NOTE! the REQUIRE group,1 fiels must be ABOVE any of the REQUIRE condA, REQUIRE condB fields.
# DATA
# DATA.WITH.MINUS1
# METADATA
# function.sourced
# TABLES.DIR
# FIGURES.DIR
# DATA.DIR
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages

#### PROGRAM SPECIFIC ####
# Get feature
cat('Running...\n')
#feature <- GROUPS[1]
#condA <- SETTINGS['condA','setting']
#condB <- SETTINGS['condB','setting']

#Set PDF
#pdf(file=paste(FIGURES.DIR, '/', 'wilcox.downsample=', DOWNSAMPLE, '.', TYPE, '.',condA, '.vs.', condB, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))

# Get conditions
#cat('Settting up conditions...\n')
#conditions <- METADATA[,feature]

# Subset conditions
#conditions <- conditions[conditions %in% c(condA, condB)]

# Make conditions factor
#conditions <- factor(conditions)

# Use only conditions that are not NA
#conditions <- conditions[!is.na(conditions)]

# Subset DATA to be those samples only that have feature == condA or feature == condB
#cat('Creating count data table...\n')
#DATA <- DATA[names(conditions),]



results <- as.matrix(colSums(DATA))

iTOLmultibar[['abundance']] <- results
RESULTS[['abundance']] <- results
RESULTS <- annotate.results(results)  

cat('Completed.\n')


# END #
