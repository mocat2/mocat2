# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC loads data into R object and plots basic overview
# PRODUCES pdf
# PRODUCES RData

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
# load.function('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages

#### PROGRAM SPECIFIC ####
# Get feature
cat('Plotting...\n')

pdf(file=paste(FIGURES.DIR, '/', 'overview.downsample=', TYPE, '..pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))

plot(sort(rowSums(DATA)), xlab='samples', ylab='total feature sum', main='Total feature sum of samples', pch=19)
hist(rowSums(DATA), ylab='samples', xlab='total feature sum', main='Total feature sum of samples', col='black')
plot(sort(rowSums(DATA)/sum(DATA)), main='Fraction of total feature sum', xlab='samples', ylab='fraction', col='blue', pch=19)
hist(rowSums(DATA)/sum(DATA), main='Fraction of total feature sum', xlab='fraction', ylab='number of samples', col='blue')
plot(sort(colSums(DATA)), xlab='features', ylab='total feature sum', main='Total feature sum of samples', col='red', pch=19)
hist(colSums(DATA), ylab='features', xlab='total features sum', main='Total sample sum of features', col='red')

dev.off()

# END #