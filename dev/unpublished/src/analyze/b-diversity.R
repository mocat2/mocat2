# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# DESC Calculates beta diversity for the selected samples
# PRODUCES xls
# PRODUCES pdf
# Available variables:
# FIGURES.DIR
# TABLES.DIR
# TYPE (eg base.raw.order.solexaqa.allbest.l45.p95)
# DATA
# DATA.WITH.MINUS1
# FUNCTION.SOURCED (equals TRUE, if this functions has been run previously for different type settings)
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load.function('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
load.function('vegan')

# Calculate beta diversity

#### PROGRAM SPECIFIC ####

# Calculate
file.result <- betadiver(DATA)                                  # Calculate the diversity

pdf(file=paste(FIGURES.DIR, '/', 'b.diversity.downsample=', TYPE, '.pdf', sep=''), width=10, height=10, title=paste('Beta Diversity of',TYPE))
plot(file.result, pch=19)
dev.off()

# Save a,b,c in different files (by adding them as different data frames in the result list)
for (counter in c(1:length(attributes(file.result)$names))) {
  letter <- attributes(file.result[counter])$names              # Get the a,b,c from $a, $b $c in the object
  file <- paste(TYPE, 'b-diversity', letter, sep='.')           # Make filename
  RESULTS[[file]] <- as.matrix(file.result[[counter]])          # Add to results list
}
#### PROGRAM SPECIFIC ####


# END #
