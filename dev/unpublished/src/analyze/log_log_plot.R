# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC create linear models for Anita's funcitonal stuff. NOTE! important notes in the R file how data has to be formatted!
# REQUIRE metadata
# REQUIRE groups|1
## REQUIRE condA
## REQUIRE condB

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


# DEBUG / CREATE SCRIPT
#setwd('/g/bork5/kultima/MM/analysis/263RefGeneCatalog.padded/anita.samples/2014Oct27_135926.prepare_data_file.ko.scaled.base/data/')
#load('session.RData')
#feature <- 'Subject-Day-Method'


# Prepare
pdf(file=paste(FIGURES.DIR, '/', 'log10-log10-plot', '.', TYPE, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))
feature <- GROUPS[1]
conditions <- METADATA[,feature]
uc <- unique(conditions)
uc <- uc[!(uc %in% NA)]

# Run
for (condition in uc) {
  subset <- DATA[METADATA[,feature] %in% condition,]
  combinations <- combn(c(1:dim(subset)[1]),2)
  
  # loop over combinations
  for (i in c(1:dim(combinations)[2])){
    subset2 <- prop.table(subset[combinations[,i],], 1)
    plot(
      x=log10(subset2[1,]), y=log10(subset2[2,]),
      xlab=rownames(subset2)[1],ylab=rownames(subset2)[2],
      main=paste("log10-log10 plot :", TYPE), pch=16
    )
    top10 <- sort(subset2[1,], decreasing = T)[1:20]
    for (j in c(1:20)) {
      text(-3, -4.6-j*0.2, paste(names(top10[j]), signif(top10[j],3), sep=': '), pos=4, col='red')
    }
  }
}


# End
dev.off()


# END #
