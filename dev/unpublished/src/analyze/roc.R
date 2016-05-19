# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.







# NOT FULLY DEVELOPED!!!







# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC plots ROC curves for individual features
# PRODUCES xls
# PRODUCES pdf
# REQUIRE metadata
# REQUIRE groups|1
# REQUIRE condA
# REQUIRE condB
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








stop("Use the ROC function in the functions file")




## Use the ROC function in the functions file





#### PROGRAM SPECIFIC ####
# Get feature
cat('Initializing...\n')
feature <- GROUPS[1]
condA <- SETTINGS['condA','setting']
condB <- SETTINGS['condB','setting']

#Set PDF
pdf(file=paste(FIGURES.DIR, '/', 'wilcox.downsample=', DOWNSAMPLE, '.', TYPE, '.',condA, '.vs.', condB, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))

# Get conditions
cat('Settting up conditions...\n')
conditions <- METADATA[,feature]

# Subset conditions
conditions <- conditions[conditions %in% c(condA, condB)]

# Make conditions factor
conditions <- factor(conditions)

# Use only conditions that are not NA
conditions <- conditions[!is.na(conditions)]

# Subset DATA to be those samples only that have feature == condA or feature == condB
cat('Creating count data table...\n')
DATA <- DATA[names(conditions),]

# ----------------------------- # wilcox test pipeline
cat('Running wilcox pipeline...\n')
result <- matrix(ncol=3, nrow=(dim(DATA)[2]))
rownames(result) <- colnames(DATA)
colnames(result) <- c('statistic', 'p.value', 'p.adj')
for (i in c(1:dim(DATA)[2])) {
  result[i,] <- c(as.numeric(unlist(wilcox.test(DATA[conditions == condA,i], DATA[conditions == condB,i], alternative="two.sided", paired=F))[1:2]), -1)
}
result[,3] <- p.adjust(result[,2], method = "fdr")
result <- result[order(result[,3], na.last=NA), ]
hist(result[,2], col='red', main='P-values of Wilcoxon Test', xlab='p-value')
hist(result[,3], col='darkgreen', main='FDR Adjusted P-values of Wilcoxon Test', xlab='adjusted p-value')
iTOLmultibar[[paste('statistic_', condA, '_vs_', condB, '_FDR_is_significant' ,sep='')]] <- result[result[,3]<=0.05,,drop=F]
RESULTS[[paste('p-values.comparison=', condA, '.vs.', condB, '.', TYPE, sep='')]] <- result
save.image(file=paste(DATA.DIR,"/session.RData",sep=""))
RESULTS <- annotate.results(RESULTS)  

# ----------------------------- # wilcox test pipeline

dev.off()

# END #
