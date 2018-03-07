# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC calculates the wilcoxon statistic for each feature between to states, eg cancer and non cancer
# PRODUCES xls
# PRODUCES pdf
# REQUIRE metadata
# REQUIRE groups|1
# REQUIRE condA
# REQUIRE condB
# REQUIRE iTOL
# SETTINGS : additional_analyses [none,plot_abundance_sum] {CHECK} {none} (Select additional analyses to perform, plot_abundance_sum generates Fig 3 in Qin et al 2012 on T2D)

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
cat('Initializing...\n')
feature <- GROUPS
condA <- names(SAMPLE.GROUPS)[1]
condB <- names(SAMPLE.GROUPS)[2]

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

# Loop over the different DATA types
for (data_type in names(DATAlist)) {
  
  cat(paste('Processing',data_type,':\n'), sep='')
  
  # Subset DATA to be those samples only that have feature == condA or feature == condB
  cat('Creating count data table...\n')
  DATA <- DATAlist[[data_type]]
  DATA <- DATA[names(conditions),]
  cat(paste('DATA has ', dim(DATA)[1], ' samples and ', dim(DATA)[2], ' features. Group by metadata column ',feature,' and conditionas are ',condA,' and ',condB,'\n', sep=''))
  
  # ----------------------------- # wilcox test pipeline
  cat('Running wilcox pipeline...\n')
  result <- matrix(ncol=7, nrow=(dim(DATA)[2]))
  rownames(result) <- colnames(DATA)
  colnames(result) <- c('statistic', 'p.value', 'p.adj.fdr', 'p.adj.bonferroni', 'p.adj.BY', 'condA_Mean', 'condB_mean')
  for (i in c(1:dim(DATA)[2])) {
    result[i,] <- c(as.numeric(unlist(wilcox.test(DATA[conditions == condA,i],
                                                  DATA[conditions == condB,i],
                                                  alternative="two.sided", paired=F))[1:2]),
                    -1,-1,-1,
                    mean(DATA[conditions == condA,i]),
                    mean(DATA[conditions == condB,i]))
  }
  result[,3] <- p.adjust(result[,2], method = "fdr")
  result[,4] <- p.adjust(result[,2], method = "bonferroni")
  result[,5] <- p.adjust(result[,2], method = "BY")
  result <- result[order(result[,3], na.last=NA), ]
  hist(result[,2], col='red', main=paste(data_type, ' - P-values of Wilcoxon Test', sep=''), xlab='p-value')
  hist(result[,3], col='blue', main=paste(data_type, ' - FDR Adjusted P-values of Wilcoxon Test', sep=''), xlab='adjusted p-value')
  hist(result[,4], col='darkgreen', main=paste(data_type, ' - Bonferroni Adjusted P-values of Wilcoxon Test', sep=''), xlab='adjusted p-value')
  hist(result[,5], col='orange', main=paste(data_type, ' - Benjamini Yekutieli Adjusted P-values of Wilcoxon Test', sep=''), xlab='adjusted p-value')
  
  # Perhaps this will be included later, but for now I will leave it out
  #iTOLmultibar[[paste(data_type, '_statistic_', condA, '_vs_', condB, '_FDR_is_significant' ,sep='')]] <- result[result[,3]<=0.05,,drop=F]
  #iTOLmultibar[[paste(data_type, '_statistic_', condA, '_vs_', condB, '_Bonferroni_is_significant' ,sep='')]] <- result[result[,4]<=0.05,,drop=F]
  #iTOLmultibar[[paste(data_type, '_statistic_', condA, '_vs_', condB, '_BY_is_significant' ,sep='')]] <- result[result[,5]<=0.05,,drop=F]
  
  RESULTS[[paste(data_type, '.p-values.comparison.', condA, '.vs.', condB, '.', TYPE, sep='')]] <- result
}
RESULTS <- annotate.results(RESULTS)  
# ----------------------------- # wilcox test pipeline

dev.off()

# END #
