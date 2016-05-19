# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC uses edgeR to identify differentially expressed features
# PRODUCES xls
# REQUIRE metadata
# REQUIRE groups|1
# SETTINGS : responses [comma separated list] {NOCHECK} {} (If you wish to use a subset of responses for the group, enter them comma separated)
# SETTINGS : edgeRnormalization [none,TMM,RLE,upperquartile] {CHECK} {TMM} (Select a method for normalization from the edgeR package that is applied to the count table)
# SUPPRESS normalize

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
load.bioconductor('multtest')
install.packages("ShotgunFunctionalizeR", contriburl="http://shotgun.math.chalmers.se")
library('ShotgunFunctionalizeR')

#### PROGRAM SPECIFIC ####
# Get feature
cat('Initializing...\n')
feature <- GROUPS[1]
responses <- unlist(strsplit(SETTINGS['responses', 'setting'],','))
responses <- gsub(pattern='^\\s*', replacement='', x=responses, perl=T)
responses <- gsub(pattern='\\s*$', replacement='', x=responses, perl=T)


# Get conditions
cat('Settting up conditions...\n')
conditions <- METADATA[,feature]

# Subset conditions
if (length(responses > 0)) {
  conditions <- conditions[conditions %in% responses]
}

# Make conditions factor
conditions <- factor(conditions)

# Use only conditions that are not NA
conditions <- conditions[!is.na(conditions)]

# Subset DATA to be those samples only that have feature == condA or feature == condB
cat('Creating count data table...\n')
DATA <- DATA[names(conditions),]

# Create count data
data <- t(DATA)

# ----------------------------- # edgeR PIPELINE
cat('Running edgeR pipeline...\n')
cat('Creating objects...\n')
y <- DGEList(counts=countTable,group=conditions)

# Normalize
if (SETTINGS['edgeRnormalization', 'setting'] != 'none') {
  cat(paste('Normalizing using ', SETTINGS['edgeRnormalization', 'setting'] ,'...\n', sep=''))
  y <- calcNormFactors(y, method=SETTINGS['edgeRnormalization', 'setting'])  
}

cat('Estimating common disp...\n')
y <- estimateCommonDisp(y)
cat('Estimating tagwise disp...\n')
y <- estimateTagwiseDisp(y)

# Get all combinations
pairs <- combn(levels(conditions),2)

for (i in 1:dim(pairs)[2]) {
  cat('Performing exact test...\n')
  et <- exactTest(y, pair=pairs[,i])
  cat('Saving results in tables...\n')
  
  RESULTS[[paste('top100p-value.comparison=', paste(et$comparison, collapse="."), '.', TYPE, sep='')]] <- topTags(et, n=100, sort.by="p.value")
  RESULTS[[paste('top100logFC.comparison=', paste(et$comparison, collapse="."), '.', TYPE, sep='')]] <- topTags(et, n=100, sort.by="logFC")
  RESULTS <- annotate.results(RESULTS)  
  
}

# ----------------------------- # edgeR PIPELINE

# END #
