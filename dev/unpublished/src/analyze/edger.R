# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC uses edgeR to identify differentially expressed features
# PRODUCES xls
# REQUIRE metadata
# REQUIRE groups|1
# REQUIRE iTOL
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
# load.function('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
load.bioconductor('edgeR')

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
countTable <- t(round(DATA))

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
  condA <- pairs[1,i]
  condB <- pairs[2,i]
  cat('Performing exact test...\n')
  et <- exactTest(y, pair=pairs[,i])
  cat('Saving results in tables...\n')
  
  # This is the only line, except the require line above that is needed for making the iTOL tree
  top <- topTags(et, n=10000, sort.by="p.value")$table
  top <- top[top[,'FDR'] < 0.05,]
  iTOLmultibar[[paste('log_fold_change_', paste(et$comparison, collapse="-"),'_FDR_is_significant' ,sep='')]] <- t(t(t(t(top))[,'logFC']))
  
  # Save results
  RESULTS[[paste('top1000p-value.comparison=', paste(et$comparison, collapse="."), '.', TYPE, sep='')]] <- topTags(et, n=1000, sort.by="p.value")
  RESULTS[[paste('top1000logFC.comparison=', paste(et$comparison, collapse="."), '.', TYPE, sep='')]] <- topTags(et, n=1000, sort.by="logFC")
  
  # ROC cruve
  pdf(file=paste(FIGURES.DIR, '/', 'edger.roc-curves.downsample=', DOWNSAMPLE, '.', TYPE, '.',condA, '.vs.', condB, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))
  ROC(rownames(top),DATA,conditions,condB)
  dev.off()
  
  # Annotate results
  RESULTS <- annotate.results(RESULTS)  
  
}

# ----------------------------- # edgeR PIPELINE

# END #
