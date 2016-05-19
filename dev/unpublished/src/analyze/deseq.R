# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC uses the DESeq package to find diferentially expressed features
# PRODUCES xls
# PRODUCES pdf
# PRODUCES csv
# REQUIRE metadata
# REQUIRE groups|1
# REQUIRE condA
# REQUIRE condB
# REQUIRE iTOL
# SUPPRESS normalize
# SETTINGS : fitType [parametric,local] {CHECK} {parametric} (Set fitType for estimateDispersions)

# Available variables:
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
cat('Loading packages...\n')
load.bioconductor('DESeq')

#### PROGRAM SPECIFIC ####
# Get feature
cat('Initializing...\n')
condA <- SETTINGS['condA','setting']
condB <- SETTINGS['condB','setting']
feature <- GROUPS[1]

# Get levels of feature
levels <- levels(as.factor(METADATA[,feature]))

# Check that number of levels is 2 and store them
#if (length(levels) != 2) {
#  cat('Levels: ')
#  cat(levels)
#  cat('\nPlease enter the 2 levels of your choice:\n')
#  condA <- readLines("stdin", n=1)    
#  condB <- readLines("stdin", n=1)
#} else {
#  condA <- levels[1]
#  condB <- levels[2]
#}

# DESeq specific variables needed to run DESeq
conditions <- METADATA[,feature]

# Subset conditions
conditions <- conditions[conditions == condA | conditions == condB]

# Use only conditions that are not NA
conditions <- conditions[!is.na(conditions)]

# Subset DATA to be those samples only that have feature == condA or feature == condB
DATA <- DATA[names(conditions),]

# Create count table
countTable <- t(round(DATA))

pdf(file=paste(FIGURES.DIR, '/', 'deseq.downsample=', DOWNSAMPLE, '.', TYPE, '.',condA, '.vs.', condB, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))

cat('Creating DESeq object...\n')
cds <- newCountDataSet(countTable, conditions)


# ----------------------------- # DESeq PIPELINE
# Pipeline
### R code from vignette source 'vignettes/DESeq/inst/doc/DESeq.Rnw'
### code chunk number 1: options
options(digits=3)
### code chunk number 10: headcounts1
#head( counts(cds) )
### code chunk number 11: estimateSizeFactors
cat('Processing [1/14]...\n')
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
### code chunk number 12: headcounts2
#head( counts( cds, normalized=TRUE ) )
### code chunk number 13: estimateDispersions
cat('Processing [2/14]...\n')
cds <- estimateDispersions( cds, fitType=SETTINGS['fitType', 'setting'] )
### code chunk number 14: str
cat('Processing [3/14]...\n')
#str( fitInfo(cds) )
### code chunk number 15: fitUntr
cat('Processing [4/14]...\n')
plotDispEsts <- function( cds )
{
  plot(
    rowMeans( counts( cds, normalized=TRUE ) ),
    fitInfo(cds)$perGeneDispEsts,
    pch = '.', log="xy" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
### code chunk number 16: figFit
cat('Processing [5/14]...\n')
plotDispEsts( cds )
### code chunk number 17: DESeq.Rnw:297-298
cat('Processing [6/14]...\n')
all(table(conditions(cds))==2)
### code chunk number 18: head
#head( fData(cds) )
### code chunk number 19: str
cat('Processing [7/14]...\n')
#str( fitInfo(cds) )
### code chunk number 20: nbt1
cat('Processing [8/14] (this will take some time)...\n')
res <- nbinomTest( cds, condA, condB)
### code chunk number 21: nbt2
#head(res)
### code chunk number 22: figDE
cat('Processing [9/14]...\n')
plotDE <- function( res )
  plot(
    res$baseMean,
    res$log2FoldChange,
    log="x", pch=20, cex=.3,
    col = ifelse( res$padj < .1, "red", "black" ) )
plotDE( res )
### code chunk number 23: histp
cat('Processing [10/14]...\n')
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
### code chunk number 24: ressig1
cat('Processing [11/14]...\n')
resSig <- res[ res$padj < 0.1, ]
### code chunk number 25: ressig2
#head( resSig[ order(resSig$pval), ] )
### code chunk number 26: ressig3
#head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
### code chunk number 27: ressig4
#head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )
### code chunk number 28: writetable (eval = FALSE)
## #not run
## write.table( res, file="results.txt" )
### code chunk number 29: ncu
cat('Processing [12/14]...\n')
#ncu <- counts( cds, normalized=TRUE )[ , conditions(cds)=="untreated" ]
### code chunk number 30: MArepl
#plot( rowMeans(ncu), log2( ncu[,2] / ncu[,1] ), pch=".", log="x" )
cat('Processing [13/14]...\n')
table <- res[res[,'pval'] < 0.1 | res[,'padj'] < 0.1,]
# ----------------------------- # DESeq PIPELINE

# Removes NA rows
table <- table[!is.na(table[,1]),]
rownames(table) <- table[,'id']
table <- table[,2:8]
cat('Processing [14/14]...\n')
file.result <- table
dev.off()
save <- paste('deseq.', TYPE, '.',  condA, '.vs.', condB, sep='')
#### PROGRAM SPECIFIC ####

#padj <- as.matrix(as.numeric(file.result[,'padj']<0.05))
#rownames(padj) <- rownames(file.result)
#colnames(padj) <- 'padj_le_0.05'

padj <- file.result[file.result[,'padj']<0.05,6:7]
if (dim(padj)[1]>0){
  padj[,1:2] <- "#063C92"
}

inf<-t(t(t(t(file.result))[!is.finite(file.result[,'log2FoldChange']),'log2FoldChange']))
inf[,1] <- "#DEB009"

# ROC cruve
#pdf(file=paste(FIGURES.DIR, '/', 'deseq.roc-curves.downsample=', DOWNSAMPLE, '.', TYPE, '.',condA, '.vs.', condB, '.pdf', sep=''), width=10, height=10, title=paste('Plots for',TYPE))
#ROC(rownames(file.result),DATA,conditions,condB)
#dev.off()

RESULTS[[save]] <- file.result   # Add to results list

iTOLcolorscale[['p_adjusted_is_significant']] <- padj
iTOLmultibar[['log2_fold_change']] <- t(t(t(t(file.result))[is.finite(file.result[,'log2FoldChange']),'log2FoldChange']))
iTOLcolorscale[['log_fold_change_is_inf_or_neg_inf']] <- inf

RESULTS <- annotate.results(RESULTS)
#### GENERIC FOR ALL FUNCTIONS ####


# END #