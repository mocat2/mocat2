# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# FUNCTIONS FOR MOCAT ANALYZE PIPELINE

# define load function
cat ('Loading functions...\n')
load.function <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep=""))) 
    cat (' OK!\n')
  } else {
    cat('Updating packages...')
    update.packages(repos="http://mirrors.softliste.de/cran/", ask=FALSE)    
    cat (' OK!\n')
    cat('Installing packages...')
    eval(parse(text=paste("install.packages('", x, "', repos=\"http://mirrors.softliste.de/cran/\")", sep=""))) 
    cat (' OK!\n')
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep="")))
    cat (' OK!\n')
  } 
}

# Define load.bioconductor
load.bioconductor <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep=""))) 
    cat (' OK!\n')
  } else {
    cat('Updating packages...')
    update.packages(repos="http://mirrors.softliste.de/cran/", ask=FALSE)    
    cat (' OK!\n')
    cat('Installing packages...')
    source("http://www.bioconductor.org/biocLite.R")
    eval(parse(text=paste("biocLite('", x, "', ask=FALSE)", sep=""))) 
    cat (' OK!\n')
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep="")))
    cat (' OK!\n')
  } 
}

# Define annotate results
annotate.results <- function(x) {
  if (LOAD.FEATUE.ANNOTATIONS == 'yes') {
    for (name in names(FEATURE.ANNOTATIONS)) {
      cat(paste('Pasting ', name, ' annotations...\n', sep=''))
      for (result in names(x)) {
        oldcolnames <- colnames(x[[result]])
        annotated <- t(t(FEATURE.ANNOTATIONS[[name]][rownames(FEATURE.ANNOTATIONS[[name]]) %in% rownames( x[[result]] ),]))
        missing.annotation <- rownames(x[[result]])[!(rownames(x[[result]]) %in% rownames(annotated))]
        not.annotated <- matrix(nrow=length(missing.annotation),ncol=dim(annotated)[2], data='')
        rownames(not.annotated) <- missing.annotation
        all <- t(t(rbind(annotated, not.annotated)))
        all <- t(t(all[,!colSums(all == '') == nrow(all)])) # removes empty columns
        if (dim(all)[2]>0){
          colnames(all) <- paste(name, colnames(all), sep=':')
        }
        x[[result]] <- cbind(as.data.frame(x[[result]]) , as.data.frame(all[rownames( x[[result]] ),]))
        colnames(x[[result]]) <- c(oldcolnames, colnames(all))
      }      
    }
  } else {
    cat('No annotations loaded.\n')
  }
  return(x)
}

# Kristoffer's downsampling
forslund.downsample <- function (x, target, lengths) {
  # 'x' is a matrix of base counts for each
  # sample, included as rows
  #
  # 'target' is the approximate target number of bases
  # after downsampling
  #
  # 'lengths' is a column vector of average read length
  # for the corresponding samples
  #
  # test case:
  #
  # data = matrix (c (150,300,450,600),nrow = 1, ncol=4)
  # downsampleReads (data,c(750),c(75))
  #
  # for each sample, a downsampling factor is defined
  # as what would be needed to bring the total down to the
  # target
  #
  # for each gene/marker gene count, replace it by
  # <read length> * (sum of <count/read length> instances of
  # uniform randomization to the downsampling factor)
  # however, if there is not an integer number of pseudoreads,
  # the final addition must use a reduced downsampling factor
  # first compute sample sizes
  sampleSizes <- rowSums (x)
  # then compute downsampling factors
  downsamplingFactors <- rep (target, length = nrow (x)) / sampleSizes
  # how many pseudoreads?
  pseudoreads <- x / lengths
  wholeReads <- floor (pseudoreads)
  fractionReads <- pseudoreads - wholeReads
  # then for each sample and gene...
  for (i in 1:nrow (x)) {
    for (j in 1:ncol (x)) {
      downsampledReads <- sum (runif (wholeReads) < downsamplingFactors [i]) + sum (runif (1) < downsamplingFactors [i] * fractionReads)
      x [i, j] <- downsampledReads * lengths [i]
    }
  }
  x
}


# Calculate ROC curve for individual features, used in deseq, edger, wilcox so far
ROC <- function(w, x, y, z) {
  # Example: ROC(rownames(result),DATA,conditions,condB)
  # w = vector of feature names to produce curves for
  # x = DATA with all the features and samples, on the format row=samples, col=features
  # y = vector with different conditions
  # z = the condition that should be considered as "disease", in this case it should normally be "y" if the two states are eg. cancer (y) or no cancer (n)
  
  load(ROCR)
  l = vector('numeric', length(y)) - 1
  l[y == z] = 1
  
  for (W in w) {
    cat (paste("Plotting ROC curve for", W, "\n"))
    p = t(DATA)[W,]
    rocr.pr = prediction(as.numeric(p), as.numeric(l), c(-1,1))
    rocr.pf = performance(rocr.pr, 'tpr', 'fpr')
    auroc <- unlist(attributes(performance(rocr.pr, 'auc'))$y.values)
    plot(rocr.pf, main=paste(W, ": AUC", signif(auroc,3)))
  }
  
  for (W in w) {
    cat (paste("Plotting combined ROC curve for", W, "\n"))
    p = t(DATA)[W,]
    rocr.pr = prediction(as.numeric(p), as.numeric(l), c(-1,1))
    rocr.pf = performance(rocr.pr, 'tpr', 'fpr')
    auroc <- unlist(attributes(performance(rocr.pr, 'auc'))$y.values)
    plot(rocr.pf, main="ROC Curves", add=T)
  }
}


# Error bars in barplots
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Generate Fig 3 in Qin et al 2012 T2D
# The result table should be on the form  name statistic p.value p.adj condA_mean condB_mean
fig3qin <- function (result) {
  signif.features <- rownames(result[result[,'p.value'] <= 0.05,])
  condB.assoc.features <- intersect(rownames(result[result[,'condB_mean'] > result[,'condA_mean'],]), signif.features)
  condA.assoc.features <- intersect(rownames(result[result[,'condB_mean'] < result[,'condA_mean'],]), signif.features)
  
  #start from sample sum
  
  # wrong
  means <- c(mean(result[condA.assoc.features, 'condA_mean']), mean(result[condB.assoc.features, 'condB_mean']))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
