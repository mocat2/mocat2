# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC Clusters samples and produces heat map
# PRODUCES pdf
# PRODUCES jpg
# REQUIRE metadata
# SETTINGS : hclust_function_sample [complete,average,single,median,centroid] {CHECK} {average} (Select clustering function for samples)
# SETTINGS : hclust_function_feature [complete,average,single,median,centroid] {CHECK} {average} (Select clustering function for feature)
# SETTINGS : dist_function_sample [cor.dist,spearman.dist,dist.euclidean,dist.maximum,dist.manhattan,dist.canberra,dist.binary,dist.minkowski] {CHECK} {dist.euclidean} (Select distance function for samples)
# SETTINGS : dist_function_feature [cor.dist,spearman.dist,dist.euclidean,dist.maximum,dist.manhattan,dist.canberra,dist.binary,dist.minkowski] {CHECK} {dist.euclidean} (Select distance function for features)
# SETTINGS : metadata_to_plot [any of the metadata columns] {NOCHECK} {} (Select which metadata columns to plot under the heatmap)

# Available variables:
# GROUPS = a vector with the selected metadata fields.
#          For the t-test script this is only one value.
#          (we know this, as we wanted only one input specified in the REQUIRE field)
#          NOTE! the REQUIRE group,1 fiels must be ABOVE any of the REQUIRE condA, REQUIRE condB fields.
# Use this to retrieve a setting: SETTINGS['specified setting above', 'setting']
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

### LOAD PACKAGES ###
load.function('gplots')
load.bioconductor('bioDist')
load.bioconductor("Heatplus")

## OOMPA is annoying because the main installation is OOMPA, but package to reuiqre is ClassDiscovery
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if (!is.installed('ClassDiscovery')) {
  load.function('mclust')
  source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
  oompainstall('ClassDiscovery')  
}
load.bioconductor('ClassDiscovery')


#### PROGRAM SPECIFIC ####
cat('Initializing...\n')


### DEFINE AND PROCESS VARIABLES ###
DATA <- t(DATA)
if (SETTINGS['metadata_to_plot', 'setting'] != '') {
  metadata.to.plot <- unlist(strsplit(SETTINGS['metadata_to_plot', 'setting'],','))
  metadata.to.plot <- gsub(pattern='^\\s*', replacement='', x=metadata.to.plot, perl=T)
  metadata.to.plot <- gsub(pattern='\\s*$', replacement='', x=metadata.to.plot, perl=T)
  annotation.all <- droplevels(as.data.frame(METADATA[, metadata.to.plot ]))
  # Removes any annotation column with fewer than two levels, this is required by annHeatmap2
  columns <- 0
  counter <- 0
  colnames <- vector()
  for (i in 1:dim(annotation.all)[2]) {
    if (length(levels(annotation.all[,i])) >= 2) {
      colnames <- c(colnames, colnames(annotation.all[i]))
      if (columns == 0) {
        columns <- columns + 1
        annotation <- as.vector(annotation.all[,i])
      } else {
        annotation <- cbind(annotation, as.vector(annotation.all[,i]))
      }   
    }
  }
  annotation <- as.matrix(annotation)
  colnames(annotation) <- colnames
  rownames(annotation) <- rownames(annotation.all)
  annotation <- droplevels(as.data.frame(annotation))
} else {
  annotation <- data.frame()
}

### DEFINE FUNCTIONS ###
# Define clustering functions
clust.sample <- function(x) {
  hclust(d=x, method=SETTINGS['hclust_function_sample', 'setting'])
}
clust.feature <- function(x) {
  hclust(d=x, method=SETTINGS['hclust_function_feature', 'setting'])
}

#Define distance functions
dist.feature <- function(x) {
  split <- unlist(strsplit( SETTINGS['dist_function_feature', 'setting'], "\\."))
  if (split[1] == 'dist') {
    eval(parse(text=paste( 'dist', "(x, method='", split[2],  "')", sep="")))
  } else {
    eval(parse(text=paste( SETTINGS['dist_function_feature', 'setting'], "(x)", sep="")))
  }
}
dist.sample <- function(x) {
  split <- unlist(strsplit( SETTINGS['dist_function_sample', 'setting'], "\\."))
  if (split[1] == 'dist') {
    eval(parse(text=paste( 'dist', "(x, method='", split[2],  "')", sep="")))
  } else {
    eval(parse(text=paste( SETTINGS['dist_function_sample', 'setting'], "(x)", sep="")))
  }
}

### PROCESSING STEPS ###

# Heatmap
cat('Generating heatmap...\n')
heatmap <- annHeatmap2(DATA,
                  ann=list(Col=list(data=annotation)),
                  dendrogram=list(
                    Row=list(
                      distfun=dist.feature,
                      clustfun=clust.feature),
                    Col=list(
                      distfun=dist.sample,
                      clustfun=clust.sample)
                  ),
                  col=blueyellow,
                  scale="row"
)

cat('Plotting heatmap...\n')
pdf(file=paste(FIGURES.DIR, '/', 'cluster.downsample=', DOWNSAMPLE, '.',  TYPE, '.pdf', sep=''), width=10, height=10, title=paste('Heatmap of',TYPE))
plot(heatmap)
dev.off()
jpeg(file=paste(FIGURES.DIR, '/', 'cluster.downsample=', DOWNSAMPLE, '.', TYPE, '.jpg', sep=''), width=1200, height=1200)
plot(heatmap)
dev.off()

# END #