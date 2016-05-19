# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# DESC Enterotyping. Clusters samples an returns PCoA plots, etc. Note that this will be done on abundances, eg sample sums are 1 (even if no normalization method is selected)
# PRODUCES pdf

# load packages
load.function('cluster')
load.function('clusterSim')
load.function('pvclust')

data <- t(DATA.WITH.MINUS1)
data <- t(t(data)/colSums(data)) # Want to use fractions, calculated BEFORE removing -1
data <- data[,!(colnames(data) == '-1'), drop=F] # remove -1
remove_noise <- 0 # not used at all, set to 0

pdf(file=paste(FIGURES.DIR, '/', 'enterotyping.', TYPE, '.pdf', sep=''), width=10, height=10, title=paste('Enterotyping of',TYPE))

### FUNCTIONS ###
cat('Loading functions...\n')
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  return(Matrix_1)
}


### EXECUTE ###
cat('Calculate JSD distances...\n')
data.dist=dist.JSD(data)
data.cluster=pam.clustering(data.dist, k=3)
nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index", main='Optimal number of clusters')

cat('Calculate optimal clusters...\n')
optimal <- which(max(nclusters, na.rm=T)==nclusters)
data.cluster=pam.clustering(data.dist, k=optimal)

# We don't remove any noise
#data.denoized=noise.removal(data, percent=remove_noise)
data.denoized <- data

cat('Running Between-Class Analysis...\n')
tmp <- t(data)
colnames(tmp) <- NULL
obs.pca=dudi.pca(data.frame(tmp), scannf=F, nf=k-1)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(rownames(obs.bet$ls)), grid=F, sub='Between-Class Analysis')

cat('Running PCoA Analysis...\n')
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=optimal)
s.class(obs.pcoa$li, fac=as.factor(rownames(obs.pcoa$li)), grid=F, sub='PCoA Analysis')

dev.off()
pdf(file=paste(FIGURES.DIR, '/', 'enterotyping.dendrogram.', TYPE, '.pdf', sep=''), width=20, height=10, title=paste('Enterotyping of',TYPE))
for (type in c('single', 'average', 'complete')) {
  cat(paste('Clustering for',type, 'linkage...\n', sep=' '))
  #OLD FUNCTION# plot(hclust(data.dist, method=type), main=paste('Sample Clustering using JSD distance and',type,'Linkage'))
  load.function(ClassDiscovery)  
  plotColoredClusters(hclust(data.dist, method=type), colnames(data), COLORS, cex = 0.7, line = -1, main=paste('Sample Clustering using JSD distance and',type,'linkage'), sub="")
  a <- pvclust(log10((data+sqrt(data*data+1))/2), method.hclust=type, nboot=100) # this is done on log10 of the data
  plot(a, main=paste('Sample Clustering using pvclust (BS=100)'))
  pvrect(a, main=paste('Sample Clustering using pvclust (BS=100)'))
}
dev.off()
