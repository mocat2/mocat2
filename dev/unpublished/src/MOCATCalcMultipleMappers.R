# This script post processes multiple mappers info files
# As input is a list of files, of which are on the format TAX tab TAX tab INTEGER

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
base <- args[1]

#base <- 'S.none'

# Required functions
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


# Load required packages
load.function('gdata')
load.function('shiny')
load.function('bigmemory')

cat('Calculate Multiple Mappers - R session :: load name files\n')
read.files <- as.matrix(read.table(paste(base, "read.files", sep="."), header=F))[,1]
insert.files <- as.matrix(read.table(paste(base, "insert.files", sep="."), header=F))[,1]
read.stats.files <- as.matrix(read.table(paste(base, "read.stats.files", sep="."), header=F))[,1]
insert.stats.files <- as.matrix(read.table(paste(base, "insert.stats.files", sep="."), header=F))[,1]
output.files <- as.matrix(read.table(paste(base, "output.files", sep="."), header=F))[,1]
samples <- as.matrix(read.table(paste(base,"samples", sep="."), header=F))[,1]


cat('Calculate Multiple Mappers - R session :: load data files\n')
read <- list()
insert <- list()
read.stats <- list()
insert.stats <- list()
for (i in 1:length(samples) ) {
  read[[samples[i]]] <- read.table(read.files[i], header=F, sep="\t", quote="", colClasses=c("character", "character", "double"))
  insert[[samples[i]]] <- read.table(insert.files[i], header=F, sep="\t", quote="", colClasses=c("character", "character", "double"))
  read.stats[[samples[i]]] <- read.table(read.stats.files[i], header=F, sep="\t", quote="", colClasses=c("character", "character", "double"))
  insert.stats[[samples[i]]] <- read.table(insert.stats.files[i], header=F, sep="\t", quote="", colClasses=c("character", "character", "double"))
}


cat('Calculate Multiple Mappers - R session :: get unique taxa\n')
all.names.read <- vector()
all.names.insert <- vector()
for (i in 1:length(samples) ) {
  all.names.read <- unique(c(all.names.read, as.vector(unique(c(read[[i]][,1], read[[i]][,2])))))
  all.names.insert <- unique(c(all.names.insert, unique(c(insert[[i]][,1], insert[[i]][,2]))))
}
all.names.read <- sort(unlist(all.names.read))
all.names.insert <- sort(all.names.insert)


cat('Calculate Multiple Mappers - R session :: generate matrices\n')
data.list.read <- list()
data.list.insert <- list()
for (i in 1:length(samples) ) {
  cat(paste('Calculate Multiple Mappers - R session :: generate matrices :: ', samples[i], '\n', sep=''))
  data.list.read[[samples[i]]] <- matrix(0, nrow=length(all.names.read), ncol=length(all.names.read))
#  data.list.read[[samples[i]]] <- big.matrix(0, nrow=length(all.names.read), ncol=length(all.names.read), type='integer',
#dimnames=list(all.names.read, all.names.read)
#                                             )

  colnames(data.list.read[[samples[i]]]) <- all.names.read
  rownames(data.list.read[[samples[i]]]) <- all.names.read


  
  for (j in 1:dim(read[[samples[i]]])[1]) {
    data.list.read[[samples[i]]][read[[samples[i]]][j, 1], read[[samples[i]]][j, 2]] <- data.list.read[[samples[i]]][read[[samples[i]]][j, 1], read[[samples[i]]][j, 2]] + read[[samples[i]]][j, 3]
  }
  
  data.list.insert[[samples[i]]] <- matrix(0, nrow=length(all.names.insert), ncol=length(all.names.insert))
#  data.list.insert[[samples[i]]] <- big.matrix(0, nrow=length(all.names.insert), ncol=length(all.names.insert), type='integer',
#                                               dimnames=list(all.names.insert, all.names.insert)
#                                               )
  colnames(data.list.insert[[samples[i]]]) <- all.names.insert
  rownames(data.list.insert[[samples[i]]]) <- all.names.insert
  for (j in 1:dim(insert[[samples[i]]])[1]) {
    data.list.insert[[samples[i]]][insert[[samples[i]]][j, 1], insert[[samples[i]]][j, 2]] <- data.list.insert[[samples[i]]][insert[[samples[i]]][j, 1], insert[[samples[i]]][j, 2]] + insert[[samples[i]]][j, 3]
  } 
}



cat('Calculate Multiple Mappers - R session :: generate final matrix\n')
sum.read <- Reduce('+', data.list.read)
sum.insert <- Reduce('+', data.list.insert)


cat('Calculate Multiple Mappers - R session :: calculating combinations\n')
top.list.read <- list()
top.list.insert <- list()
comb.read <- combn(colnames(sum.read),2)
comb.insert <- combn(colnames(sum.insert),2)
names.read <- paste(comb.read[1,], comb.read[2,], sep=':')
names.insert <- paste(comb.insert[1,], comb.insert[2,], sep=':')


for (i in 1:length(samples) ) {
  cat(paste('Calculate Multiple Mappers - R session :: find top multiple mappers :: ', samples[i], '\n', sep=''))
  
  # reads
  temp <- data.list.read[[samples[i]]]
  diag(temp) <- 0
  un.sum <- sort(unmatrix(temp)[names.read], decreasing=T)
  top.list.read[[samples[i]]] <- un.sum[un.sum >= 5]
  
  # inserts
  temp <- data.list.insert[[samples[i]]]
  diag(temp) <- 0
  un.sum <- sort(unmatrix(temp)[names.insert], decreasing=T)
  top.list.insert[[samples[i]]] <- un.sum[un.sum >= 5] 
}
  

cat('Calculate Multiple Mappers - R session :: find top multiple mappers for all samples\n')
# reads
temp <- sum.read
diag(temp) <- 0
un.sum <- sort(unmatrix(temp)[names.read], decreasing=T)
top.read <- un.sum[un.sum >= 5]

# inserts
temp <- sum.insert
diag(temp) <- 0
un.sum <- sort(unmatrix(temp)[names.insert], decreasing=T)
top.insert <- un.sum[un.sum >= 5]

# Set additional variables
sum.insert.interactions <- 0
sum.read.interactions <- 0
for (i in samples) {
  sum.read.interactions <- sum.read.interactions + read.stats[[i]][3,3]
  sum.insert.interactions <- sum.insert.interactions + insert.stats[[i]][3,3]
}


cat('Calculate Multiple Mappers - R session :: saving objects\n')
reads <- list()
inserts <- list()
reads[['top']] <- top.read
inserts[['top']] <- top.insert 
reads[['top_samples']] <- top.list.read
inserts[['top_samples']] <- top.list.read
reads[['sum']] <- sum.read
inserts[['sum']] <- sum.insert 
reads[['sum_samples']] <- data.list.read
inserts[['sum_samples']] <- data.list.read
reads[['stats']] <- read.stats
inserts[['stats']] <- insert.stats
reads[['interactions']]<- sum.read.interactions
inserts[['interactions']]<- sum.insert.interactions

save(samples, reads, inserts, file=paste(output.files[1], '.RData', sep=''))
reads_sum <- reads[['sum']]
inserts_sum <- inserts[['sum']]
save(reads_sum, inserts_sum, file=paste(output.files[1], '.sum.matrix.RData', sep=''))

cat('Calculate Multiple Mappers - R session :: completed\n')
