# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# DESC Calculates the Shannon and Simpson measures for each sample using Shini's method. This should only be applied on curated.species when these are mOTUs AND on scaled insert counts. Any pre-processing is actualy discarded because the data is loaded from the files in the COGs folder. It's recommended to use sclaed inserts and mOTUs.
# PRODUCES xls

# Available variables:
# DATA
# DATA.WITH.MINUS1
# function.sourced
# THE DATA matrix has the format rows=samples, columns=data for each taxonomic level

# Available functions:
# load.function('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
load.function('vegan')

# This variable will be used
# MOTU_BASE_NAME

# Initialize
cat('Initilizing...\n')
sample.num <- dim(DATA)[1]
OTUcal.spe <- list()
OTUcal.spe.nz <- list()
MGs <- c('COG0012',
         'COG0016',
         'COG0018',
         'COG0172',
         'COG0215',
         'COG0495',
         'COG0525',
         'COG0533',
         'COG0541',
         'COG0552')

cat('Loading COGs...\n')
for (i in MGs){
  cat(paste('   ', i, '...\n', sep=''))
  OTUcal.spe[[i]] <- read.table(gzfile(paste(MOTU_BASE_NAME, i, 'gz', sep='.')), header=T, sep="\t", check.names=F, quote="")
  OTUcal.tmp.colnames <- colnames(OTUcal.spe[[i]])
  OTUcal.tmp <- cbind(OTUcal.spe[[i]][,1],floor(OTUcal.spe[[i]][-1]))
  colnames(OTUcal.tmp) <- OTUcal.tmp.colnames
  OTUcal.spe[[i]] <- OTUcal.tmp
  OTUcal.spe[[i]] <- OTUcal.spe[[i]][,-1]
  OTUcal.spe.nz[[i]] <- OTUcal.spe[[i]][rowSums(OTUcal.spe[[i]]) !=0,]
}

cat('Calculating...\n')
#Eco-indexes without downsampling
OTUcal.spe.nz.t <- list()
OTUcal.spe.nz.t.pres <- list()
OTUcal.spe.nz.t.relf <- list()
OTUcal.spe.nz.t.relspe <- list()
OTUcal.spe.nz.t.relsam <- list()
OTUcal.spe.nz.t.rich <- list()
OTUcal.spe.nz.t.chao1 <- list()
OTUcal.spe.nz.t.shannon.entropy <- list()
OTUcal.spe.nz.t.shannon.diversity <- list()
OTUcal.spe.nz.t.simpson.diversity <- list()
OTUcal.spe.nz.t.shannon.even <- list()
OTUcal.spe.nz.t.simpson.even <- list()
for (i in MGs){
  OTUcal.spe.nz.t[[i]] <- t(OTUcal.spe.nz[[i]])
  OTUcal.spe.nz.t.chao1[[i]] <- estimateR(floor(OTUcal.spe.nz.t[[i]]))["S.chao1",]
  OTUcal.spe.nz.t.pres[[i]] <- apply(OTUcal.spe.nz.t[[i]] > 0, 2, sum)
  OTUcal.spe.nz.t.relf[[i]] <- sort(100*(apply(OTUcal.spe.nz.t[[i]] > 0, 2, sum)/nrow(OTUcal.spe.nz.t[[i]]))) #relative frequencies of species in all samples
  OTUcal.spe.nz.t.relspe[[i]] <- decostand(OTUcal.spe.nz.t[[i]],"total",MARGIN=2) #relative abundances per species; check: apply(OTUcal.spe.nz.t.relspe[[i]],2,sum)
  OTUcal.spe.nz.t.relsam[[i]] <- decostand(OTUcal.spe.nz.t[[i]],"total",MARGIN=1) #relative abundances per sample; check: apply(OTUcal.spe.nz.t.relsam[[i]],1,sum)
  OTUcal.spe.nz.t.rich[[i]] <- specnumber(OTUcal.spe.nz.t[[i]]) #same as: apply(OTUcal.spe.nz.t[[i]] > 0, 1, sum)
  OTUcal.spe.nz.t.shannon.entropy[[i]] <- diversity(OTUcal.spe.nz.t[[i]],index="shannon")
  OTUcal.spe.nz.t.shannon.diversity[[i]] <- exp(OTUcal.spe.nz.t.shannon.entropy[[i]])
  OTUcal.spe.nz.t.simpson.diversity[[i]] <- diversity(OTUcal.spe.nz.t[[i]],index="inv")
  OTUcal.spe.nz.t.shannon.even[[i]] <- OTUcal.spe.nz.t.shannon.diversity[[i]]/OTUcal.spe.nz.t.rich[[i]]
  OTUcal.spe.nz.t.simpson.even[[i]] <- OTUcal.spe.nz.t.simpson.diversity[[i]]/OTUcal.spe.nz.t.rich[[i]]
}

#were OTUcal.spe.nz is the floored, scaled insert count matrix after removing zero rows, and
#rownames(MGs) are the 10 marker genes. For Shannon diversity index, you actually only need OTUcal.spe.nz.t.shannon.entropy.

#After that, I do:

mat <- matrix(unlist(OTUcal.spe.nz.t.shannon.entropy),10,sample.num,byrow=T)
colnames(mat) <- names(unlist(OTUcal.spe.nz.t.shannon.entropy[[i]]))
rownames(mat) <- MGs
mean.shannon.index <- apply(mat,2,mean)
mean.shannon.index.sd <- apply(mat,2,sd)

#where mean.shannon.index will be Shanon diversity indices, and
#mean.shannon.index.sd the stdevs across 10 marker genes.

# Paste into results matrix
result <- cbind(mean.shannon.index, mean.shannon.index.sd)
colnames(result) <- c('Shannon Index', 'SD Shannon Index')

#### GENERIC FOR ALL FUNCTIONS ####
save <- 'a-diversity.shannon'
RESULTS[[save]] <- result  # Add to results list
#### GENERIC FOR ALL FUNCTIONS ####

# END #
