# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC create linear models for Anita's funcitonal stuff. NOTE! important notes in the R file how data has to be formatted!
# REQUIRE metadata
## REQUIRE groups|1
## REQUIRE condA
## REQUIRE condB

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



########################################## 
# Packages and libraries
##########################################
load.function(lme4)
load.function(ggplot2)
load.function(reshape)
load.bioconductor("vsn")
load.function(ape)
load.function(stringr)


##############################
# TESTING, code form Anita to load the original data form her MS
# workingDirectory = '/g/bork5/avoigt/TimeSeries'
# setwd(workingDirectory)
# AnnotatedotuData <- read.table('mm.402.motu.is.annotated.mOTU.clusters.gz',quote="", header=T, row.names=1, check.names=FALSE, sep='\t', comment='#')
# Samples <- read.table('All_samples',quote="", header=F, check.names=FALSE, sep='\t', comment='#') # Load Samples to be included, has to be saved as vector, (no row.names=1 otherwise it expects data in column 1) # Header has to be F otherwise it takes a sample as header
# Samples <- as.vector(Samples[,1]) #  Make a vector
# subset <- t(AnnotatedotuData[,colnames(AnnotatedotuData) %in% Samples])            # Select wanted samples in table based on vector, both ways work
# Day7_samples <- names(which(sapply(row.names(subset), function(x) {return(strsplit(x,'-')[[1]][3])})==7))
# Day7_samples <- Day7_samples[!Day7_samples == 'minniemouse-11-7-0']                         # Remove the samples which has no method comparison
# Day7_samples <- Day7_samples[!Day7_samples == 'daisy-11-7-0']
# Day7_samples <- Day7_samples[!Day7_samples == 'alien-9-7-0']                         # Remove the samples which has no method comparison
# Day7_samples <- Day7_samples[!Day7_samples == 'halbarad-9-7-0']
# Day7_samples <- Day7_samples[!Day7_samples == 'scavenger-9-7-0']                         # Remove the samples which has no method comparison
# Day7_samples <- Day7_samples[!Day7_samples == 'peacemaker-9-7-0']
# Day7_samples <- Day7_samples[!Day7_samples == 'bugkiller-9-7-0']
# Day7 <-as.data.frame(subset(subset, row.names(subset) %in% Day7_samples))   # Take all replicates from fraction table
# DATAo <- DATA
# DATA <- Day7
##############################



########################################## 
# Data
##########################################
cat ('Running script...\n')
subset <- rownames(DATA) # the samples to include, for now all of them

# REMOVE FOR FINAL VERSION
subset <- names(which(sapply(subset, function(x) {return(strsplit(x,'-')[[1]][3])})==392)) # RIGHT NOW WE WORK ONLY ON DAY 7

subdata <- DATA[rownames(DATA) %in% subset,] # get subset
prop <- prop.table(as.matrix(subdata), 1) # make proportions

# NOTE 1
# The names has to be on the format XXX-A-B-C, where XXX is the name of the individual.
# because here we select species to test perindividual somehow.
# from Anita's file: Check for Species that are overlapping whithin 1 indivudual and unique accross individuals, better than just all the overlapping species
list <- sapply(subset, function(x) {return(strsplit(x,'-')[[1]][1])})
toTest <- list()
i <- 1
for (samples in unique(list)){
  toTest[[i]] <- names(which(colSums( prop[list %in% samples,] ==0)==0))
  i <- i + 1
}
toTest <- unique(unlist(toTest)) # here we get a final list of feaures to test
prop <- prop[,toTest] # here we create the subset, these are the only features we look at

# NOTE 2
# Paul said this could be done, but in the final manuscript Anita didn't run this.
# "Paul: asin only should be fine... Variance stabilization:"
# prop <- asin(prop)

prop <- as.data.frame(prop) # cerate data frame so we can add protocoal and day information
prop$protocol <- sapply(rownames(prop),function(x) {return(strsplit(x,'-')[[1]][2])})
prop$patient <- sapply(rownames(prop),function(x) {return(strsplit(x,'-')[[1]][1])})
combinations <- combn(unique(prop$protocol), 2) # create combinations to test



########################################## 
# Model
##########################################
for (j in 1:dim(combinations)[2]) {
  subprop <- prop[prop$protocol %in% combinations[,j],] # get subset
  meltsub <- melt(subprop, id=c("patient","protocol")) # melt it
  names(meltsub) <- c("Patient", "Protocol", "Features", "Count") # rename 
  features = levels(meltsub$Features)               # take all (unique) motus using levels
  N <- length(features)                           # How many motus are there
  pS <- rep("NA", N)                           # repeat/put NA N times (number of motus), later we can see whether all p-values have been computed
  i <- 1                                       # start with 1 to loop through motus
  for (a_feature in features) {                                                                          # loops through all motus
    current_data <- subset (meltsub, Features == a_feature)                                                  # variable = species column, 
    current_data$Protocol <- factor (current_data$Protocol)                                                 # make numbers to factors so that it is not taken as numeric
    current_data$Patient <- factor (current_data$Patient)
    if (length (unique (current_data$Count)) > 1) {                                                    # all those species that are different from each other meaning not all with the same (0) value
      modelFull <- lmer (data = current_data, Count ~ Protocol + Patient - 1 + (1 | Patient))        # full outer model with fixed effect for patient (without patient) and protocol plus random effect for Patient
      modelReduced <- lmer (data = current_data, Count ~  Patient - 1 + (1 | Patient))               # reduced model without the factor that we want to look at, to see whether this factor has an effect
      aP <- anova (modelReduced, modelFull)$"Pr(>Chisq)"[2]                                       # takes the p value from the anova test output
      pS [i] <- aP                                                                                # write the anova output into ps for all motus
    }
    i <- i + 1      
  } # ignore the warnings?
  pS_fdr <- p.adjust (pS, method = "BH")                        # Does Benjamini hochberg korrektion,Makes a list of all p values that are significant but not NA
  Z <- cbind (features [(!is.na (pS_fdr)) & pS_fdr <= 2],        # Gives all the species with fdr <2 that should be tested further
              pS_fdr [(!is.na (pS_fdr)) & pS_fdr <= 2], 
              pS [(!is.na (pS_fdr)) & pS_fdr <= 2]) 
  
  pS_by <- p.adjust (pS, method = "BY")                        # Does Benjamini hochberg korrektion,Makes a list of all p values that are significant but not NA
  Z2 <- cbind (features [(!is.na (pS_by)) & pS_by <= 2],        # Gives all the species with fdr <2 that should be tested further
               pS_by [(!is.na (pS_by)) & pS_by <= 2], 
               pS [(!is.na (pS_by)) & pS_by <= 2]) 
  
  pS_bonferroni <- p.adjust (pS, method = "bonferroni")                        # Does Benjamini hochberg korrektion,Makes a list of all p values that are significant but not NA
  Z3 <- cbind (features [(!is.na (pS_bonferroni)) & pS_bonferroni <= 2],        # Gives all the species with fdr <2 that should be tested further
               pS_bonferroni [(!is.na (pS_bonferroni)) & pS_bonferroni <= 2], 
               pS [(!is.na (pS_bonferroni)) & pS_bonferroni <= 2]) 
  
  
  Z <- cbind(paste(combinations[,j], collapse="-"), Z)
  colnames(Z) <- c('Protocal X vs Y', 'Feature', 'Benjamini Hochberg (BH) corrected P-value', 'P-value')
  
  Z2 <- cbind(paste(combinations[,j], collapse="-"), Z2)
  colnames(Z2) <- c('Protocal X vs Y', 'Feature', 'Benjamini Yekutieli (BY) corrected P-value', 'P-value')
  
  Z3 <- cbind(paste(combinations[,j], collapse="-"), Z3)
  colnames(Z3) <- c('Protocal X vs Y', 'Feature', 'Bonferroni corrected P-value', 'P-value')
  
  RESULTS[[ paste('Protocol_BH_correction_', paste(combinations[,j], collapse="-") , sep="") ]] <- Z
  RESULTS[[ paste('Protocol_BY_correction_', paste(combinations[,j], collapse="-") , sep="") ]] <- Z2
  RESULTS[[ paste('Protocol_bonferroni_correction_', paste(combinations[,j], collapse="-") , sep="") ]] <- Z3
  
  cat(paste("\n\n============= OUTPUT : PROTOCOL ",
            paste(combinations[,j], collapse="-"),
            " MODEL FULL", " =============\n", sep="") )
  print(summary (modelFull))
  
  cat(paste("\n\n============= OUTPUT : PROTOCOL ",
            paste(combinations[,j], collapse="-"),
            " MODEL REDUCES", " =============\n", sep="") )
  print(summary(modelReduced))
  
  
  ########################################## 
  # Check significant Species
  ##########################################
  for (a_feature in c(Z[,2], Z2[,2], Z3[,2])) {                                                                          # loops through all motus
    current_data <- subset (meltsub, Features == a_feature)                                                  # variable = species column, 
    current_data$Protocol <- factor (current_data$Protocol)                                                 # make numbers to factors so that it is not taken as numeric
    current_data$Patient <- factor (current_data$Patient)
    
    modelFull <- lmer (data = current_data, Count ~ Protocol + Patient - 1 + (1 | Patient))        # full outer model with fixed effect for patient (without patient) and protocol plus random effect for Patient
    modelReduced <- lmer (data = current_data, Count ~  Patient - 1 + (1 | Patient))               # reduced model without the factor that we want to look at, to see whether this factor has an effect
    
    residSS.Full = sum((current_data$Count - predict(modelFull) )^2) 
    residSS.Red = sum((current_data$Count - predict(modelReduced) )^2) 
    totalSS = var(current_data$Count)*(length(current_data$Count)-1)
    
    expVar.Full = 1 - residSS.Full/totalSS
    expVar.Red = 1 - residSS.Red/totalSS
    
    n = dim(current_data)[1]      # Sample size
    pFull = 2                     # Regessoren 2 fuer Protocol und Patient
    pRed = 1                      # Regressor fuer Patient
    
    Adj.Full = 1-(1 - expVar.Full)* (n-1)/(n-pFull-1)
    Adj.Red = 1-(1 - expVar.Red)* (n-1)/(n-pRed-1)
    aP <- anova (modelReduced, modelFull)                                      # takes the p value from the anova test output
    
    cat(paste("\n\n============= OUTPUT : PROTOCOL",
              paste(combinations[,j], collapse="-"),
              "CHECK SIGNIFICANT SPECIES: ", a_feature, " =============\n", sep="") )
    print (aP)
    
    x <- c(a_feature,i, expVar.Full, expVar.Red, Adj.Full, Adj.Red, current_data$Count)
    y <- as.vector(apply( as.matrix(current_data[,1:2]) , 1 , paste , collapse = "-" ))
    names(x) <- c('Feature', 'i', 'expVar.Full', 'expVar.Red', 'Adj.Full', 'Adj.Red', y)
    tmp <- str_replace_all(string=a_feature, pattern=" ", repl="_")
    #RESULTS[[paste('check_signf_feature_', tmp, '_', paste(combinations[,j], collapse="-"), sep='') ]] <- x
    i <- i + 1      
  }
  
} # end combinations

# END #
