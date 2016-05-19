# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# This script does the following:
# - Installs packages if needed, loads them
# - Loads data from pre-created MOCAT SETTINGS files
# - Runs selected external scripts
# - Summarizes and prints

# To create new functions, please use the scripts in the
# MOCAT/unpublished/src/analyze/ folder as examples, and also save new scripts there.


# ALSO THIS SCRIPT CAN TAKE SETTINGS AND REQUIRE AS REQUIREMEWNTS
# SETTINGS : filter [numeric] {NOCHECK} {0} (Specify filter cutoff for min sum over samples, if 0 disabled)
# SETTINGS : filter_fraction [numeric] {NOCHECK} {0} (Specify filter for top fraction of samples to keep, eg. 0.1 keeps top 10% samples, if 0 disabled)
# SETTINGS : filter_low_abundant_features [numeric] {NOCHECK} {0} (If >0, features that are above X, in fraction, in at least 1 sample will be kept)
# SETTINGS : filter_prevalence_of_features [numeric] {NOCHECK} {0} (If >0, features have to have abundance >0 in at least X fraction of samples, if run on genes, use 0.01 to specify that genes have to have abundance >0 in at least 1% of the samples)
# SETTINGS : normalize [none,decostand:total,decostand:max,decostand:freq,decostand:normalize,decostand:standardize,decostand:pa,decostand:hellinger,decostand:log] {CHECK} {none} (Select a method for normalization)


# LOAD FILES
cat ('Loading settings...\n')
SETTINGS <- as.matrix(read.table('R.settings', header=T, row.names=1, sep="\t"))
load.functions <- SETTINGS['load_functions','setting']
DATA.DIR <- SETTINGS['data_dir','setting']

# Load functions
source(load.functions)
cat ('Initializing history and saving commands...\n')
load.function('TeachingDemos')
txtStart(file=paste(DATA.DIR,"/session.History",sep=""), cmdfile=paste(DATA.DIR,"/session.Rerun",sep=""))

cat ('Continuing loading settings...\n')
WORKING.DIR <- SETTINGS['working_dir','setting']
TABLES.DIR <- SETTINGS['tables_dir','setting']
FIGURES.DIR <- SETTINGS['figures_dir','setting']
DOWNSAMPLE  <- SETTINGS['downsample','setting']
LOAD.GROUPS <- SETTINGS['load_groups','setting']
LOAD.METADATA <- SETTINGS['load_metadata','setting']
LOAD.FEATUE.ANNOTATIONS <- SETTINGS['load_feature_annotations','setting']

setwd(WORKING.DIR)
files <- as.matrix(read.table('R.files', header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char=""))
functions <- as.matrix(read.table('R.functions', quote="", header=T, row.names=1, sep="\t", comment.char=""))

# LOAD FILES IF REQUIRED
if (LOAD.GROUPS == 'yes') {
  cat ('Loading groups...\n')
  GROUPS <- as.matrix(read.table('R.groups'))[,1]
}
if (LOAD.METADATA != 'no') {
  cat ('Loading metadata...\n')
  METADATA.FULL <- as.matrix(read.table(LOAD.METADATA, quote="", header=T, row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char=""))  
}


# Load feature annotations
FEATURE.ANNOTATIONS <- list()
if (LOAD.FEATUE.ANNOTATIONS == 'yes') {
  cat ('Loading feature annotations:\n')
  FEATURE.ANNOTATION.FILES <- as.matrix(read.table('R.feature_annotations'))
  for (i in c(1:dim(FEATURE.ANNOTATION.FILES)[1])) {
    cat (paste('        Loading ', FEATURE.ANNOTATION.FILES[i,1], ' annotations...\n', sep=''))
    FEATURE.ANNOTATIONS[[ FEATURE.ANNOTATION.FILES[i,1] ]] <- as.matrix(read.table(FEATURE.ANNOTATION.FILES[i,2],
                                                                                   header=F, row.names=1, sep="\t", quote="", check.names=F, fill=T,
                                                                                   na.strings="NA", comment.char="", nrows = as.numeric(FEATURE.ANNOTATION.FILES[i,4]),
                                                                                   colClasses=rep('character', as.numeric(FEATURE.ANNOTATION.FILES[i,3]))))
  }
}


# CALCULATIONS

# Prepare
RESULTS <- list()

# Process function by function
for (execute in functions[,'location']){
  FUNCTION.SOURCED <- FALSE
  
  # Process file, one by one
  for (index in c(1:dim(files)[1])) {
    
    # Load info
    FILE <- files[index, 'file']
    TYPE <- files[index, 'all']
    ZIPPED <- files[index, 'zipped']
    
    if (files[index, 'tax_level'] == 'mOTU') {
      path.split <- unlist(strsplit(FILE, '/'))
      name <- path.split[length(path.split)]
      name <- sub('.gz$', '', name, perl=T)
      path <- paste(paste(path.split[-length(path.split)], collapse='/'), '/COGs/', name, sep='') 
      MOTU_BASE_NAME <- path
    }
    
    # Print status
    cat(paste(" == PROCESSING ", TYPE, " ==\n", sep=''))
    
    # Load data
    cat('Loading data...\n')
    if (ZIPPED == 'no') {
      DATA.WITH.MINUS1 <- t(as.matrix(read.table(FILE, quote="", header=T, row.names=1, sep="\t", check.names=F),drop=F))
    } else if (ZIPPED == 'yes') {
      DATA.WITH.MINUS1 <- t(as.matrix(read.table(gzfile(FILE), quote="", header=T, row.names=1, sep="\t", check.names=F),drop=F))
    } else {
      stop('File not specified as zipped or not zipped. Fatal error.')
    }
      
    mode(DATA.WITH.MINUS1) <- 'numeric'
    DATA <- DATA.WITH.MINUS1[,!(colnames(DATA.WITH.MINUS1) == '-1' | colnames(DATA.WITH.MINUS1) == 'mapped' | colnames(DATA.WITH.MINUS1) == 'unassigned' | colnames(DATA.WITH.MINUS1) == 'sum_annotated' | colnames(DATA.WITH.MINUS1) == 'sum_not_annotated' | colnames(DATA.WITH.MINUS1) == 'sum_not_annotated_and_annotated' | colnames(DATA.WITH.MINUS1) == 'mapped_inserts_or_bases' | colnames(DATA.WITH.MINUS1) == 'total_inserts_or_bases'), drop=F]
    #DATA.ORIGINAL <- DATA # here we could save a copy of the data before any processing
    
    # Downsample
    # SOME NOTES: Downsampling is done on the DATA frame without -1.
    # Assume we have actual species, then the minus isn't included, which make the most sense.
    # Assume we have OTUs then it doesn't matter, the DATA and DATA.WITH.MINUS1 will be identical.    
    if (DOWNSAMPLE == 'none') {
    } else if (DOWNSAMPLE == 'rrarefy.min') {
      load.function('vegan')
      DATA <- floor(DATA)
      cat(paste("Downsampling ", TYPE, " using rrarefy.min...\n", sep=''))
      DATA <- rrarefy(DATA, min(rowSums(floor(DATA))))
    } else if (DOWNSAMPLE == 'rrarefy.top90percent') {
      load.function('vegan')
      top90 <- names(sort(rowSums(floor(DATA))))[c((ceiling( dim(DATA)[1] / 10)):dim(DATA)[1])]
      DATA <- DATA[top90,, drop=F]
      DATA <- floor(DATA)
      cat(paste("Downsampling ", TYPE, " using rrarefy.top90percent...\n", sep=''))
      DATA <- rrarefy(DATA, min(rowSums(floor(DATA))))
    } else if (DOWNSAMPLE == 'forslund') {
      if (LOAD.METADATA == 'no') {
        stop ('Using downsampling method forslund requires your script to load METADATA with a column called average_read_length.\n')
      }
      if ('average_read_length' %in% colnames(METADATA.FULL) == FALSE) {
        stop ('METADATA is loaded, but using downsampling method \'forslund\' requires your script to load METADATA with a column called average_read_length.\n')
      }
      cat(paste("Downsampling ", TYPE, " using \'forslund\'...\n", sep=''))
      samples.with.values <- names(!is.na(METADATA.FULL[rownames(DATA),'average_read_length'])) [as.vector(!is.na(METADATA.FULL[rownames(DATA),'average_read_length']))]
      DATA <- forslund.downsample(DATA[samples.with.values,], as.numeric(SETTINGS['downsample_size','setting']), as.numeric(METADATA.FULL[samples.with.values,'average_read_length']))
    } else if (DOWNSAMPLE == 'rrarefy') {
      load.function('vegan')      
      DATA <- DATA[rowSums(DATA)>as.numeric(SETTINGS['downsample_size','setting']),, drop=F]      
      DATA <- floor(DATA)
      cat(paste("Downsampling ", TYPE, " using rrarefy...\n", sep=''))
      DATA <- rrarefy(DATA, as.numeric(SETTINGS['downsample_size','setting']))
    } else {
      stop ('Unknown downsampling method selected.\n')
    } 
    
    # Normalize
    if ('normalize' %in% rownames(SETTINGS)) {
      if (SETTINGS['normalize', 'setting'] != 'none') {
        cat(paste('Normalizing using ',SETTINGS['normalize', 'setting'],'...\n', sep=''))
        normalize.method <- unlist(strsplit(SETTINGS['normalize', 'setting'],':'))
        if (normalize.method[1] == 'decostand') {
          load.function('vegan')
          DATA <- decostand(DATA, normalize.method[2], 1)
        }
      }
    }
    
    # Filter out columns (features) with low total sums across samples
    if ('filter' %in% rownames(SETTINGS)) {
      if (as.numeric(SETTINGS['filter', 'setting']) == 0) {
      } else if (as.numeric(SETTINGS['filter', 'setting']) > 0) {
        cat('Filtering using integer and min sum over samples...\n')
        DATA <- DATA[,colSums(DATA) >= as.numeric(SETTINGS['filter', 'setting']), drop=F]
      } else {
        stop('Filter settings are incorrect.\n')
      }
    }
    
    # Filter using a fraction of samples
    if ('filter_fraction' %in% rownames(SETTINGS)) {
      if (as.numeric(SETTINGS['filter_fraction', 'setting']) == 0) {
      } else if (as.numeric(SETTINGS['filter_fraction', 'setting']) > 0 && as.numeric(SETTINGS['filter_fraction', 'setting']) < 1) {
        cat('Filtering using a fraction of samples...\n')
        # Filter fraction
        use <- names(sort(rowSums(floor(DATA))))[c((ceiling( ((1-(as.numeric(SETTINGS['filter_fraction', 'setting'])) )*dim(DATA)[1]))):dim(DATA)[1])]
        DATA <- DATA[use,, drop=F]
      } else {
        stop('Filter_fraction settings are incorrect.\n')
      }
    }
    
    # Low feature abundance filter
    if ('filter_low_abundant_features' %in% rownames(SETTINGS)) {
      if (as.numeric(SETTINGS['filter_low_abundant_features', 'setting']) == 0) {
      } else {
        cat('Filtering using a low abundance feature filter...\n')
        max <- apply( DATA / rowSums(DATA) , 2, max)
        ids <- which(max >= as.numeric(SETTINGS['filter_low_abundant_features', 'setting']))
        oldDim <- dim(DATA)[2]
        DATA <- DATA[,ids, drop=F]
        newDim <- dim(DATA)[2]
        system(paste("echo 'Total features: ", oldDim, "\nFeatures kept after low abundance feature filter: ", newDim, "\nFraction features kept: ", signif(newDim/oldDim,3), "\nSetting: ", as.numeric(SETTINGS['filter_low_abundant_features', 'setting']), "\n' > ", WORKING.DIR, "/feature.low.abundance.filter.log", sep=''))
        cat(paste("   Fraction features kept: ", signif(newDim/oldDim,3), "\n", sep=""))
      }
    }
    
    # Prevalence filtering
    if('filter_prevalence_of_features' %in% rownames(SETTINGS)){
      pf <- as.numeric(SETTINGS['filter_prevalence_of_features', 'setting'])
      if (pf == 0) {
      } else if (pf > 0) {
        cat('Filtering using a prevalence feature filter...\n')
        minSamples <- signif(pf*dim(DATA)[1],0)
        ids <- (colSums(DATA > 0) >= minSamples)
        oldDim <- dim(DATA)[2]
        DATA <- DATA[,ids, drop=F]
        newDim <- dim(DATA)[2]
        system(paste("echo 'Total features: ", oldDim, "\nFeatures kept after prevalence feature filter: ", length(ids), "\nFraction features kept: ", signif(newDim/oldDim,3), "\nSetting: ", pf, "\nMin samples feature should be >0 in:", minSamples, "\n' > ", WORKING.DIR, "/feature.prevalance.filter.log", sep=''))
        cat(paste("   Fraction features kept: ", signif(newDim/oldDim,3), "\n", sep=""))
      } else {
        stop("Unknown filter_prevalence_of_features setting!\n")
      }
    }
    
    # Load metadata for only samples left
    if (LOAD.METADATA != 'no') {
      
      save.image(file=paste(DATA.DIR,"/session.RData",sep=""))
      
      cat('Loading metadata...\n')
      METADATA <- METADATA.FULL[rownames(DATA),,drop=F]
      colnames(METADATA) <- colnames(METADATA.FULL)
    }
    
    # Automatically provide DATA split defined by groups
    if (exists("GROUPS")) {
      if (length(GROUPS == 1)) {
        cat ('Group length is 1, creating SAMPLE.GROUPS...\n')
        levels <- levels(as.factor(METADATA[,GROUPS]))
        SAMPLE.GROUPS <- list()
        for (i in 1:length(levels)) {
          SAMPLE.GROUPS[[levels[i]]] <- rownames(METADATA)[METADATA[,GROUPS] == levels[i]]
        }
      }
    }
    
    # prepare iTOL/iPATH functionality
    if ('iTOL' %in% rownames(SETTINGS)) {
      cat('Using iTOL/iPATH data...\n')
      iTOLmultibar <- list()
      iTOLcolorscale <- list()
    }
    
    # Execute the function
    cat('Executing function...\n')
    DATA.BEFORE.EXE <- DATA # This is used for iTOL below
    source(execute)
    FUNCTION.SOURCED <- TRUE
  }
}



# Print results to tables
cat('Saving results...\n')
if (length(attributes(RESULTS)$names) > 0) {
  for (counter in c(1:length(attributes(RESULTS)$names))) {
    filename <- attributes(RESULTS[counter])$names  
    write.csv(RESULTS[[counter]], file=paste(TABLES.DIR,"/",filename, ".table",sep=""), quote=FALSE, row.names=T)
  }
}


if ( exists("iTOLmultibar") || exists("iTOLcolorscale") ) {
  cat('Saving iTOL/iPATH data...\n')  
  included.species <- vector()
  if (exists("iTOLmultibar")) {
    if (length(attributes(iTOLmultibar)$names) > 0) {
      for (counter in c(1:length(attributes(iTOLmultibar)$names))) {
        filename <- attributes(iTOLmultibar[counter])$names  
        write.csv(iTOLmultibar[[counter]], file=paste(DATA.DIR,"/",filename, ".multibar.data.iTOL",sep=""), quote=FALSE, row.names=T)
        included.species <- c(included.species, row.names(iTOLmultibar[[counter]]))
      }
    }
  }
  
  if (exists("iTOLcolorscale")) {
    if (length(attributes(iTOLcolorscale)$names) > 0) {
      for (counter in c(1:length(attributes(iTOLcolorscale)$names))) {
        filename <- attributes(iTOLcolorscale[counter])$names  
        write.csv(iTOLcolorscale[[counter]], file=paste(DATA.DIR,"/",filename, ".colorscale.data.iTOL",sep=""), quote=FALSE, row.names=T)
        included.species <- c(included.species, row.names(iTOLcolorscale[[counter]]))
      }
    }
  }
  
  # Get somewhat abundant taxa
  #load.function('vegan')
  #cs <- colSums(DATA)
  #ncs <- t(decostand(cs, 'total', 2))
  # Gives you top abundant species, using the DATA matrix as it was just before the program may have changed it
  #top.species <- names(sort(rowSums(floor(t(DATA.BEFORE.EXE)))))[c((ceiling( ((1-0.50)*dim(DATA.BEFORE.EXE)[2]))):dim(DATA.BEFORE.EXE)[2])]
  max <- apply( DATA.BEFORE.EXE / rowSums(DATA.BEFORE.EXE), 2, max)
  top.species <- names(which(max >= 0.001))
  all.species <- unique(c(included.species, top.species))
  abundances <- log10(colSums(DATA.BEFORE.EXE+1)[all.species])
  write.csv(abundances, file=paste(DATA.DIR,"/abundances.data.iTOL",sep=""), quote=FALSE, row.names=T)
  
  # Write to file
  write.table(as.matrix(all.species), file=paste(DATA.DIR,"/taxa.iTOL",sep=""), sep="\t", quote=F, row.names=F, col.names=F )
}




# Save R data and exit
cat('Saving R session...\n')
save.image(file=paste(DATA.DIR,"/session.RData",sep=""))
save(DATA, file=paste(DATA.DIR,"/DATA.RData",sep=""))
cat('R SESSION COMPLETED SUCCESFULLY\n')
txtStop()


# Done
