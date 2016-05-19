# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright Jens ROat Kultima, EMBL 2015
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC Runs Friedman test for Anita's day 7 samples
# REQUIRE metadata


########################
# Packages and libraries
########################
#load.function(ggplot2) # This function is MOCAT specfic
#load.function(gridExtra)
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(vsn)
library(miscTools)

###########
# Load data
###########
# IMPORTANT NOTE! From MOCAT, data is by default stored in DATA

FIGURES.DIR <- './'
#FILTER_CUTOFF <- 0.000001
#TRANSFORMATION <- 'log10' # arcsin, DESeq2_log
#FILTER_TYPE <- 'keep with value in at least 1 sample' # 'keep with value in all samples'

FILTER_TYPE <- 'keep_features_above_X_in_Y_samples'
FILTER_MIN_SAMPLES <- 3


TEST <- 'anova'


#for (type2 in c('mOTU')){
for (type2 in c('mOTU', 'module', 'pathway', 'ko', 'cog')){
  
  for (TRANSFORMATION in c('arcsin')) {  #DESeq2_log
    #for (TRANSFORMATION in c('none')) {
    for (FILTER_CUTOFF in c(0.0001)) {
      #for (FILTER_CUTOFF in c( 0.0001)) {
      
      
      
      # Note that this is a fake link to the insert scaled for the mOTU
      file <- paste('/g/bork5/mocat/MM/functional.profiles/263RefGeneCatalog.padded/screened.screened.adapter.on.hg19.solexaqa.allbest.l45.p95/anita.samples.functional.profile.screened.screened.adapter.on.hg19.on.263RefGeneCatalog.padded.solexaqa.allbest.l45.p95.base.norm.', type2, sep='')
      if(TRANSFORMATION == 'DESeq2_vst' & type2 != 'mOTU') {
        file <- paste('/g/bork5/mocat/MM/functional.profiles/263RefGeneCatalog.padded/screened.screened.adapter.on.hg19.solexaqa.allbest.l45.p95/anita.samples.functional.profile.screened.screened.adapter.on.hg19.on.263RefGeneCatalog.padded.solexaqa.allbest.l45.p95.insert.raw.', type2, sep='')
      }
      
      if (type2 == 'mOTU') {
        file <- '/g/bork5/avoigt/TimeSeries/Data/mm.402.motu.is.annotated.mOTU.clusters.gz'
      }
      
      DATA <- t(read.table(file,quote="", header=T, row.names=1, check.names=FALSE, sep='\t', comment='#', skip = 4))
      if (type2 != 'mOTU') {
        
        DATA <- DATA[,c(3:dim(DATA)[2])] # Needed if we have funcitonal data
      }
      
      samples <- read.table('/g/bork1/kultima/projects/anita_functional_analysis/files/samples',quote="", header=F, check.names=FALSE, sep='\t', comment='#') # Load Samples to be included, has to be saved as vector, (no row.names=1 otherwise it expects data in column 1) # Header has to be F otherwise it takes a sample as header
      samples <- as.vector(samples[,1]) #  Make as vector
      DATA_SUBSET <- DATA[rownames(DATA) %in% samples,] # Select wanted samples in table based on vector
      
      
      # This is also loaded from MOCAT
      if(type2 == 'mOTU') {
        type <- 'Species'
      }
      if(type2 == 'cog') {
        type <- 'Orthologous Genes'
      }
      if(type2 == 'ko') {
        type <- 'KEGG Orthologous Genes'
      }
      if(type2 == 'module') {
        type <- 'KEGG Modules'
      }
      if(type2 == 'pathway') {
        type <- 'KEGG Pathways'
      }
      if(type2 == 'gene') {
        type <- 'Genes'
      }
      
      ##########
      # Open PDF
      ##########
      cat (paste('Running : ', type,TRANSFORMATION,FILTER_CUTOFF,TEST,'\n', sep='.'))
      PDF <- paste(type,TRANSFORMATION,FILTER_TYPE,FILTER_MIN_SAMPLES,FILTER_CUTOFF,TEST,'pdf', sep='.')
      #PDF=''
      if (PDF != '') {
        pdf(file=PDF, width=15, height=20)
      }
      
      
      ######################
      # Low abundance filter
      ######################
      
      if (FILTER_TYPE == 'keep with value in at least 1 sample'){
        max <- apply( DATA_SUBSET / rowSums(DATA_SUBSET) , 2, max)
        ids <- which(max >= FILTER_CUTOFF)
        oldDim <- dim(DATA_SUBSET)[2]
        DATA_SUBSET <- DATA_SUBSET[,ids, drop=F]
        newDim <- dim(DATA_SUBSET)[2]
        cat(paste("1. Fraction features kept: ",newDim/oldDim, "\n", sep=""))
      }
      if (FILTER_TYPE == 'keep with value in all samples'){
        stop() # not implemented currently, probably not meaningful
      }
      if (FILTER_TYPE == 'keep_features_above_X_in_Y_samples'){
        oldDim <- dim(DATA_SUBSET)[2]
        DATA_SUBSET <- DATA_SUBSET[,colSums(DATA_SUBSET / rowSums(DATA_SUBSET) >= FILTER_CUTOFF) >= FILTER_MIN_SAMPLES]
        newDim <- dim(DATA_SUBSET)[2]
        cat(paste(" 2. Fraction features kept: ",newDim/oldDim, "\n", sep=""))
      }
      
      
      
      ######################################
      # Rename species with weird annotation
      ######################################
      colnames(DATA_SUBSET) <- gsub(' ','_',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('/','_',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('-','_',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('\\[','',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub(']','',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('\\(','',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub(')','',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('\'','',colnames(DATA_SUBSET))
      colnames(DATA_SUBSET) <- gsub('\\.','',colnames(DATA_SUBSET))
      
      
      ###########
      # Fractions
      ###########
      DATA_PROP <- data.frame(prop.table(as.matrix(DATA_SUBSET), 1)) # Calculates fractions
      DATA_PROP_PT <- DATA_PROP
      
      ################
      # Transformation
      ################
      
      if ( TRANSFORMATION == 'arcsin' ){
        DATA_PROP <- (asin(sqrt(DATA_PROP_PT)))
        meanSdPlot(t(DATA_PROP), main=paste(type, '(arcsin transformation)'))
      }
      if ( TRANSFORMATION == 'log_wolfgang' ){
        lgf <- function(x){
          (x+sqrt(x*x+1))/2
        }
        DATA_PROP <- lgf(DATA_PROP_PT)
      }
      if ( TRANSFORMATION == 'log10' ){
        DATA_PROP <- log10(DATA_PROP_PT+1e-6)
        meanSdPlot(t(DATA_PROP), main=paste(type, '(log10 transformation)'))
      }
      if ( TRANSFORMATION == 'DESeq2_vst' ){
        
        DATA_PROP2 <- NULL
        DATA_PROP2$Patient <- sapply(rownames(DATA_PROP),function(x) {return(strsplit(x,'-')[[1]][1])})
        DATA_PROP2$Protocol <- sapply(rownames(DATA_PROP),function(x) {return(strsplit(x,'-')[[1]][2])})
        DATA_PROP2$Day <- sapply(rownames(DATA_PROP),function(x) {return(strsplit(x,'-')[[1]][3])})
        DATA_PROP2$Replicate <- sapply(rownames(DATA_PROP),function(x) {return(strsplit(x,'-')[[1]][4])})
        DATA_PROP2 <- as.data.frame(DATA_PROP2)
        
        DATA_PROP3 <- list()
        DATA_PROP4 <- NULL
        
        # ALL
        #         countD <- as.matrix(t(round(DATA_SUBSET)))
        #         DESEQ_DATA <- DESeqDataSetFromMatrix(countData = countD,
        #                                              colData = DATA_PROP2,
        #                                              design = ~ 1)
        #         DESEQ_DATA_VST <- varianceStabilizingTransformation(DESEQ_DATA, blind=TRUE)
        #         DATA_PROP <- t(assay(DESEQ_DATA_VST))
        #         plot(main=paste(patient, protocol, day) ,  ( as.numeric(as.matrix(DATA_PROP[2,])) )  , as.matrix((round(DATA_SUBSET[2,])))    )
        
        
        #         for (patient in as.vector(unique(DATA_PROP2$Patient))) {
        #           for (protocol in c('11', '12', '13') ) {
        for (day in c('7', '392') ) {
          X <- paste(day, sep='-')
          cD <- (DATA_PROP2[grepl(X, rownames(DATA_SUBSET)), ])
          #cD$Patient <- factor(, levels=X)
          
          countD <- as.matrix(t(round(DATA_SUBSET[ grepl(X, rownames(DATA_SUBSET)), ,drop=F])), drop=F)
          DESEQ_DATA <- DESeqDataSetFromMatrix(countData = countD,
                                               colData = cD,
                                               design = ~ 1)
          DESEQ_DATA_VST <- varianceStabilizingTransformation(DESEQ_DATA, blind=TRUE)
          #DESEQ_DATA_LOG <- rlogTransformation(DESEQ_DATA)
          DATA_PROP3[[X]] <- t(assay(DESEQ_DATA_VST))
          DATA_PROP4 <- rbind(DATA_PROP4, t(assay(DESEQ_DATA_VST)) )
          #plot(main=paste(patient, protocol, day) ,  ( as.numeric(as.matrix(DATA_PROP3[[patient]])) )  , as.matrix((round(DATA_SUBSET[grepl(X, rownames(DATA_SUBSET)), ])))    )
          meanSdPlot(t( DATA_PROP3[[X]] ), main=X)
        }
        #           }
        DATA_PROP <- DATA_PROP4
        meanSdPlot(t(DATA_PROP), main=paste(type, '(VST)'))
        #      cds <- estimateSizeFactors( DESEQ_DATA )
        #     cds <- estimateDispersions( cds )
        #      plotDispEsts(cds)
      }
      
      
      
      
      
      
      
      ################
      # Loop over days
      ################
      FOR_PLOTTING <- list()
      LINE <- list()
      for (DAY in c(7, 392)) {
        subset_names <- names(which(sapply(row.names(DATA_SUBSET), function(x) {return(strsplit(x,'-')[[1]][3])})==DAY)) # extract names for this day
        SUBSET <- data.frame(DATA_PROP[subset_names,]) # get data subset
        SUBSET_PT <- data.frame(DATA_PROP_PT[subset_names,]) # get data subset
        
        ########################################
        # Take average, of replicates if day 392
        ########################################
        if (DAY == 7) {
          SUBSET_7 <- SUBSET
          SUBSET_7_PT <- SUBSET_PT
        }
        if (DAY == 392) {
          patient_protocol <- unique(paste(sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][1])}), sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][2])}), sep='-'))
          SUBSET_MEAN <- NULL
          new_rownames <- NULL
          for (pp in patient_protocol) {
            SUBSET_MEAN <- rbind(SUBSET_MEAN, t(colMedians(SUBSET[grepl(pp, rownames(SUBSET)),,drop=F])))
            new_rownames <- c(new_rownames, grep(pp, rownames(SUBSET), value = T)[1])
          }
          rownames(SUBSET_MEAN) <- new_rownames
          SUBSET <- data.frame(SUBSET_MEAN)
          SUBSET_392 <- SUBSET
          
          SUBSET_MEAN_PT <- NULL
          new_rownames <- NULL
          for (pp in patient_protocol) {
            SUBSET_MEAN_PT <- rbind(SUBSET_MEAN_PT, t(colMedians(SUBSET_PT[grepl(pp, rownames(SUBSET_PT)),,drop=F])))
            new_rownames <- c(new_rownames, grep(pp, rownames(SUBSET_PT), value = T)[1])
          }
          rownames(SUBSET_MEAN_PT) <- new_rownames
          SUBSET_PT <- data.frame(SUBSET_MEAN_PT)
          SUBSET_392_PT <- SUBSET_PT          
        }
        
        
        #####################################
        # Identifiy Protocol, Patient and Day
        #####################################
        SUBSET$Patient <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][1])})
        SUBSET$Protocol <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][2])})
        SUBSET$Day <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][3])})
        SUBSET$Replicate <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][4])})
        
        
        ######################
        # Create empty vectors
        ######################
        p_values <- vector(length = ncol(SUBSET)-4) # create empty vectors to store data
        q_ratios <- vector(length = ncol(SUBSET)-4) # create empty vectors to store data
        
        
        ####################
        # Loop over features
        ####################
        for (c in 1:(ncol(SUBSET)-4)) { # does not consider the last 4 columns that include the patient, protocol, day
          
          if (TEST == 'anova') {
            r <- aov(formula( paste(colnames(SUBSET)[c], '~ Patient + Protocol' )), SUBSET)
            r <- unlist(summary(r))
            p_values[c] <- r[14]
            q_ratios[c] <- r[5]/r[4]
          }
          
          if (TEST == 'friedman') {
            r <- friedman.test(formula( paste(colnames(SUBSET)[c], '~ Protocol | Patient' )), SUBSET) # run test
            p_values[c] <- r$p.value # save values
            #if (r$statistic == 0) {
            q_ratios[c] <- r$statistic # save values
            #} else {
            #  q_ratios[c] <- 1/r$statistic # save values
            #}
          } 
        }
        
        #################
        # Adjust p values
        #################
        p_adjusted2 <- p.adjust(p_values, method = "BH") # adjust p-values
        p_adjusted <- p_values
        
        
        # Calculate the vertical line, one per day
        ret <- -1
        seq <- 1
        while (ret == -1) {
          seq <- seq * 0.99
          pvalue <- c(seq, p_values)
          pvalueA <- p.adjust(pvalue, method = "BH")
          if (pvalueA[1] <= 0.05) {
            ret <- seq
          }
        }
        LINE[[DAY]] <- ret
        
        ################################
        # Create data frame for plotting
        ################################
        df <- as.data.frame(p_adjusted) 
        df$q_ratios <- q_ratios
        df$p_adjusted2 <- p_adjusted2
        rownames(df) <- paste(colnames(SUBSET), ' (', DAY, ')', sep='')[1:(dim(df)[1])]
        FOR_PLOTTING[[DAY]] <- df
        
      }
      
      
      
      
      
      
      
      
      ##################
      # Merge dataframes
      ##################
      df <- rbind(cbind(FOR_PLOTTING[[7]], 7))
      colnames(df)[4] <- 392            
      df <- rbind(df, cbind(FOR_PLOTTING[[392]], 392))
      colnames(df)[4] <- 'day'
      df <- df[!is.na(df$p_adjusted2),]
      #df[!(df$p_adjusted2 <= 0.05 & df$q_ratios > 1),'day'] <- 'C'
      df$day <- factor(df$day, levels=c('7', '392', 'C'))
      df$alpha <- 1
      df[!(df$p_adjusted2 <= 0.05 & df$q_ratios > 1),'alpha'] <- 0.3
      df$fill <- df$day       
      df[(df$p_adjusted2 <= 0.05 & df$q_ratios > 1),'fill'] <- 'C'
      df$size <- 4
      df[(df$p_adjusted2 <= 0.05 & df$q_ratios > 1),'size'] <- 6
      df$shape <- 17
      #df[(df$p_adjusted2 <= 0.05 & df$q_ratios > 1),'shape'] <- 17
      df[df$day==392,'shape'] <- 16
      
      ######
      # plot
      ######
      
      
      COLORS <- c('#546bab', '#87889c', '#ED1C24', '#44B7E9', '#A5CF51', '#FBB040', '#6FC59B', '#A97C50', '#231F20', '#BCBEC0')
      names(COLORS) <- c('7','392', 'alien', 'bugkiller', 'peacemaker', 'scavenger', 'halbarad', '11', '12', '13')
      
      
      #########################
      # Plot dots and densities
      #########################
      
      # placeholder plot - prints nothing at all
      empty <- ggplot()+geom_point(aes(1,1), colour="white") +
        theme(                              
          plot.background = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
        )
      
      # Right hand density
      plot_right <- ggplot() + geom_histogram(aes(df$q_ratios, fill=factor(df$day) ), alpha=0.5) + coord_flip() + theme(legend.position = "none")  +
        theme(                              
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          plot.margin=unit(c(0,0,0,-0.5), "cm")
        ) + scale_fill_manual(values = COLORS) + scale_x_log10()
      #+ xlim(min(df[,2]),max(df[,2]))
      
      # Top density
      plot_top <- ggplot() + geom_histogram(aes(df$p_adjusted, fill=factor(df$day) ), alpha=0.5) + theme(legend.position = "none")  +
        theme(                              
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          plot.margin=unit(c(0,0,-0.5,0), "cm")
        ) + scale_fill_manual(values = COLORS) + scale_x_log10()
      #+ xlim(0,1)
      
      # Main
      d7UL <- signif(sum(df$p_adjusted2 <= 0.05 & df$q_ratios > 1 & df$day == 7, na.rm = T)/sum(!is.na(FOR_PLOTTING[[7]]$p_adjusted))*100,2)
      d7BL <- signif(sum(df$p_adjusted2 < 0.05 & df$q_ratios <= 1 & df$day == 7, na.rm = T)/sum(!is.na(FOR_PLOTTING[[7]]$p_adjusted))*100,2)
      d7UR <- signif(sum(df$p_adjusted2 > 0.05 & df$q_ratios > 1 & df$day == 7, na.rm = T)/sum(!is.na(FOR_PLOTTING[[7]]$p_adjusted))*100,2)
      d7BR <- signif(sum(df$p_adjusted2 > 0.05 & df$q_ratios <= 1 & df$day == 7, na.rm = T)/sum(!is.na(FOR_PLOTTING[[7]]$p_adjusted))*100,2)
      
      d392UL <- signif(sum(df$p_adjusted2 <= 0.05 & df$q_ratios > 1 & df$day == 392, na.rm = T)/sum(!is.na(FOR_PLOTTING[[392]]$p_adjusted))*100,2)
      d392BL <- signif(sum(df$p_adjusted2 < 0.05 & df$q_ratios <= 1 & df$day == 392, na.rm = T)/sum(!is.na(FOR_PLOTTING[[392]]$p_adjusted))*100,2)
      d392UR <- signif(sum(df$p_adjusted2 > 0.05 & df$q_ratios > 1 & df$day == 392, na.rm = T)/sum(!is.na(FOR_PLOTTING[[392]]$p_adjusted))*100,2)
      d392BR <- signif(sum(df$p_adjusted2 > 0.05 & df$q_ratios <= 1 & df$day == 392, na.rm = T)/sum(!is.na(FOR_PLOTTING[[392]]$p_adjusted))*100,2)
      
      if(is.na(d7BL)){d7BL<-0}
      if(is.na(d7BR)){d7BR<-0}
      if(is.na(d7UL)){d7UL<-0}
      if(is.na(d7UR)){d7UR<-0}
      
      if(is.na(d392BL)){d392BL<-0}
      if(is.na(d392BR)){d392BR<-0}
      if(is.na(d392UL)){d392UL<-0}
      if(is.na(d392UR)){d392UR<-0}
      
      # Get min and max x and scale
      roundDown <- function(x) 10^floor(log10(x))
      roundUp <- function(x) 10^ceiling(log10(x))
      i <- log10(roundDown(min(df$p_adjusted, na.rm = T)))
      j <- log10(roundUp(max(df$p_adjusted, na.rm = T)))
      x_steps <- 10^(seq(from = i-1, to = j, by = 1))
      i <- log10(roundDown(min(df$q_ratios, na.rm = T)))
      j <- log10(roundUp(max(df$q_ratios, na.rm = T)))
      y_steps <- 10^(seq(from = i, to = j, by = 1))
      
      #x_steps <- 10^(seq(from = -6, to = 0, by = 1))
      #y_steps <- 10^(seq(from = -6, to = 0, by = 1))
      
      
      # make fancy converter
      fancy_scientific <- function(l) {
        l <- format(l, scientific = TRUE)
        l <- gsub("^(.*)e", "e", l)
        l <- gsub("e", "10^", l)
        l <- gsub("\\+", "", l)
        parse(text=l)
      }
      
      main <- ggplot(df, aes(p_adjusted,q_ratios),label=df$name,) + geom_point( aes(color=df$day,alpha=0.9), size=6, shape=df$shape)  + #df$day
        ylab("Ratio of variance explained: Method effects over subject effects\n") +
        xlab("\nFDR adjusted p-values") +
        geom_hline(yintercept=1, colour='black', size = 1, linetype="dashed") +
        #geom_vline(xintercept=0.05, colour='red', size = 1, linetype="dotted") +
        geom_vline(xintercept=LINE[[7]], colour=COLORS['7'], size = 1, linetype="dashed") +
        geom_vline(xintercept=LINE[[392]], colour=COLORS['392'], size = 1, linetype="dashed") +
        theme_bw() +
        theme(
          # plot.margin=unit(c(0,-0.5,0,0), "cm"), only add this if adding histogram
          legend.position="none",
          axis.text=element_text(size=22),
          axis.title=element_text(size=22),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()       
        ) +
        scale_fill_manual(values = COLORS) + scale_colour_manual(values = COLORS) +
        #annotate("text", label=paste('N=',dim(FOR_PLOTTING[[7]])[1], sep='' ), x = max(LINE[[7]], LINE[[392]], na.rm = T), y = 2 * min(df$q_ratios,na.rm = T), size = 8, colour = COLORS['7'], hjust = 0) +
        #annotate("text", label=paste('N=',dim(FOR_PLOTTING[[392]])[1], sep='' ), x = max(LINE[[7]], LINE[[392]], na.rm = T), y = min(df$q_ratios,na.rm = T), size = 8, colour = COLORS['392'], hjust = 0) +
        annotate("text", label=paste('N=',sum(!is.na(FOR_PLOTTING[[7]]$p_adjusted)), sep='' ), x = x_steps[1], y = 2 * min(df$q_ratios,na.rm = T), size = 8, colour = COLORS['7'], hjust = 0) +
        annotate("text", label=paste('N=',sum(!is.na(FOR_PLOTTING[[392]]$p_adjusted)), sep='' ), x = x_steps[1], y = min(df$q_ratios,na.rm = T), size = 8, colour = COLORS['392'], hjust = 0) +
        
        #         annotate("text", label=paste(d7UL, '%  ', sep=''), x = LINE[[7]], y = 1.9, size = 8, colour = COLORS['7'], hjust = 1) +
        #         annotate("text", label=paste(d392UL, '%  ', sep=''), x = LINE[[392]], y = 4, size = 8, colour = COLORS['392'], hjust = 1) +
        #         annotate("text", label=paste("   ", d7UR, '%', sep=''), x = LINE[[7]], y = 1.9, size = 8, colour = COLORS['7'], hjust = 0) +
        #         annotate("text", label=paste("   ", d392UR, '%', sep=''), x = LINE[[392]], y = 4, size = 8, colour = COLORS['392'], hjust = 0) +
        annotate("text", label=paste(d7UL, '%  ', sep=''), x = LINE[[7]], y = max(df$q_ratios), size = 8, colour = COLORS['7'], hjust = 1) +
        annotate("text", label=paste(d392UL, '%  ', sep=''), x = LINE[[392]], y = max(df$q_ratios), size = 8, colour = COLORS['392'], hjust = 1) +
        annotate("text", label=paste("   ", d7UR, '%', sep=''), x = LINE[[7]], y = max(df$q_ratios), size = 8, colour = COLORS['7'], hjust = 0) +
        annotate("text", label=paste("   ", d392UR, '%', sep=''), x = LINE[[392]], y = max(df$q_ratios), size = 8, colour = COLORS['392'], hjust = 0) +
        
        #annotate("text", label=paste(d7BL, '%  ', sep=''), x = LINE[[7]], y = 0.55, size = 8, colour = COLORS['7'], hjust = 1) +
        #annotate("text", label=paste(d392BL, '%  ', sep=''), x = LINE[[392]], y = 0.25, size = 8, colour = COLORS['392'], hjust = 1) +
        #annotate("text", label=paste("  ", d7BR, '%', sep=''), x = LINE[[7]], y = 0.55, size = 8, colour = COLORS['7'], hjust = 0) +
        #annotate("text", label=paste("  ", d392BR, '%', sep=''), x = LINE[[392]], y = 0.25, size = 8, colour = COLORS['392'], hjust = 0) +    
        # xlim(0,1) + ylim(min(df[,2]),max(df[,2])) + Reminant of when we didn't use log scale and had histograms
        scale_x_continuous(breaks=c(x_steps), trans="log10", labels=fancy_scientific, limits=c(min(x_steps), max(x_steps))) +
        scale_y_continuous(breaks=c(y_steps), trans="log10", labels=fancy_scientific) +
        #theme(panel.border = element_rect(colour = "black"))
        scale_shape(solid = TRUE)
      
      
      
      # Fix scales
      gB <- ggplotGrob(main)
      gA <- ggplotGrob(plot_top)
      gC <- ggplotGrob(empty)
      gD <- ggplotGrob(plot_right)
      maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
      maxHeight = grid::unit.pmax(gB$heights[2:5], gD$heights[2:5])
      gA$widths[2:5] <- as.list(maxWidth)
      gB$widths[2:5] <- as.list(maxWidth)
      gB$heights[2:5] <- as.list(maxHeight)
      gD$heights[2:5] <- as.list(maxHeight)
      
      # Plot
      
      
      
      
      
      
      
      #########################
      # Plot correlation graphs
      #########################
      
      # Construct data
      COMBINED <- rbind(SUBSET_7_PT, SUBSET_392_PT)
      patients <- unique(sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][1])}))
      protocols <- sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][2])})
      days <- sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][3])})
      both <- paste(protocols, days, sep='-')
      both_u <- unique(both)
      both_u7 <-grep('-7', both_u, value=T)
      both_u392 <-grep('-392', both_u, value=T)
      
      # Loop over protocols
      COR <- NULL
      for (d in c(7,392)){ #for (d in c(7, 392)){
        if(d == 7) {
          both_U <- both_u7
        } else {
          both_U <- both_u392
        }
        for (i in c(1:length(both_U))) {
          for (j in c(i+1:length(both_U))){
            for (patient in patients) {
              x <- COMBINED[grepl(paste(patient,both_U[i], sep='-'), rownames(COMBINED)),]
              y <- COMBINED[grepl(paste(patient,both_U[j], sep='-'), rownames(COMBINED)),]
              p1 <- strsplit(both_U[i],'-')[[1]][1]
              p2 <- strsplit(both_U[j],'-')[[1]][1]
              if(dim(x)[1] > 0 & dim(y)[1] > 0) {
                
                if(p1 == '11' & p2 == '12') {
                  V <- 'A'
                }
                if(p1 == '11' & p2 == '13' & d == 7) {
                  V <- 'B'
                }
                if(p1 == '11' & p2 == '13' & d == 392) {
                  V <- 'D'
                }
                if(p1 == '12' & p2 == '13') {
                  V <- 'C'
                }
                COR <- rbind(COR, c(V, d, patient, cor(  (as.numeric(x)),  (as.numeric(y)),  method = "spearman"), NA))
              }
            }
          }
        }
      }
      
      # Loop over days
      for (i in c(11,13)) {
        for (patient in patients) {
          x <- COMBINED[grepl(paste(patient,i,7, sep='-'), rownames(COMBINED)),]
          y <- COMBINED[grepl(paste(patient,i,392, sep='-'), rownames(COMBINED)),]
          if(dim(x)[1] > 0 & dim(y)[1] > 0) {
            COR <- rbind(COR, c(paste('Day 7 vs Day 392'), '5', patient, cor(  (as.numeric(x)),  (as.numeric(y)),  method = "spearman") , NA))
          }
        }
      }
      
      # Loop over individuals
      for (protocol in unique(protocols)) {
        for (day in unique(days)) {
          s <- COMBINED[grepl(paste("",protocol, day, "", sep='-'), rownames(COMBINED)),]
          for (i in c(1:dim(s)[1])) {
            for (j in c(i+1:dim(s)[1])){
              COR <- rbind(COR,   c(       'Subject vs\nSubject',        day,          protocol,        cor((as.numeric(s[i,])) ,(as.numeric(s[j,])) ,method = "spearman"),       paste(rownames(s)[i], 'vs', rownames(s)[j]) ))
            }
          }
        }
      }
      
      COR <- COR[!is.na(COR[,4]),] 
      COR <- data.frame(COR)
      colnames(COR) <- c('protocol', 'day', 'patient', 'value', 'Ind')
      COR$value <- as.numeric(as.vector(COR$value))
      COR$f <- factor( rep( paste('', type, "", sep='' ) , dim(COR),dim(COR)[1]))
      COR$day <- factor(COR$day, levels=c('7', '392', '5', '66'))
      COR$shape <- 17
      COR[COR$day == 392,'shape'] <- 16
      COR[COR$day == 5,'shape'] <- 18
      COR$size <- 6
      COR[COR$day == 5,'size'] <- 8
      
      # Create funciton for median line
      f <- function(x, height) {
        ans <- median(x)
        data.frame(ymin = ans-height/2, ymax = ans+height/2, y = ans)
      }
      
      # save plot
      plot2 <- ggplot(COR, aes(factor(protocol), as.numeric(value))) + geom_jitter(size=COR$size, shape=COR$shape, position = position_jitter(width = 0.2), aes(color=(COR$patient))) + 
        stat_summary(fun.data = f, geom = "crossbar", height = 0.005,
                     colour = NA, fill = "black", width = 0.5, alpha = 1)+facet_grid(.~f) + theme_bw() +
        theme(legend.position="none") + 
        theme(axis.text.x=element_text(size=16),
              axis.text.y=element_text(size=22),
              axis.title=element_text(size=22),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()       
        ) +
        theme(strip.text.x = element_text(size = 28, angle = 0)) + #, face="bold"
        xlab('') + ylab('Rank correlation within subjects\n') +
        scale_fill_manual(values = COLORS) + scale_colour_manual(values = COLORS) +
        scale_x_discrete(labels = c(
          expression(atop("Frozen vs RNALater ","(RT, 1w)")),
          expression(atop("Frozen vs RNALater",("+4-10"~degree*C~", 1w"))),
          expression(atop("RNALater (RT, 1w)","vs "("+4-10"~degree*C~", 1w"))),
          expression(atop("Frozen vs RNALater ",("+4"~degree*C~", 24h"))),
          expression(atop("Day 7 vs","Day 392")),
          expression(atop("Subject vs","Subject"))
        )) + ylim(0, 1.01)
      
      
      ###################
      # Finally plot both
      ###################
      # with histograms, doesn't work properly yet with log scales
      #grid.arrange(plot2,
      #             arrangeGrob(gA, gC,gB, gD, ncol=2, nrow=2, widths=c(6, 1), heights=c(1, 6)),
      #             nrow=2, heights=c(0.5,1)
      grid.arrange(plot2, main, nrow=2, heights=c(0.5,1))
      if (PDF != '') {
        dev.off()
      }
      
      SVG <- paste(type,TRANSFORMATION,FILTER_TYPE,FILTER_MIN_SAMPLES,FILTER_CUTOFF,TEST,'svg', sep='.')
      svg(filename = SVG, width = 15, height=20)
      grid.arrange(plot2, main, nrow=2, heights=c(0.5,1))
      dev.off()
      
      ##################
      # Print some stats
      ##################
      #pr <- df[df$p_adjusted2 <= 0.05,]
      pr <- df
      write.csv(pr, file=paste(type,TRANSFORMATION,FILTER_CUTOFF,FILTER_TYPE,FILTER_MIN_SAMPLES,TEST,'txt', sep='.'))
      save.image(file = paste(type,TRANSFORMATION,FILTER_CUTOFF,FILTER_TYPE,FILTER_MIN_SAMPLES,TEST,'RData', sep='.'))
    }
  }
}

