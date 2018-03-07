# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC calculates the wilcoxon statistic for each feature between to states, eg cancer and non cancer
# PRODUCES xls
# PRODUCES pdf
# REQUIRE metadata
# REQUIRE groups|1
# REQUIRE condA
# REQUIRE condB
# REQUIRE iTOL
# SETTINGS : additional_analyses [none,plot_abundance_sum] {CHECK} {none} (Select additional analyses to perform, plot_abundance_sum generates Fig 3 in Qin et al 2012 on T2D)

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
# load('X') - This loads any required packages. Please use this function and not require.

# How to save results
# PREFERABLY SAVE THE RESULT IN THE file.result vector, which should have one entry per sample
# THEN THE GENERIC PART BELOW will take care of saving into the results list and everything will
# work out just great!

# load packages
load.function('ggplot2')
load.function('reshape')
load.function('wq')

#### PROGRAM SPECIFIC ####
# This is the second part of the wilcox_heatmap pipeline
# It's time ot make some plots


cat ('Loading data...\n')
args <- commandArgs(trailingOnly = TRUE)
DATAlistH <-DATAlist
SETTINGS2 <- SETTINGS
FILE <- SETTINGS['load_file','setting']
load(FILE)
SETTINGS1 <- SETTINGS
SETTINGS <- SETTINGS2
DATA.DIR <- SETTINGS['data_dir','setting']
WORKING.DIR <- SETTINGS['working_dir','setting']
TABLES.DIR <- SETTINGS['tables_dir','setting']
FIGURES.DIR <- SETTINGS['figures_dir','setting']
for (name in names(DATAlistH)) {
  DATAlist[[name]] <- DATAlistH[[name]]
}
rm(DATAlistH)

# Let's get all the levels we want to plot for
LEVELS <- gsub('.horizontal', '', names(DATAlist)[grep('.horizontal', names(DATAlist))])
g1 <- SAMPLE.GROUPS[[1]]
g2 <- SAMPLE.GROUPS[[2]]
g1n <- names(SAMPLE.GROUPS)[[1]]
g2n <- names(SAMPLE.GROUPS)[[2]]

for (level in LEVELS) {
  if (level != 'gene') {
    
    get <- paste('significant_', level, sep='')
    file <- SETTINGS[get,'setting']
    significant <- read.table(file, header=T, quote="", row.names=1, sep=",", check.names=F, na.strings="NA", comment.char="")
    selected <- significant[significant$p.adj.fdr <= 0.05,]
    
    cat(paste('Processing', level, '\n'))
    horizontal <- paste(level,'.horizontal', sep='')
    d <- DATAlist[[level]] # get data
    d <- d[,rownames(selected)]
    
    # these lines splits d into two parts and sorts them according to total abundances across samples
    #     d1 <- d[g1,sort(colSums(d), decreasing = F, index.return=T)$ix] # sort it
    #     d1 <- d1[sort(rowSums(d1), dec=T, index.return=T)$ix,] # sort it again
    #     d2 <- d[g2,sort(colSums(d), decreasing = F, index.return=T)$ix] # sort it
    #     d2 <- d2[sort(rowSums(d2), dec=T, index.return=T)$ix,] # sort it again
    
    # These lines keeps the order of significance, with most significant on top, but reorders the samples
    d1 <- d[g1,]
    d1 <- d1[sort(rowSums(d1), dec=T, index.return=T)$ix,] # sort it again
    d2 <- d[g2,]
    d2 <- d2[sort(rowSums(d2), dec=T, index.return=T)$ix,] # sort it again
    
    
    
    d <- rbind(d1, 0, d2)
    dm <- melt(d)
    dm$X1 <- factor(dm$X1, levels=rownames(d))
    dm$X2 <- factor(dm$X2, levels=colnames(d))
    
    h <- DATAlist[[horizontal]] # get data
    h1 <- h[rownames(d1),colnames(d1)]
    h2 <- h[rownames(d2),colnames(d2)]
    h <- rbind(h1, 0, h2)
    hm <- melt(h)
    hm$X1 <- factor(dm$X1, levels=rownames(d))
    hm$X2 <- factor(dm$X2, levels=colnames(d))
    
    dm$X1 <- paste((dm$X1), '-C', sep='')
    #dm$X1 <- factor(dm$X1, levels=paste(rep(rownames(d), each=2), c('-C', '-H'), sep=''))
    hm$X1 <- paste((hm$X1), '-H', sep='')
    #hm$X1 <- factor(hm$X1, levels=paste(rep(rownames(d), each=2), c('-C', '-H'), sep=''))
    
    dm <- cbind(dm, 'c')
    colnames(dm) <- c('X1', 'X2', 'value', 'type')
    hm <- cbind(hm, 'h')
    colnames(hm) <- c('X1', 'X2', 'value', 'type')
    cm <- rbind(dm, hm)
    
    # rescale
    cm[cm$type == 'c','value'] <- cm[cm$type == 'c','value'] * ( max(cm[cm$type == 'h','value']) / max(cm[cm$type == 'c','value']) )
    cm[cm$type == 'c','value'] <- -cm[cm$type == 'c','value']
    
    cm$X1 <- factor(cm$X1, levels=paste(rep(rownames(d), each=2), c('-C', '-H'), sep=''))
    cm_levels <- paste(rep(rownames(d), each=2), c('-C', '-H'), sep='')
    C <- ggplot(cm, aes(X1, X2)) +
      geom_tile(aes(fill = value), colour = "white") +
      scale_fill_gradient2(low ="darkorange4", mid = "white", high = "steelblue4") +
      scale_y_discrete(limits=unique(cm$X2)) +
      scale_x_discrete(limits=cm_levels) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, color="black")) +
      #theme(legend.justification=c(0,1), legend.position=c(0,0.90)) +
      theme(legend.justification=c(0,1), legend.position="none") +
      xlab(paste(g1n, ' - ', g2n, ' | orange=coverage, blue=hosrizontal coverage', sep='')) +
      theme(plot.title = element_text(size = rel(3))) +
      theme(axis.text.x = element_text(size = rel(3))) + 
      theme(axis.title.x = element_text(size = rel(3))) +
      theme(axis.text.y = element_text(size = rel(3))) + 
      theme(axis.title.y = element_text(size = rel(3))) +
      theme(legend.text = element_text(size = rel(0))) +
      theme(strip.text = element_text(size = rel(0))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title.y=element_blank())
    
    # These plots are still good to have  
    D <- ggplot(dm, aes(X1, X2)) +
      geom_tile(aes(fill = value), colour = "white") +
      scale_fill_gradient(low = "white", high = "darkorange4") +
      scale_y_discrete(limits=unique(dm$X2)) +
      scale_x_discrete(limits=unique(dm$X1)) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, color="black")) +
      #theme(legend.justification=c(0,1), legend.position=c(0,0.90)) +
      theme(legend.justification=c(0,1), legend.position="none") +
      xlab(paste(g1n, ' - ', g2n, ' | orange=coverage', sep='')) +
      theme(plot.title = element_text(size = rel(3))) +
      theme(axis.text.x = element_text(size = rel(3))) + 
      theme(axis.title.x = element_text(size = rel(3))) +
      theme(axis.text.y = element_text(size = rel(3))) + 
      theme(axis.title.y = element_text(size = rel(3))) +
      theme(legend.text = element_text(size = rel(0))) +
      theme(strip.text = element_text(size = rel(0))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title.y=element_blank())
    
    H <- ggplot(hm, aes(X1, X2)) +
      geom_tile(aes(fill = value), colour = "white") +
      scale_fill_gradient(low = "white", high = "steelblue4") +
      scale_y_discrete(limits=unique(hm$X2)) +
      scale_x_discrete(limits=unique(hm$X1)) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, color="black")) +
      #theme(legend.justification=c(0,1), legend.position=c(0,0.90)) +
      theme(legend.justification=c(0,1), legend.position="none") +
      xlab(paste(g1n, ' - ', g2n, ' | blue=horizontal coverage', sep='')) +
      theme(plot.title = element_text(size = rel(3))) +
      theme(axis.text.x = element_text(size = rel(3))) + 
      theme(axis.title.x = element_text(size = rel(3))) +
      theme(axis.text.y = element_text(size = rel(3))) + 
      theme(axis.title.y = element_text(size = rel(3))) +
      theme(legend.text = element_text(size = rel(0))) +
      theme(strip.text = element_text(size = rel(0))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title.y=element_blank())
    
    cat(paste('Plotting', level, '\n'))
    jpeg(filename=paste(FIGURES.DIR, '/', 'figure.',level,'.jpeg', sep=''), width = dim(d)[1]*250, height=(4000+dim(selected)[1]*100))
    multiplot(C, D, H, cols=1)
    dev.off()
    
    
    ################### experimental #######################
    # cancer modules
    #   cancer_modules <- c('M00080', 'M00063', 'M00060', 'M00320', 'M00064')
    #   exp <- cm[cm$X2 %in% cancer_modules,]
    #   E <- ggplot(exp, aes(X1, X2)) +
    #     geom_tile(aes(fill = value), colour = "white") +
    #     scale_fill_gradient2(low ="darkorange4", mid = "white", high = "steelblue4") +
    #     scale_y_discrete(limits=unique(exp$X2)) +
    #     scale_x_discrete(limits=cm_levels) +
    #     theme(axis.text.x = element_text(angle = 25, hjust = 1, color="black")) +
    #     #theme(legend.justification=c(0,1), legend.position=c(0,0.90)) +
    #     theme(legend.justification=c(0,1), legend.position="none") +
    #     xlab(paste(g1n, ' - ', g2n, ' | orange=coverage, blue=hosrizontal coverage', sep='')) +
    #     theme(plot.title = element_text(size = rel(3))) +
    #     theme(axis.text = element_text(size = rel(1))) + 
    #     theme(axis.title = element_text(size = rel(1))) +
    #     theme(legend.text = element_text(size = rel(0))) +
    #     theme(strip.text = element_text(size = rel(0))) +
    #     theme(axis.title.x = element_text(size = rel(1))) +
    #     theme(axis.title.y = element_text(size = rel(0.5))) +
    #     theme(axis.text.y = element_text(size = rel(0.5))) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #     theme(axis.title.y=element_blank())
    ########################################################
    #   layOut(
    #     list(C, 1, 1:2),
    #     list(D, 2, 1),
    #     list(H, 2, 2)
    #   )
  } # end if not gene
} # end levels









