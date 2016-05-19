# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright Jens ROat Kultima, EMBL 2015
# This code is released under GNU GPL v3.

# THESE FIELDS DETERMINE REQUIRED INPUT AND DESCRIPTION
# DESC Runs Friedman test for Anita's day 7 samples
# REQUIRE metadata


########################
# Packages and libraries
########################
load.function(ggplot2) # This function is MOCAT specfic
load.function(gridExtra)
load.function(gtable)

###########
# Load data
###########
# IMPORTANT NOTE! From MOCAT, data is by default stored in DATA
# DATA <- t(read.table('/g/bork5/avoigt/TimeSeries/Data/mm.402.motu.is.annotated.mOTU.clusters.gz',quote="", header=T, row.names=1, check.names=FALSE, sep='\t', comment='#'))
samples <- read.table('/g/bork1/kultima/projects/anita_functional_analysis/files/samples',quote="", header=F, check.names=FALSE, sep='\t', comment='#') # Load Samples to be included, has to be saved as vector, (no row.names=1 otherwise it expects data in column 1) # Header has to be F otherwise it takes a sample as header
samples <- as.vector(samples[,1]) #  Make as vector
DATA_SUBSET <- DATA[rownames(DATA) %in% samples,] # Select wanted samples in table based on vector

# This is also loaded from MOCAT
#TYPE <- "base.norm.pathway.solexaqa.allbest.l45.p95"
type <- strsplit(TYPE, '\\.')[[1]][3]
if(type == 'mOTU') {
  type <- 'Species'
}
if(type == 'cog') {
  type <- 'COGs'
}
if(type == 'ko') {
  type <- 'KOs'
}
if(type == 'module') {
  type <- 'Modules'
}
if(type == 'pathway') {
  type <- 'Pathways'
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




################
# Loop over days
################
FOR_PLOTTING <- list()
for (DAY in c(7, 392)) {
  subset_names <- names(which(sapply(row.names(DATA_SUBSET), function(x) {return(strsplit(x,'-')[[1]][3])})==DAY)) # extract names for this day
  SUBSET <- data.frame(DATA_PROP[subset_names,]) # get data subset
  
  
  ########################################
  # Take average, of replicates if day 392
  ########################################
  if (DAY == 7) {
    SUBSET_7 <- SUBSET
  }
  if (DAY == 392) {
    patient_protocol <- unique(paste(sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][1])}), sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][2])}), sep='-'))
    SUBSET_MEAN <- NULL
    new_rownames <- NULL
    for (pp in patient_protocol) {
      SUBSET_MEAN <- rbind(SUBSET_MEAN, t(colMeans(SUBSET[grepl(pp, rownames(SUBSET)),,drop=F])))
      new_rownames <- c(new_rownames, grep(pp, rownames(SUBSET), value = T)[1])
    }
    rownames(SUBSET_MEAN) <- new_rownames
    SUBSET <- data.frame(SUBSET_MEAN)
    SUBSET_392 <- SUBSET
  }
  
  
  #####################################
  # Identifiy Protocol, Patient and Day
  #####################################
  SUBSET$Patient <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][1])})
  SUBSET$Protocol <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][2])})
  SUBSET$Day <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][3])})
  SUBSET$Replicate <- sapply(rownames(SUBSET),function(x) {return(strsplit(x,'-')[[1]][4])})
  
  
  #############################
  # Potentially filter features
  #############################
  SUBSET <- SUBSET[,colSums(SUBSET > 0) > 0 | colnames(SUBSET)=='Replicate'] # Include only features with at least 1 non-zero value
  
  ######################
  # Create empty vectors
  ######################
  p_values <- vector(length = ncol(SUBSET)-4) # create empty vectors to store data
  q_ratios <- vector(length = ncol(SUBSET)-4) # create empty vectors to store data
  
  
  ####################
  # Loop over features
  ####################
  for (c in 1:(ncol(SUBSET)-4)) { # does not consider the last 4 columns that include the patient, protocol, day
    subset <- SUBSET[,c(c,(dim(SUBSET)[2]-3):dim(SUBSET)[2])] # world's ugliest way of extracting last three columns and the specific species we want to test
    colnames(subset)[1] <- 'Input'
    r <- friedman.test(formula( paste(colnames(SUBSET)[c], '~ Protocol | Patient' )), SUBSET) # run test
    p_values[c] <- r$p.value # save values
    q_ratios[c] <- r$statistic # save values
  } 
  
  
  #################
  # Adjust p values
  #################
  p_adjusted <- p.adjust(p_values, method = "BH") # adjust p-values
  
  
  ################################
  # Create data frame for plotting
  ################################
  df <- as.data.frame(p_adjusted) 
  df$q_ratios <- q_ratios
  FOR_PLOTTING[[DAY]] <- df
  
}

##################
# Merge dataframes
##################
df <- rbind(cbind(FOR_PLOTTING[[7]], 7))
colnames(df)[3] <- 392            
df <- rbind(df, cbind(FOR_PLOTTING[[392]], 392))
colnames(df)[3] <- 'day'


######
# plot
######

# This only works within MOCAT
pdf(file=paste(FIGURES.DIR, '/', 'fdr.ratio.pdf', sep=''), width=10, height=10)

# Set colors
day7_color <- 'red'
day392_color <- 'blue'
day_diff_color <- 'black'
colors <- c(day7_color, day392_color, day_diff_color)


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
plot_right <- ggplot() + geom_density(aes(df$q_ratios, fill=factor(df$day) ), alpha=0.5) + coord_flip() + theme(legend.position = "none") + xlim(min(df[,2]),max(df[,2])) +
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
  ) + scale_fill_manual(values = colors)

# Top density
plot_top <- ggplot() + geom_density(aes(df$p_adjusted, fill=factor(df$day) ), alpha=0.5) + theme(legend.position = "none") + xlim(0,1) +
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
  ) + scale_fill_manual(values = colors)

# Main
main <- ggplot(df, aes(p_adjusted,q_ratios),label=df$name,) + geom_point( aes(colour=factor(df$day)), size=3) +
  ylab("Ratio of variance explained: method effects over subject effects") +
  xlab("FDR adjusted p-values") +
  geom_hline(yintercept=1, colour='red', size = 1, linetype="dotted") +
  geom_vline(xintercept=0.05, colour='red', size = 1, linetype="dotted") +
  theme(plot.margin=unit(c(0,-0.5,0,0), "cm"),
        legend.position="none"
  ) + xlim(0,1) + ylim(min(df[,2]),max(df[,2])) + scale_fill_manual(values = colors) + scale_colour_manual(values = colors)
  

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
COMBINED <- rbind(SUBSET_7, SUBSET_392)
patients <- unique(sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][1])}))
protocols <- sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][2])})
days <- sapply(rownames(COMBINED),function(x) {return(strsplit(x,'-')[[1]][3])})
both <- paste(protocols, days, sep='-')
both_u <- unique(both)
both_u7 <-grep('-7', both_u, value=T)
both_u392 <-grep('-392', both_u, value=T)

# Loop over protocols
COR <- NULL
for (d in c(7, 392)){
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
          p1 <- gsub('11', 'Frozen', p1)
          p1 <- gsub('12', 'RNALater (RT)', p1)
          p1 <- gsub('13', 'RNALater (4 C)', p1)
          p2 <- gsub('11', 'Frozen', p2)
          p2 <- gsub('12', 'RNALater (RT)', p2)
          p2 <- gsub('13', 'RNALater (4 C)', p2)
          
          COR <- rbind(COR, c(paste('',p1,'vs\n', p2), d, patient, cor(as.numeric(x),as.numeric(y),method = "spearman")))
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
      COR <- rbind(COR, c(paste('DAY 7 vs 392'), '5', patient, cor(as.numeric(x),as.numeric(y),method = "spearman")))
    }
  }
}
COR <- data.frame(COR)
colnames(COR) <- c('protocol', 'day', 'patient', 'value')
COR$value <- as.numeric(as.vector(COR$value))
COR$f <- factor(rep(type,dim(COR),dim(COR)[1]))

# Create funciton for median line
f <- function(x, height) {
  ans <- median(x)
  data.frame(ymin = ans-height/2, ymax = ans+height/2, y = ans)
}

# Plot
plot2 <- ggplot(COR, aes(factor(protocol), as.numeric(value))) + geom_jitter(size=8, aes(color=factor(COR$day))) + 
  stat_summary(fun.data = f, geom = "crossbar", height = 0.008,
               colour = NA, fill = "black", width = 0.8, alpha = 1)+facet_grid(.~f) + theme_bw() +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)
        ) +
  theme(strip.text.x = element_text(size = 20, face="bold", angle = 0)) +
  xlab('') + ylab('Rank correlation within subjects') +
  scale_fill_manual(values = colors) + scale_colour_manual(values = colors)




grid.arrange(gA, gC,gB, gD, ncol=2, nrow=3, widths=c(4, 1), heights=c(1, 4))

par(mfrow=c(2,1))


grid.arrange(plot2,
             arrangeGrob(gA, gC,gB, gD, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)),
             nrow=2, heights=c(1,1)
)

grid.draw(a)


