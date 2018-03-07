# The idea is that using the previous scripts in the
# wilcox heatmap package we have identified a number of
# e.g. modules that we'd like to look closer at, and that's
# what this script does.Currently it's semi focused on module selection.

# specify the file to load, it should be a file
# that is the result of step two of this package.
FILE <- '/g/bork1/kultima/projects/m3/testset/maxmin10/TMP/out/2015Jun16_195233/data/session.RData'
load(FILE)

# THESE NEED TO BE ADDED TO SETTINGS
ko2pathway <- read.table('/g/bork/kultima/db/kegg/ko2pathway.txt')
ko2module <- read.table('/g/bork/kultima/db/kegg/ko2module.txt')
functional_map <- '/g/bork1/kultima/MOCAT/data/779.CC.RefGeneCatalog.padded.1-2.functional.map'
load_path <- '/g/bork1/kultima/projects/m3/testset/maxmin10/TMP/out/2015Jun16_195233/data/'
# this need to be run /g/bork8/kultima/GIT/MOCAT/master/public/src/MOCATFraction.pl -in gene.gene -out gene.gene.fraction


# Select a module
module <- 'M00080'

# needed function
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
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

require(reshape)
require(ggplot2)

# because this is a module, we can check the associated pathways
# we could either go via the genes that are significant, or do it via the KEGG map
# I prefer the first
file <- SETTINGS['significant_gene', 'setting']
genes <- read.table(file, header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char="", fill = T)
selected_genes <- genes[grepl(module, genes$module),] # get genes
pathways <- sort(table(selected_genes$pathway)) # get possible
biggest_pathway <- names(pathways[length(pathways)]) # get bigegst

# this here goes via the significant genes, but let's use the KEGG maps instead
#selected_pathways <- genes[grepl(biggest_pathway, genes$pathway),] # after finding biggest pathway, let's get the other modules of that pathway
all_kos_in_biggest_oathway <- unique(as.vector(ko2pathway[ko2pathway[,2] %in% biggest_pathway,1]))
# then let's get the modules assoc with those KOs
all_modules_in_biggest_oathway <- unique(as.vector(ko2module[ko2module[,1] %in% all_kos_in_biggest_oathway,2]))
# let's get the list of significant kos and modules
file <- SETTINGS['significant_ko','setting']
significant_ko <- read.table(file, header=T, quote="", row.names=1, sep=",", check.names=F, na.strings="NA", comment.char="")
file <- SETTINGS['significant_module','setting']
significant_module <- read.table(file, header=T, quote="", row.names=1, sep=",", check.names=F, na.strings="NA", comment.char="")
file <- SETTINGS['significant_pathway','setting']
significant_pathway <- read.table(file, header=T, quote="", row.names=1, sep=",", check.names=F, na.strings="NA", comment.char="")
signif_ko <- rownames(significant_ko[significant_ko$p.adj.fdr <= 0.05,])
here_ko <- as.vector(ko2pathway[ko2pathway[,2] %in% biggest_pathway,1])
signif_ko <- intersect(signif_ko, here_ko)

signif_module <- rownames(significant_module[significant_module$p.adj.fdr <= 0.05,])
here_module <- unique(as.vector(ko2module[ko2module[,1] %in% signif_ko,2]))
signif_module <- intersect(signif_module, here_module)

signif_pathway <- rownames(significant_pathway[significant_pathway$p.adj.fdr <= 0.05,])
# let's get the insignificant ones
insignif_ko <- all_kos_in_biggest_oathway[!all_kos_in_biggest_oathway %in% signif_ko]
insignif_module <- all_modules_in_biggest_oathway[!all_modules_in_biggest_oathway %in% signif_module]

# here we have to extract all the genes for the specific pathway and also the the map
system( paste("grep -v '^#' ",  load_path, "gene.gene.fraction | head -1 > ", SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.h', sep=''  ))
system( paste( "grep -w ", biggest_pathway, " ", functional_map, " | cut -f 1 | sed 's/$/\t/' | fgrep -f - ", load_path, '/gene.gene.fraction > ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.1', sep='' ) )
system( paste( "cat ", SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.h ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.1 > ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, sep='') )
system( paste( "grep -w ", biggest_pathway, " ", functional_map, " > ", SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.map', sep='') )
system( paste( "head -1 ", functional_map, " > ", SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.map.header', sep='') )
system( paste( "grep -w ", biggest_pathway, " ", functional_map, " | cut -f 1 | sed 's/$/\t/' | fgrep -f - ", load_path, '/gene.gene.horizontal > ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.2', sep='' ) )
system( paste( "cat ", SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.h ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.2 > ', SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.hor', sep='') )

gene_map <- paste(SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.map', sep='')
gene_map_header <- paste(SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.map.header', sep='')
all_genes <- paste(SETTINGS['data_dir', 'setting'], "/", biggest_pathway, sep='')
all_genes_horizontal <- paste(SETTINGS['data_dir', 'setting'], "/", biggest_pathway, '.hor', sep='')
# load the gene map
gene_map <- read.table(gene_map, header=F, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA")
gene_map_header <- colnames(read.table(gene_map_header, header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char=""))
colnames(gene_map) <- gene_map_header
# load the genes
all_genes <- read.table(all_genes, header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA")
all_genes_horizontal <- read.table(all_genes_horizontal, header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA")


# this could probablyy be skipped -->
# tmp <- as.vector(unique(selected_pathways$module))
# tmp <- tmp[tmp != '']
# all_modules_in_top_pathway <- genes[genes$module %in% tmp,]
# tmp <- unique(all_modules_in_top_pathway$module)
# unique_modules <- as.vector(tmp[tmp != ''])
# unique_modules <- unique(unlist(strsplit(unique_modules,',')))
# tmp <- unique(all_modules_in_top_pathway$ko)
# unique_kos <- as.vector(tmp[tmp != ''])
# unique_kos <- unique(unlist(strsplit(unique_kos,',')))
# <-- this could probablyy be skipped
# ok, so now we have:
# 1. a list of unique modules in the biggest pathway: unique_modules
# 2. the biggest pathway: biggest_pathway
# 3. The complete list of genes that all this boils down to: all_modules_in_top_pathway
#    and I guess these together take up 1000s of genes
# 4. list of unique KOs: unique_kos
# Let's try and plot this in some kind of nice way...

# Load some data
ko <- t(DATAlist[['ko']])
module <- t(DATAlist[['module']])
pathway <- t(DATAlist[['pathway']])
gene <- t(DATAlist[['gene']])
Hgene <- t(DATAlist[['gene.horizontal']])
Hko <- t(DATAlist[['ko.horizontal']])
Hmodule <- t(DATAlist[['module.horizontal']])
Hpathway <- t(DATAlist[['pathway.horizontal']])

# It's going to take longer to do it sample, by sample, but I think it'll be easier
# it's currently too complicatedto get it down to gene level, let's do it for pathway, module and ko

# these are not part of this pathway
exclude_ko <- setdiff(unique(unlist(strsplit(as.vector(gene_map$ko),","))), rownames(ko))
exclude_module <- setdiff(unique(unlist(strsplit(as.vector(gene_map$module),","))), rownames(module))

PLUS <- 0
INDEX <- 0
MAX <- 0
SAMPLES <- 0
CENTER <- NULL
NAMES <- NULL
A4 <- list()
COLOR <- list()
AiA2 <- list()
PLOTS <- list()
cat(paste('Processing (this wil take a few minutes)', '\n'))
for (sample in c(g1[1:2], g2[1:2])){
  cat (paste('Processing', sample), '\n')
  ALLINALL <- NULL
  ALLINALL2 <- NULL
  a2 <- NULL
  a2b <- NULL
  a3 <- NULL
  a4 <- NULL
  INDEX <- INDEX + 1
  #sample <- 'CCIS31434951ST-20-0'
  #gene_abundance <- gene[rownames(all_modules_in_top_pathway), sample]
  #set <- cbind(all_modules_in_top_pathway, gene_abundance)
  
  # we need to process all KOs and modules, and keep track of which ones are significant
  
  # all this is probably not used...
  #     s_ko <- ko[rownames(ko) %in% signif_ko,sample,drop=F] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  #     ns_ko <- ko[rownames(ko) %in% insignif_ko,sample,drop=F]
  #     Hs_ko <- Hko[rownames(Hko) %in% rownames(s_ko),sample,drop=F] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  #     Hns_ko <- Hko[rownames(Hko) %in% rownames(ns_ko),sample,drop=F]
  #     a_ko <- rbind(s_ko, ns_ko)
  #     Ha_ko <- rbind(Hs_ko, Hns_ko)
  #     
  #     s_mod <- module[rownames(module) %in% signif_module,sample,drop=F] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  #     ns_mod <- module[rownames(module) %in% insignif_module,sample,drop=F]
  #     Hs_mod <- Hmodule[rownames(Hmodule) %in% rownames(s_mod),sample,drop=F] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  #     Hns_mod <- Hmodule[rownames(Hmodule) %in% rownames(ns_mod),sample,drop=F]
  #     a_mod <- rbind(s_mod, ns_mod)
  #     Ha_mod <- rbind(Hs_mod, Hns_mod)
  #     
  g <- rownames(gene_map)[grepl( paste(signif_ko, collapse = "|"), gene_map$ko)]
  s_gene <- all_genes[g,sample,drop=F]
  g <- rownames(gene_map)[grepl( paste(insignif_ko, collapse = "|"), gene_map$ko)]
  ns_gene <- all_genes[g,sample,drop=F]
  Hs_gene <- all_genes_horizontal[rownames(s_gene),sample,drop=F]
  Hns_gene <- all_genes_horizontal[rownames(ns_gene),sample,drop=F]
  
  
  if(biggest_pathway %in% signif_pathway) {
    pathway_s <- TRUE
    whether_significant = 'significant'
  } else {
    pathway_s <- FALSE
    whether_significant = 'not significant'
  }
  
  p1 <- cbind(s_gene, pathway[biggest_pathway,sample], Hs_gene, Hpathway[biggest_pathway,sample], 'A')
  rownames(p1) <- paste('SIGNIF', sample, biggest_pathway, rownames(p1))
  colnames(p1) <- c('GENE-C', 'PATHWAY-C', 'GENE-H', 'PATHWAY-H', 'GENE-SIGNIF')
  
  p2 <- cbind(ns_gene, pathway[biggest_pathway,sample], Hns_gene, Hpathway[biggest_pathway,sample], 'B')
  rownames(p2) <- paste('NON-SIGNIF', sample, biggest_pathway, rownames(p2))
  colnames(p2) <- c('GENE-C', 'PATHWAY-C', 'GENE-H', 'PATHWAY-H', 'GENE-SIGNIF')
  
  
  all <- rbind(p1, p2)
  map <- gene_map[c(rownames(s_gene), rownames(ns_gene)),]
  # let's loop over the genes, this may be slow, but I know what's going on
  for (cg in c(1:dim(map)[1])) { # loop first over significant and then non significant
    l <- gene_map[cg,]
    k <- unlist(strsplit(as.vector(l$module), ",")) # loop over each module
    kk <- unlist(strsplit(as.vector(l$ko), ","))
    for (k1 in k) {
      if(!k1 %in% exclude_module){ # only proceed if this module is in this pathway
        v <- module[k1,sample]
        Hv <- Hmodule[k1,sample]
        for (k2 in kk) {
          if(!k2 %in% exclude_ko){ # only proceed if this ko is in this pathway
            vv <- ko[k2,sample]
            Hvv <- Hko[k2,sample]
            if(k1 %in% signif_module){
              MS <- 'A'
            } else {
              MS <-'B'
            }
            if(k2 %in% signif_ko){
              KS <- 'A'
            } else {
              KS <- 'B'
            }
            ALLINALL <- rbind(ALLINALL,
                              cbind(all[cg,], k1, k2, v, Hv, vv, Hvv, MS, KS))
          }
        }
      }
    }
  } # end loop over each gene line
  
  
  
  
  
  #ALLINALL$index <- c(1:dim(ALLINALL)[1])
  colnames(ALLINALL) <- c('GENE-C', 'PATHWAY-C', 'GENE-H', 'PATHWAY-H', 'GENE-SIGNIF', 'MODULE', 'KO', 'MODULE-C', 'MODULE-H', 'KO-C', 'KO-H', 'MODULE-SIGNIF', 'KO-SIGNIF')
  ALLINALL2 <- ALLINALL[, c('GENE-C', 'PATHWAY-C', 'GENE-H', 'PATHWAY-H', 'GENE-SIGNIF', 'MODULE-C', 'MODULE-H', 'KO-C', 'KO-H', 'MODULE-SIGNIF', 'KO-SIGNIF')]
  rownames(ALLINALL2) <- paste(rownames(ALLINALL2), ALLINALL$KO, ALLINALL$MODULE)
  #O2 <- paste(ORDER, as.vector(ALLINALL$KO), as.vector(ALLINALL$MODULE))
  
  
  
  ORDER <- rownames(ALLINALL2)[sort(paste(ALLINALL$`MODULE-SIGNIF`, ALLINALL$MODULE, ALLINALL$`KO-SIGNIF`, ALLINALL$KO), index.return=T)$ix]
  ALLINALL2 <- ALLINALL2[ORDER,]
  ALLINALL2$`GENE-SIGNIF` <- as.vector(ALLINALL2$`GENE-SIGNIF`)
  ALLINALL2$`MODULE-SIGNIF` <- as.vector(ALLINALL2$`MODULE-SIGNIF`)
  ALLINALL2$`KO-SIGNIF` <- as.vector(ALLINALL2$`KO-SIGNIF`)
  ALLINALL2[ALLINALL2 == 'B'] <- 0 #before if we use the two color scale, this would be 0 and below be NA
  ALLINALL2[ALLINALL2 == 'A'] <- 100
  ALLINALL2b <- ALLINALL2
  ALLINALL2$index <- c((1):(dim(ALLINALL2)[1]))
  ALLINALL2b$N <- rownames(ALLINALL2) # because we use index, I don't use N anymore
  
  # rescale
  # this has to be done with the regards to the max in this sample for each category
  # they should all just be multiplied by 100 because we go from fractions to %
  # NO, we rescale so that the highest coverage value would be similar to the highest horizontal value
  # NO, again, we have to do it on per individual category
  #max <- max(c(ALLINALL2$`GENE-C`, ALLINALL2$`KO-C`, ALLINALL2$`MODULE-C`, ALLINALL2$`PATHWAY-C`),na.rm = T)
  #rescale <- 100/max
  #  ALLINALL2$`GENE-C` <- -ALLINALL2$`GENE-C` * max(ALLINALL2$`GENE-H`,na.rm = T) / max(ALLINALL2$`GENE-C`, na.rm=T)
  #  ALLINALL2$`MODULE-C` <- -ALLINALL2$`MODULE-C` * max(ALLINALL2$`MODULE-H`, na.rm=T) / max(ALLINALL2$`MODULE-C`, na.rm=T)
  #  ALLINALL2$`KO-C` <- -ALLINALL2$`KO-C` * max(ALLINALL2$`KO-H`, na.rm=T) / max(ALLINALL2$`KO-C`, na.rm=T)
  #  ALLINALL2$`PATHWAY-C` <- -ALLINALL2$`PATHWAY-C` * max(ALLINALL2$`PATHWAY-H`, na.rm=T) / max(ALLINALL2$`PATHWAY-C`, na.rm=T)
  
  ALLINALL2$`GENE-C` <- ALLINALL2$`GENE-C` * 100 / max(ALLINALL2$`GENE-C`, na.rm=T)
  ALLINALL2$`MODULE-C` <- ALLINALL2$`MODULE-C` * 100 / max(ALLINALL2$`MODULE-C`, na.rm=T)
  ALLINALL2$`KO-C` <- ALLINALL2$`KO-C` * 100 / max(ALLINALL2$`KO-C`, na.rm=T)
  ALLINALL2$`PATHWAY-C` <- ALLINALL2$`PATHWAY-C` * 100 / max(ALLINALL2$`PATHWAY-C`, na.rm=T)
  
  # ALLINALL2$`GENE-C` <- -ALLINALL2$`GENE-C` * rescale
  # ALLINALL2$`MODULE-C` <- -ALLINALL2$`MODULE-C` * rescale
  # ALLINALL2$`KO-C` <- -ALLINALL2$`KO-C` * rescale
  # ALLINALL2$`PATHWAY-C` <- -ALLINALL2$`PATHWAY-C` * rescale
  
  
  # melt
  #a2 <- melt(ALLINALL2, id='N')
  a2 <- melt(ALLINALL2, id='index')
  a2b <- melt(ALLINALL2b, id='N')
  a2$value <- as.numeric(a2$value)
  
  sp <- t(matrix( unlist(strsplit(a2b$N, " ")) , nrow = 6))
  a2$mod <- sp[,6]
  a2$ko <- sp[,5]
  
  # we need to make the colors distinguished between different levels
  a2$C <- NA
  a2[a2$variable == 'PATHWAY-C','C'] <- biggest_pathway
  a2[a2$variable == 'PATHWAY-H','C'] <- biggest_pathway
  a2[a2$variable == 'MODULE-C','C'] <- a2[a2$variable == 'MODULE-C','mod']
  a2[a2$variable == 'MODULE-H','C'] <- a2[a2$variable == 'MODULE-H','mod']
  a2[a2$variable == 'MODULE-SIGNIF','C'] <- a2[a2$variable == 'MODULE-SIGNIF','mod']
  a2[a2$variable == 'KO-C','C'] <- a2[a2$variable == 'KO-C','ko']
  a2[a2$variable == 'KO-H','C'] <- a2[a2$variable == 'KO-H','ko']
  a2[a2$variable == 'KO-SIGNIF','C'] <- a2[a2$variable == 'KO-SIGNIF','ko']
  a2[a2$variable == 'GENE-C','C'] <- a2[a2$variable == 'GENE-C','ko']
  a2[a2$variable == 'GENE-H','C'] <- a2[a2$variable == 'GENE-H','ko']
  a2[a2$variable == 'GENE-SIGNIF','C'] <- a2[a2$variable == 'GENE-SIGNIF','ko']
  
  # this is an attempt to change heights...
  
  a2$value <- a2$value/100
  
  a3 <- a2
  a3b <- a2
  a3$variable <- gsub('-C', '-C2', a3$variable)
  a3$variable <- gsub('-H', '-H2', a3$variable)
  a3$variable <- gsub('-SIGNIF', '-SIGNIF2', a3$variable)
  a3$C <- a3$variable
  a3$C <- gsub('.*-', '', a3$C, perl=T)
  a3b$value <- 1
  a4 <- rbind(a3b,a3)
  
  a4$H <- 1
  
  #a4[a4$variable == 'GENE-C', 'H'] <- 0.5
  
  a4$Y <- NA
  
  a4[a4$variable == 'PATHWAY-C','Y'] <- 11
  a4[a4$variable == 'PATHWAY-C2','Y'] <- 10.4
  a4[a4$variable == 'PATHWAY-H2','Y'] <- 9.4
  
  a4[a4$variable == 'MODULE-C','Y'] <- 8.2
  a4[a4$variable == 'MODULE-C2','Y'] <- 7.6
  a4[a4$variable == 'MODULE-H2','Y'] <- 6.6
  a4[a4$variable == 'MODULE-SIGNIF2','Y'] <- 6
  
  a4[a4$variable == 'KO-C','Y'] <- 5.2
  a4[a4$variable == 'KO-C2','Y'] <- 4.6
  a4[a4$variable == 'KO-H2','Y'] <- 3.6
  a4[a4$variable == 'KO-SIGNIF2','Y'] <- 3
  
  a4[a4$variable == 'GENE-C','Y'] <- 2.2
  a4[a4$variable == 'GENE-C2','Y'] <- 1.6
  a4[a4$variable == 'GENE-H2','Y'] <- 0.6
  a4[a4$variable == 'GENE-SIGNIF2','Y'] <- 0
  
  
  a4[a4$variable == 'PATHWAY-C','H'] <- 0.2
  a4[a4$variable == 'PATHWAY-C2','H'] <- 1
  a4[a4$variable == 'PATHWAY-H2','H'] <- 1
  
  a4[a4$variable == 'MODULE-C','H'] <- 0.2
  a4[a4$variable == 'MODULE-C2','H'] <- 1
  a4[a4$variable == 'MODULE-H2','H'] <- 1
  a4[a4$variable == 'MODULE-SIGNIF2','H'] <- 0.2
  
  a4[a4$variable == 'KO-C','H'] <- 0.2
  a4[a4$variable == 'KO-C2','H'] <- 1
  a4[a4$variable == 'KO-H2','H'] <- 1
  a4[a4$variable == 'KO-SIGNIF2','H'] <- 0.2
  
  a4[a4$variable == 'GENE-C','H'] <- 0.2
  a4[a4$variable == 'GENE-C2','H'] <- 1
  a4[a4$variable == 'GENE-H2','H'] <- 1
  a4[a4$variable == 'GENE-SIGNIF2','H'] <- 0.2
  
  
  # to get the colors in right order
  a4$C <- gsub("C2", ' Coverage', a4$C)
  a4$C <- gsub("H2", ' Horiz. Cov.', a4$C)
  a4$C <- gsub("SIGNIF2", ' Significance', a4$C)
  a4$C <- factor(a4$C)
  COLORS <- c('black', 'black', 'black', (ggplotColours(length(levels(a4$C))-3)))
  
  # we need to get the order of the colors
  # get pathway
  col_order1 <- as.vector((a4[a4$variable == 'PATHWAY-C','C'])[1])
  col_order2 <- (as.vector((a4[a4$variable == 'MODULE-C','C'])))
  col_order2 <- col_order2[!duplicated(col_order2)]
  col_order3 <- (as.vector((a4[a4$variable == 'KO-C','C'])))
  col_order3 <- col_order3[!duplicated(col_order3)]
  col_order <- c(' Coverage', ' Horiz. Cov.', ' Significance', col_order1, col_order2, col_order3)
  #current_order <- levels(a4$C)
  names(COLORS) <- col_order
  a4$C <- factor(a4$C, levels=col_order)
  #COLORS <- COLORS [ current_order ]
  
  # order4
  #ORDER2 <- rev(c('PATHWAY-C', 'PATHWAY-H', 'MODULE-C', 'MODULE-H', 'MODULE-SIGNIF', 'KO-C', 'KO-H', 'KO-SIGNIF', 'GENE-C', 'GENE-H', 'GENE-SIGNIF'))
  ORDER4 <- rev(c('PATHWAY-C', 'PATHWAY-C2', 'PATHWAY-H2', 'MODULE-C','MODULE-C2', 'MODULE-H2', 'MODULE-SIGNIF2', 'KO-C','KO-C2', 'KO-H2', 'KO-SIGNIF2', 'GENE-C','GENE-C2', 'GENE-H2', 'GENE-SIGNIF2'))
  
  #ORDER2 <- rev(c('PATHWAY-C', 'PATHWAY-H', 'MODULE-C', 'MODULE-H', 'KO-C', 'KO-H', 'GENE-C', 'GENE-H'))
  # plot
  
  # save data to keep order from first sample
  if(INDEX == 1) {
    COLORS_MAIN <- COLORS
    MIN <- dim(ALLINALL2)[1]*0.1
  }
  PLUS <- PLUS + dim(ALLINALL2)[1] + 200
  MAX <- dim(ALLINALL2)[1] + MAX
  SAMPLES <- SAMPLES + 1
  NAMES <- c(NAMES, sample)
  
  A4[[sample]] <- a4
  COLOR[[sample]] <- COLORS
  AiA2[[sample]] <- ALLINALL2
  
  size <- 7
  
  PLOTS[[sample]] <- ggplot() +
    geom_tile(data=a4, aes(index, Y, fill=C, alpha=value, height=H), width=1, colour=NA) +
    scale_fill_manual(values = COLORS, name="") +
    xlim(-dim(ALLINALL2)[1]*0.1, dim(ALLINALL2)[1]) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(), 
          axis.title=element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    scale_alpha(guide = 'none') +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 0, label = "Significance", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 3, label = "Significance", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 6, label = "Significance", family='Courier', fontface='bold', size=size) +
    
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 0.6, label = "Horiz. Cov.", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 3.6, label = "Horiz. Cov.", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 6.6, label = "Horiz. Cov.", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 9.4, label = "Horiz. Cov.", family='Courier', fontface='bold', size=size) +
    
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 1.6, label = "Coverage", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 4.6, label = "Coverage", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 7.6, label = "Coverage", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 10.4, label = "Coverage", family='Courier', fontface='bold', size=size) +
    
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 11, label = paste(biggest_pathway, sep=''), family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 8.2, label = "Modules", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 5.2, label = "KOs", family='Courier', fontface='bold', size=size) +
    annotate("text", x = -dim(ALLINALL2)[1]*0.08, y = 2.2, label = "KOs", family='Courier', fontface='bold', size=size) +
    
    annotate("text", x = dim(ALLINALL2)[1]/2, y = -0.2, label = "< Genes >", family='Courier', fontface='bold', size=size) +
    annotate("text", x = dim(ALLINALL2)[1]/2, y = 2.8, label = "< KOs >", family='Courier', fontface='bold', size=size) +
    annotate("text", x = dim(ALLINALL2)[1]/2, y = 5.8, label = "< Modules >", family='Courier', fontface='bold', size=size) +
    annotate("text", x = dim(ALLINALL2)[1]/2, y = 8.8, label = paste("<", biggest_pathway,">"), family='Courier', fontface='bold', size=size) +
    
    annotate("text", x = dim(ALLINALL2)[1]/2, y = 11.4,
             label = paste('Sample ', sample,'\nOverview of Pathway ',
                           biggest_pathway, '\n(', whether_significant ,'ly differentially abundant)', sep=''),
             family='Courier', fontface='bold', size=size) +
    
    theme(
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.background=element_rect(fill="white", colour=NA),
      legend.key=element_rect(colour="white"),
      legend.title=element_text(size=rel(1), face="bold", family='Courier', hjust=0),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
      strip.background=element_blank()
    ) +  
    theme(legend.text = element_text(size = rel(1.5), colour = "black"))
  
  
} # end loop over samples



jpeg(filename = 'final.jpeg', width = 1000*1.5 * SAMPLES, height = 3000*1.5)
multiplot(PLOTS[[1]],PLOTS[[2]],PLOTS[[3]],PLOTS[[4]], cols = length(PLOTS) )
dev.off()

# OLD ONES
# png(filename = 'b.png', width = 1000, height = 3000)
# ggplot() +
#   geom_tile(data=a2, aes(N, variable, alpha = factor(value), fill=factor(paste(C))), colour = NA) +
#   #scale_fill_gradient2(low ="darkorange4", mid = "white", high = "steelblue4", na='red') +
#   scale_x_discrete(limits=(rownames(ALLINALL2))) +
#   scale_y_discrete(limits=unique(ORDER2))
# #scale_alpha_manual(values = sort(unique(a2$value)) )
# dev.off()
# 
# # this plots in the two different colors and uses negative value above for the horiztonal
# png(filename = 'e.png', width = 1000, height = 3000)
# ggplot() +
#   geom_tile(data=a4, aes(N, variable, fill=factor(C), alpha=value, height=H), colour=NA) +
#   #scale_fill_gradient2(low ="darkorange4", mid = "white", high = "steelblue4", na='red') +
#   scale_x_discrete(limits=(rownames(ALLINALL2))) +
#   scale_y_discrete(limits=unique(ORDER4)) 
# dev.off()
#  geom_tile(data=bmh, aes(X1, X2, fill = value*2), colour = 'red', na.value='black')

# theme(axis.text.x = element_text(angle = 25, hjust = 1, color="black")) +
#   #theme(legend.justification=c(0,1), legend.position=c(0,0.90)) +
#   theme(legend.justification=c(0,1), legend.position="none") +
#   xlab(paste(g1n, ' - ', g2n, ' | orange=coverage', sep='')) +
#   theme(plot.title = element_text(size = rel(3))) +
#   theme(axis.text.x = element_text(size = rel(3))) + 
#   theme(axis.title.x = element_text(size = rel(3))) +
#   theme(axis.text.y = element_text(size = rel(3))) + 
#   theme(axis.title.y = element_text(size = rel(3))) +
#   theme(legend.text = element_text(size = rel(0))) +
#   theme(strip.text = element_text(size = rel(0))) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.title.y=element_blank()) +
#   theme(axis.ticks = element_blank(), axis.text.x = element_blank())











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
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
