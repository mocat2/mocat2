
# THIS DATA COMES FROM THE server.R
selected_types <- as.vector(input$type)
print_species <- input$printSpecies
exclude_phylum <- input$excludePhylum
exclude_genus <- input$excludeGenus
exclude_species <- input$excludeSpecies
X <- strsplit(as.vector(input$module), " \\(")
selected_modules <- sapply(  X  , "[[", 1)
X <- strsplit(as.vector(input$KO), " \\(")
selected_kos <- sapply(  X   , "[[", 1)
X <- strsplit(as.vector(input$pathway), " \\(")
selected_pathways <- sapply(  X  , "[[", 1)
#selected_pathways <- input$pathway
#selected_modules <- as.vector(input$module)
#selected_kos <- as.vector(input$KO)
# THIS DATA COMES FROM THE server.R

# input <- list()
# input$module <- 'NO MODULE'
# input$KO <- 'KO00677'
# selected_pathways <- 'ko00540'

# GENE_WIDTH <- 20
# GENE_HEIGHT <- 10
# HEIGHT_KO_TAXA <- 3
# KO_HEIGHT <- 50
# MOD_HEIGHT <- 30
# PATHWAY_HEIGHT <- 30
# Y_START_GENES <- 0
# SPACE_BW_GENE_ROWS <- 40
# SPACE_FOR_EACH_TAXA <- 8
# SPACE_BW_KOS <- 10
# SPACE_BW_PATHWAYS <- 3000
# SPACE_BW_MODULES <- 1000
# SPACE_BW_GROUPS <- 3000
# MAX_GENES <- 50
# KO_ABOVE_GENE <- 40
# MOD_ABOVE_GENE <- 80
# PATHWAY_ABOVE_GENE <- 120
# SPACE_BW_KO_COV_HOR <- 0.2
# selected_types <- c('Gene coverages',
#                     'Horizontal gene coverages',
#                     'KO coverages',
#                     'Horizontal KO coverages',
#                     'Module coverages',
#                     'Horizontal module coverages',
#                     'Gene species abundances',
#                     'Gene genus abundances',
#                     'Gene phylum abundances',
#                     'Pathway species abundances',
#                     'Pathway genus abundances',
#                     'Pathway phylum abundances',
#                     'KO species abundances',
#                     'KO genus abundances',
#                     'KO phylum abundances',
#                     'Module species abundances',
#                     'Module genus abundances',
#                     'Module phylum abundances')
# print_species <- FALSE

GENE_WIDTH <- as.numeric(as.vector(input$GENE_WIDTH))
GENE_HEIGHT <- as.numeric(as.vector(input$GENE_HEIGHT))
HEIGHT_KO_TAXA <- 3
KO_HEIGHT <- as.numeric(as.vector(input$KO_HEIGHT))
MOD_HEIGHT <- as.numeric(as.vector(input$MOD_HEIGHT))
Y_START_GENES <- 0
SPACE_BW_GENE_ROWS <- as.numeric(as.vector(input$SPACE_BW_GENE_ROWS))
SPACE_FOR_EACH_TAXA <- as.numeric(as.vector(input$SPACE_FOR_EACH_TAXA))
SPACE_BW_KOS <- as.numeric(as.vector(input$SPACE_BW_KOS))
SPACE_BW_GROUPS <- as.numeric(as.vector(input$SPACE_BW_GROUPS))
SPACE_BW_MODULES <- as.numeric(as.vector(input$SPACE_BW_MODULES))
MAX_GENES <- as.numeric(as.vector(input$MAX_GENES))
KO_ABOVE_GENE <- as.numeric(as.vector(input$KO_ABOVE_GENE))
MOD_ABOVE_GENE <- as.numeric(as.vector(input$MOD_ABOVE_GENE))
PATHWAY_ABOVE_GENE <- as.numeric(as.vector(input$PATHWAY_ABOVE_GENE))
SPACE_BW_KO_COV_HOR <- as.numeric(as.vector(input$SPACE_BW_KO_COV_HOR))

TO_PLOT <- NULL
ko_original <- ko
Hko_original <- Hko
module_original <- module
Hmodule_original <- Hmodule
pathway_original <- pathway
Hpathway_original <- pathway
all_genes_original <- all_genes
all_genes_horizontal_original <- all_genes_horizontal

X_START <- 0
X_STOP <- 0

for (J in c(1:length(names(SAMPLE.GROUPS)))){
  
  to_plot_hor_coverage_gene_ALL <- NULL
  to_plot_coverage_gene_ALL <- NULL
  to_plot_coverage_ko_ALL <- NULL
  to_plot_hor_coverage_ko_ALL <- NULL
  to_plot_coverage_mod_ALL <- NULL
  to_plot_hor_coverage_mod_ALL <- NULL
  to_plot_coverage_pathway_ALL <- NULL
  to_plot_hor_coverage_pathway_ALL <- NULL
  to_plot_ko_s_ALL <- NULL
  to_plot_ko_g_ALL <- NULL
  to_plot_ko_p_ALL <- NULL
  to_plot_mod_s_ALL <- NULL
  to_plot_mod_g_ALL <- NULL
  to_plot_mod_p_ALL <- NULL
  to_plot_pathway_s_ALL <- NULL
  to_plot_pathway_g_ALL <- NULL
  to_plot_pathway_p_ALL <- NULL
  
  ko <- ko_original
  Hko <- Hko_original
  module <- module_original
  Hmodule <- Hmodule_original
  pathway <- pathway_original
  Hpathway <- pathway_original
  all_genes <- all_genes_original
  all_genes_horizontal <- all_genes_horizontal_original
  
  
  if( length(as.vector(input[[paste0("gr", J)]])  )) {
    selected_samples_here <- 'CCIS33816588ST-4-0'
    selected_samples_here <- as.vector(input[[paste0("gr", J)]]) 
    print(selected_samples_here)
    if(length(selected_samples_here) > 1) { # we take the average
      cat('Calculating average for samples.\n')
      sample <- paste('Average for', names(SAMPLE.GROUPS)[J])
      
      
      all_genes <- as.data.frame(t(t(rowSums(  all_genes[,selected_samples_here]  ))))
      all_genes_horizontal <- as.data.frame(t(t(rowSums(  all_genes_horizontal[,selected_samples_here]  ))))
      ko <- as.data.frame(t(t(rowSums(  ko[,selected_samples_here]  ))))
      Hko <- as.data.frame(t(t(rowSums(  Hko[,selected_samples_here]  ))))
      module <- as.data.frame(t(t(rowSums(  module[,selected_samples_here]  ))))
      Hmodule <- as.data.frame(t(t(rowSums(  Hmodule[,selected_samples_here]  ))))
      pathway <- as.data.frame(t(t(rowSums(  pathway[,selected_samples_here]  ))))
      Hpathway <- as.data.frame(t(t(rowSums(  Hpathway[,selected_samples_here]  ))))
      print(dim(ko))
      colnames(all_genes)[1] <- sample
      colnames(all_genes_horizontal)[1] <- sample
      colnames(ko)[1] <- sample
      colnames(Hko)[1] <- sample
      colnames(module)[1] <- sample
      colnames(Hmodule)[1] <- sample
      colnames(pathway)[1] <- sample
      colnames(Hpathway)[1] <- sample
    } else {
      sample <- selected_samples_here
    }
    
    # let's remvoe genes to excluded phylum
    TAXA <- as.matrix(TAXA)
    TAXA[is.na(TAXA)] <- 'Unknown'
    TAXA <- as.data.frame(TAXA)
    remove <- NULL
    
    cat('================\n')
    print(dim(all_genes))
    if(length(exclude_phylum)>0){
      cat('remove phyla\n')
      remove <- rownames(TAXA)[ TAXA[,3] %in% as.vector(exclude_phylum),drop=F ]
      if('Unknown' %in% as.vector(exclude_phylum)) {
        remove <- c(remove, rownames(all_genes)[!rownames(all_genes) %in% rownames(TAXA)])
      }
    }
    
    print(dim(all_genes))
    
    if(length(exclude_genus)>0){
      cat('remove genus\n')
      print(exclude_genus)
      print(dim(all_genes))
      remove <- c(remove, rownames(TAXA)[ TAXA[,7] %in% as.vector(exclude_genus),drop=F ])
      if('Unknown' %in% as.vector(exclude_genus)) {
        remove <- c(remove, rownames(all_genes)[!rownames(all_genes) %in% rownames(TAXA)])
      }
      print(dim(all_genes))
    }
    
    print(dim(all_genes))
    
    if(length(exclude_species)>0){
      cat('remove phyla\n')
      remove <- c(remove, rownames(TAXA)[ TAXA[,8] %in% as.vector(exclude_species),drop=F ])
      if('Unknown' %in% as.vector(exclude_species)) {
        remove <- c(remove, rownames(all_genes)[!rownames(all_genes) %in% rownames(TAXA)])
      }
    }
    
    
    if(length(remove)>0) {
      print(length(remove))
      
      all_genes <- as.data.frame(all_genes[ -which(rownames(all_genes) %in% unique(remove)), ,drop=F])
      all_genes_horizontal <- as.data.frame(all_genes_horizontal[ -which(rownames(all_genes) %in% unique(remove)), ,drop=F])
    }
    
    
    print(dim(all_genes))
    
    for (biggest_pathway in selected_pathways[1]) {
      
      current_pathway <- biggest_pathway
      print (current_pathway)
      to_plot_coverage_gene_current_pathway <- NULL
      PATHWAY_START <- X_START
      
      if (current_pathway == 'NO PATHWAY') {
        PATHWAY_SUM <- mean(ko[current_kos, sample]) * PATHWAY_HEIGHT/ max(pathway) # i guess technically this should include the NO MODULE too, but I'm cutting a coern but not including that in the max calcs
        PATHWAY_SUM_HOR <- mean(Hko[current_kos, sample]) * PATHWAY_HEIGHT/ max(Hpathway)
      } else {
        PATHWAY_SUM <- pathway[current_pathway, sample] * PATHWAY_HEIGHT/ max(pathway)
        PATHWAY_SUM_HOR <- Hmodule[current_pathway, sample] * PATHWAY_HEIGHT/ max(Hpathway)
      }
      
      
      # first we have to convert the horrible gene table
      # our pathway is: biggest_pathway
      # all the KOs: all_kos_in_biggest_pathway
      # all the modules: all_modules_in_biggest_pathway
      # all genes coverages: all_genes
      # all genes horizontal coverages: all_genes_horizontal
      
      # get nonzero genes
      all_genes_sample <- all_genes[all_genes[,sample,drop=F] > 0, sample, drop=F]
      all_genes_sample <- all_genes_sample * 100/max(all_genes) # need to rescale genes
      
      # first we order the modules with the most KOs
      all_kos_in_biggest_pathway <- unique(as.vector(ko2pathway[ko2pathway[,2] %in% as.vector(biggest_pathway),1]))
      all_modules_in_biggest_pathway <- unique(as.vector(ko2module[ko2module[,1] %in% all_kos_in_biggest_pathway,2]))
      order_module <- rev(sort(table(ko2module[ko2module[,2] %in% all_modules_in_biggest_pathway,2])))
      if ('NO MODULE' %in% names(order_module)) {
        order_module <- order_module[c(names(order_module[order_module > 0 & !names(order_module) == 'NO MODULE']), 'NO MODULE')]
      } else {
        order_module <- order_module[names(order_module[order_module > 0])]
      }
      
      
      # this is the step that before took a lot of time, and was annoying
      # we need to get a gene to ko map
      #       map <- cbind(rownames(gene_map), as.vector(gene_map$ko))
      #       map2a <- map[!grepl(',', map[,2]),]
      #       map2b <- map[grepl(',', map[,2]),]
      #       map2b[,1] <- paste(map2b[,1], map2b[,1], sep=',')
      #       map3a <- unlist(strsplit(map2a[,1], ","))  
      #       map3b <- unlist(strsplit(map2a[,2], ","))
      #       map4 <- cbind(map3a, map3b)
      #       gene2ko <- rbind(map2a, map4)
      #       gene2ko <- gene2ko[gene2ko[,1] %in% rownames(all_genes_sample),] # get only genes >0
      #       colnames(gene2ko) <- c('gene', 'ko')
      #       
      # then we need to get the order of the KOs
      
      
      order_ko <- rev(sort(table(gene2ko[gene2ko[,2] %in% unique(gene2ko[,2]),2]))) # here we pick only form KOs in this sample
      order_ko <- order_ko[order_ko > 0]
      order_ko <- order_ko[names(order_ko) %in% as.vector(ko2pathway[ko2pathway[,2] %in% biggest_pathway,1])]
      
      # combine the modules and kos to get a first order of them
      #order_module <- c(order_module, 0) # we have to add this for the KOs that don't have module
      #names(order_module) <- c(names(order_module)[1:length(order_module)-1], 'NO MODULE')
      #names(order_module) <- names(order_module)
      ones_we_have <- ko2module[ko2module[ ,1] %in% rownames(order_ko),]
      #ones_we_dont_have <- cbind(rownames(order_ko)[!rownames(order_ko) %in% ones_we_have[,1]], 'NO MODULE')
      #our_ko2module <- rbind(ones_we_have, ones_we_dont_have)
      our_ko2module <- ones_we_have
      
      # loop over modules to construct start and stops for genes, KOs and modules, note that we will exclude 0 abundant genes
      #for (current_module in names(order_module)[names(order_module) %in% selected_modules]) {
      #for (current_module in names(order_module[names(order_module) %in% selected_modules] )) {
      for (current_module in selected_modules[selected_modules %in% names(order_module)] ) {
        print (current_module)
        to_plot_coverage_gene_current_mod <- NULL
        current_kos <- as.vector(our_ko2module[our_ko2module[,2] == current_module,1]) # get the KOs
        current_kos <- rev(sort(order_ko[current_kos]))
        MOD_START <- X_START
        if (current_module == 'NO MODULE') {
          MODULE_SUM <- mean(ko[current_kos, sample]) * MOD_HEIGHT/ max(module) # i guess technically this should include the NO MODULE too, but I'm cutting a coern but not including that in the max calcs
          MODULE_SUM_HOR <- mean(Hko[current_kos, sample]) * MOD_HEIGHT/ max(Hmodule)
        } else {
          MODULE_SUM <- module[current_module, sample] * MOD_HEIGHT/ max(module)
          MODULE_SUM_HOR <- Hmodule[current_module, sample] * MOD_HEIGHT/ max(Hmodule)
        }
        
        
        
        # loop over the KOs and see how many genes are above zero and save these in a table
        #         for (current_ko in names(current_kos[names(current_kos) %in% selected_kos])) {
        for (current_ko in selected_kos[selected_kos %in% names(order_ko)]) {
          current_genes <- unique(gene2ko[gene2ko[,2] %in% current_ko,1,drop=F])
          current_genes <- all_genes_sample[as.vector(current_genes$gene), sample, drop=F]
          current_genes <- current_genes[sort(current_genes[,1], index=T, dec=T)$ix,,drop=F] # this order genes by abundance
          # ---> this order genes by taxa, good luck understanding :)
          gene_order <- TAXA[rownames(current_genes),,drop=F]
          gene_order <- gene_order[ order(gene_order[,3], gene_order[,7], gene_order[,8]), ,drop=F]
          gene_order <- gene_order[intersect(rownames(gene_order), rownames(current_genes)),,drop=F]
          genes_taxa <- gene_order
          genes_taxa[,1] <- as.vector(genes_taxa[,1])
          genes_taxa[,2] <- as.vector(genes_taxa[,2])
          genes_taxa[,3] <- as.vector(genes_taxa[,3])
          genes_taxa[,4] <- as.vector(genes_taxa[,4])
          genes_taxa[,5] <- as.vector(genes_taxa[,5])
          genes_taxa[,6] <- as.vector(genes_taxa[,6])
          genes_taxa[,7] <- as.vector(genes_taxa[,7])
          genes_taxa[,8] <- as.vector(genes_taxa[,8])
          genes_taxa[,9] <- as.vector(genes_taxa[,9])
          genes_taxa[,10] <- as.vector(genes_taxa[,10])
          genes_taxa[is.na(genes_taxa)] <- 'Unknown'
          gene_order <- c(rownames(gene_order),
                          rownames(current_genes[setdiff(rownames(current_genes), rownames(gene_order)),,drop=F]))
          current_genes <- current_genes[gene_order,,drop=F]
          # <---      
          KO_START <- X_START
          X_RANGE <- 1 + MAX_GENES * GENE_WIDTH
          X_STOP <- X_START + X_RANGE
          X_SEQUENCE <- seq(from = X_START+GENE_WIDTH/2, to = X_STOP-GENE_WIDTH/2, by = GENE_WIDTH)
          ncols <- length(X_SEQUENCE)
          nrows <- ceiling(dim(current_genes)[1]/ncols)
          X_POSITIONS <- rep(X_SEQUENCE, nrows)[1:dim(current_genes)[1]]
          Y_POSITIONS <- NULL
          X_START <- X_STOP + SPACE_BW_KOS # time to set the new start
          for (c in nrows:1) {
            Y_POSITIONS <- c(Y_POSITIONS, rep((GENE_HEIGHT+SPACE_BW_GENE_ROWS/2)*c, ncols))
          }
          Y_POSITIONS <- Y_POSITIONS[1:dim(current_genes)[1]]
          Y_POSITIONS <- Y_POSITIONS - max(Y_POSITIONS)
          
          # fill the data frame
          V2 <- (rownames(current_genes))
          CCG <- current_genes[,1]
          CCGH <- all_genes_horizontal[rownames(current_genes),sample]
          if(length(rownames(current_genes))==0) {
            V2 <- 'NA'
          }
          if(dim(current_genes)[1]==0) {
            CCG <- 0
          }
          if(length(all_genes_horizontal[rownames(current_genes),sample])==0) {
            CCGH <- 0 
          }
          
          to_plot_coverage_gene <- as.data.frame(
            cbind(sample,
                  V2, # names
                  current_ko, 
                  current_module,
                  CCG, # coverage here has been rescaled above
                  GENE_WIDTH,
                  1, # this is alpha, used for later
                  X_POSITIONS, 
                  Y_POSITIONS, # this shifts the heights to a common Y axis
                  'Unknown', # these are the three taxonomic annotations
                  'Unknown',
                  'Unknown',
                  'GENE',
                  'Gene coverages'
            ))
          to_plot_hor_coverage_gene <- as.data.frame(
            cbind(sample,
                  V2, # names
                  current_ko, 
                  current_module,
                  CCGH, # horizontal coverage
                  GENE_WIDTH, 
                  1, # this is alpha, used for later
                  X_POSITIONS, 
                  Y_POSITIONS, 
                  'Unknown', # these are the three taxonomic annotations
                  'Unknown',
                  'Unknown',
                  'GENE',
                  'Horizontal gene coverages'
            ))
          to_plot_hor_coverage_gene[,c(10:13)] <-
            as.matrix(to_plot_hor_coverage_gene[,c(10:13)])
          to_plot_coverage_gene[,c(10:13)] <-
            as.matrix(to_plot_coverage_gene[,c(10:13)])
          rownames(to_plot_hor_coverage_gene) <- to_plot_hor_coverage_gene$V2
          rownames(to_plot_coverage_gene) <- to_plot_coverage_gene$V2
          to_plot_hor_coverage_gene[rownames(genes_taxa),10] <- as.vector(genes_taxa[,3])
          to_plot_hor_coverage_gene[rownames(genes_taxa),11] <- as.vector(genes_taxa[,7])
          to_plot_hor_coverage_gene[rownames(genes_taxa),12] <- as.vector(genes_taxa[,8])
          to_plot_coverage_gene[rownames(genes_taxa),10] <- as.vector(genes_taxa[,3])
          to_plot_coverage_gene[rownames(genes_taxa),11] <- as.vector(genes_taxa[,7])
          to_plot_coverage_gene[rownames(genes_taxa),12] <- as.vector(genes_taxa[,8])
          to_plot_coverage_gene_ALL <- rbind(to_plot_coverage_gene_ALL, to_plot_coverage_gene)
          to_plot_coverage_gene_current_mod <- rbind(to_plot_coverage_gene_current_mod, to_plot_coverage_gene)
          to_plot_coverage_gene_current_pathway <- rbind(to_plot_coverage_gene_current_pathway, to_plot_coverage_gene)
          to_plot_hor_coverage_gene_ALL <- rbind(to_plot_hor_coverage_gene_ALL, to_plot_hor_coverage_gene)
          
          # time to add the KO
          to_plot_coverage_ko <- as.data.frame(
            cbind(sample,
                  'not applicable', # names
                  current_ko, 
                  current_module,
                  ko[current_ko,sample]*KO_HEIGHT/max(ko), # horizontal coverage
                  X_RANGE, 
                  1, # this is alpha, used for later
                  as.numeric(as.vector( KO_START+X_RANGE/2 )), 
                  KO_ABOVE_GENE, 
                  'Unknown', # these are the three taxonomic annotations
                  'Unknown',
                  'Unknown',
                  current_ko,
                  'KO coverages'
            ))
          to_plot_coverage_ko[,c(10:13)] <-
            as.matrix(to_plot_coverage_ko[,c(10:13)])
          to_plot_coverage_ko_ALL <- rbind(to_plot_coverage_ko_ALL, to_plot_coverage_ko)
          # HORIZONTAL
          to_plot_hor_coverage_ko <- to_plot_coverage_ko # HOR
          to_plot_hor_coverage_ko[,14] <- 'Horizontal KO coverages'
          to_plot_hor_coverage_ko[,5] <- Hko[current_ko,sample]*KO_HEIGHT/max(Hko) # HOR
          to_plot_hor_coverage_ko_ALL <- rbind(to_plot_hor_coverage_ko_ALL, to_plot_hor_coverage_ko) # HOR
          
          # add the taxonomic breakdown for each KO
          #print (dim(to_plot_coverage_gene))
          
          #if(dim(to_plot_coverage_gene)[1]>0) {}
          
          x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene[,5])), by=list(as.vector(to_plot_coverage_gene[,10])), FUN=sum, na.rm=T )
          y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; P <- y[pick,]
          P <- rbind(P, x[x[,1]=='Unknown',] )
          P$x <- (P$x*1/sum(P$x) * X_RANGE)
          P$x_pos <- cumsum(P$x)-P$x/2 + KO_START
          x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene[,5])), by=list(as.vector(to_plot_coverage_gene[,11])), FUN=sum, na.rm=T )
          y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; G <- y[pick,]
          G <- rbind(G, x[x[,1]=='Unknown',] )
          G$x <- (G$x*1/sum(G$x) * X_RANGE)
          G$x_pos <- cumsum(G$x)-G$x/2 + KO_START 
          x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene[,5])), by=list(as.vector(to_plot_coverage_gene[,12])), FUN=sum, na.rm=T )
          y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; S <- y[pick,]
          S <- rbind(S, x[x[,1]=='Unknown',] )
          S$x <- (S$x*1/sum(S$x) * X_RANGE)
          S$x_pos <- cumsum(S$x)-S$x/2 + KO_START 
          S2 <- S
          if(!print_species){
            S2[S2[,1] != 'Unknown',1] <- 'ANNOTATED SPECIES' # fix to remove species colors
          }
          l1 <- 2.5
          l2 <- 1.5
          l3 <- 0.5
          if (sum(selected_types %in% 'KO species abundances') == 1 &
                sum(selected_types %in% 'KO genus abundances') == 1 &
                sum(selected_types %in% 'KO phylum abundances') == 1) {
            l1 <- 2.5
            l2 <- 1.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 0 &
                sum(selected_types %in% 'KO genus abundances') == 1 &
                sum(selected_types %in% 'KO phylum abundances') == 1) {
            l1 <- 2.5
            l2 <- 1.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 1 &
                sum(selected_types %in% 'KO genus abundances') == 0 &
                sum(selected_types %in% 'KO phylum abundances') == 1) {
            l1 <- 1.5
            l2 <- 1.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 1 &
                sum(selected_types %in% 'KO genus abundances') == 1 &
                sum(selected_types %in% 'KO phylum abundances') == 0) {
            l1 <- 1.5
            l2 <- 0.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 1 &
                sum(selected_types %in% 'KO genus abundances') == 0 &
                sum(selected_types %in% 'KO phylum abundances') == 0) {
            l1 <- 0.5
            l2 <- 0.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 0 &
                sum(selected_types %in% 'KO genus abundances') == 1 &
                sum(selected_types %in% 'KO phylum abundances') == 0) {
            l1 <- 0.5
            l2 <- 0.5
            l3 <- 0.5
          }
          if (sum(selected_types %in% 'KO species abundances') == 0 &
                sum(selected_types %in% 'KO genus abundances') == 0 &
                sum(selected_types %in% 'KO phylum abundances') == 1) {
            l1 <- 0.5
            l2 <- 0.5
            l3 <- 0.5
          }
          to_plot_ko_s <- as.data.frame(
            cbind(sample,
                  'not applicable', # names
                  current_ko, 
                  current_module,
                  HEIGHT_KO_TAXA, # horizontal coverage
                  S$x, # width
                  1, # this is alpha, used for later
                  S$x_pos, # X position
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR, # this is the Y position
                  'Unset', # these are the three taxonomic annotations
                  'Unset',
                  S[,1],
                  S2[,1], # this is the color
                  'KO species abundances',
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR # Y_original
            ))
          to_plot_ko_s_ALL <- rbind(to_plot_ko_s_ALL, to_plot_ko_s)
          to_plot_ko_g <- as.data.frame(
            cbind(sample,
                  'not applicable', # names
                  current_ko, 
                  current_module,
                  HEIGHT_KO_TAXA, # horizontal coverage
                  G$x, # width
                  1, # this is alpha, used for later
                  G$x_pos, # X position
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR, # this is the Y position
                  'Unset', # these are the three taxonomic annotations
                  G[,1],
                  'Unset',
                  G[,1], # this is the color
                  'KO genus abundances',
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR # Y_original
            ))
          to_plot_ko_g_ALL <- rbind(to_plot_ko_g_ALL, to_plot_ko_g)    
          to_plot_ko_p <- as.data.frame(
            cbind(sample,
                  'not applicable', # names
                  current_ko, 
                  current_module,
                  HEIGHT_KO_TAXA, # horizontal coverage
                  P$x, # width
                  1, # this is alpha, used for later
                  P$x_pos, # X position
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR, # this is the Y position
                  P[,1], # these are the three taxonomic annotations
                  'Unset',
                  'Unset',
                  P[,1], # this is the color
                  'KO phylum abundances',
                  KO_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR # Y_original
            ))
          to_plot_ko_p_ALL <- rbind(to_plot_ko_p_ALL, to_plot_ko_p)
          
          
          
        } # end loop over KO
        
        
        # MOD ################################################
        # time to add the MOD
        to_plot_coverage_mod <- as.data.frame(
          cbind(sample,
                'not applicable', # names
                current_ko, 
                current_module,
                MODULE_SUM, # horizontal coverage
                X_START-MOD_START+1, 
                1, # this is alpha, used for later
                as.numeric(as.vector( MOD_START+(X_START-MOD_START+1)/2 )), 
                MOD_ABOVE_GENE, 
                'Unknown', # these are the three taxonomic annotations
                'Unknown',
                'Unknown',
                current_module,
                'Module coverages'
          ))
        to_plot_coverage_mod[,c(10:13)] <-
          as.matrix(to_plot_coverage_mod[,c(10:13)])
        to_plot_coverage_mod_ALL <- rbind(to_plot_coverage_mod_ALL, to_plot_coverage_mod)
        # HORIZONTAL
        to_plot_hor_coverage_mod <- to_plot_coverage_mod # HOR
        to_plot_hor_coverage_mod[,14] <- 'Horizontal module coverages'
        to_plot_hor_coverage_mod[,5] <- MODULE_SUM_HOR # HOR
        to_plot_hor_coverage_mod_ALL <- rbind(to_plot_hor_coverage_mod_ALL, to_plot_hor_coverage_mod) # HOR
        
        # add the taxonomic breakdown for each MOD
        MOD_RANGE <- X_START-MOD_START+1
        x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_mod[,5])), by=list(as.vector(to_plot_coverage_gene_current_mod[,10])), FUN=sum, na.rm=T )
        y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; P <- y[pick,]
        P <- rbind(P, x[x[,1]=='Unknown',] )
        P$x <- (P$x*1/sum(P$x) * MOD_RANGE)
        P$x_pos <- cumsum(P$x)-P$x/2 + MOD_START
        x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_mod[,5])), by=list(as.vector(to_plot_coverage_gene_current_mod[,11])), FUN=sum, na.rm=T )
        y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; G <- y[pick,]
        G <- rbind(G, x[x[,1]=='Unknown',] )
        G$x <- (G$x*1/sum(G$x) * MOD_RANGE)
        G$x_pos <- cumsum(G$x)-G$x/2 + MOD_START 
        # fix to remove species colors
        x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_mod[,5])), by=list(as.vector(to_plot_coverage_gene_current_mod[,12])), FUN=sum, na.rm=T )
        y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; S <- y[pick,]
        S <- rbind(S, x[x[,1]=='Unknown',] )
        S$x <- (S$x*1/sum(S$x) * MOD_RANGE)
        S$x_pos <- cumsum(S$x)-S$x/2 + MOD_START 
        S2 <- S
        if(!print_species){
          S2[S2[,1] != 'Unknown',1] <- 'ANNOTATED SPECIES'
        }
        l1 <- 2.5
        l2 <- 1.5
        l3 <- 0.5
        if (sum(selected_types %in% 'Module species abundances') == 1 &
              sum(selected_types %in% 'Module genus abundances') == 1 &
              sum(selected_types %in% 'Module phylum abundances') == 1) {
          l1 <- 2.5
          l2 <- 1.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 0 &
              sum(selected_types %in% 'Module genus abundances') == 1 &
              sum(selected_types %in% 'Module phylum abundances') == 1) {
          l1 <- 2.5
          l2 <- 1.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 1 &
              sum(selected_types %in% 'Module genus abundances') == 0 &
              sum(selected_types %in% 'Module phylum abundances') == 1) {
          l1 <- 1.5
          l2 <- 1.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 1 &
              sum(selected_types %in% 'Module genus abundances') == 1 &
              sum(selected_types %in% 'Module phylum abundances') == 0) {
          l1 <- 1.5
          l2 <- 0.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 1 &
              sum(selected_types %in% 'Module genus abundances') == 0 &
              sum(selected_types %in% 'Module phylum abundances') == 0) {
          l1 <- 0.5
          l2 <- 0.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 0 &
              sum(selected_types %in% 'Module genus abundances') == 1 &
              sum(selected_types %in% 'Module phylum abundances') == 0) {
          l1 <- 0.5
          l2 <- 0.5
          l3 <- 0.5
        }
        if (sum(selected_types %in% 'Module species abundances') == 0 &
              sum(selected_types %in% 'Module genus abundances') == 0 &
              sum(selected_types %in% 'Module phylum abundances') == 1) {
          l1 <- 0.5
          l2 <- 0.5
          l3 <- 0.5
        }
        to_plot_mod_s <- as.data.frame(
          cbind(sample,
                'not applicable', # names
                current_ko, 
                current_module,
                HEIGHT_KO_TAXA, # horizontal coverage
                S$x, # width
                1, # this is alpha, used for later
                S$x_pos, # X position
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR, # this is the Y position
                'Unset', # these are the three taxonomic annotations
                'Unset',
                S[,1],
                S2[,1], # this is the color
                'Module species abundances',
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR # Y_original
          ))
        to_plot_mod_s_ALL <- rbind(to_plot_mod_s_ALL, to_plot_mod_s)
        to_plot_mod_g <- as.data.frame(
          cbind(sample,
                'not applicable', # names
                current_ko, 
                current_module,
                HEIGHT_KO_TAXA, # horizontal coverage
                G$x, # width
                1, # this is alpha, used for later
                G$x_pos, # X position
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR, # this is the Y position
                'Unset', # these are the three taxonomic annotations
                G[,1],
                'Unset',
                G[,1], # this is the color
                'Module genus abundances',
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR # Y_original
          ))
        to_plot_mod_g_ALL <- rbind(to_plot_mod_g_ALL, to_plot_mod_g)    
        to_plot_mod_p <- as.data.frame(
          cbind(sample,
                'not applicable', # names
                current_ko, 
                current_module,
                HEIGHT_KO_TAXA, # horizontal coverage
                P$x, # width
                1, # this is alpha, used for later
                P$x_pos, # X position
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR, # this is the Y position
                P[,1], # these are the three taxonomic annotations
                'Unset',
                'Unset',
                P[,1], # this is the color
                'Module phylum abundances',
                MOD_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR # Y_original
          ))
        to_plot_mod_p_ALL <- rbind(to_plot_mod_p_ALL, to_plot_mod_p)
        # MOD ################################################
        
        X_END_PATHWAY <- X_START
        X_START <- X_START + SPACE_BW_MODULES
      } # end loop over module





# PATHWAY ################################################
# time to add the PATHWAY
to_plot_coverage_pathway <- as.data.frame(
  cbind(sample,
        'not applicable', # names
        current_ko, 
        current_pathway,
        PATHWAY_SUM, # horizontal coverage
        X_END_PATHWAY-PATHWAY_START+1, 
        1, # this is alpha, used for later
        as.numeric(as.vector( PATHWAY_START+(X_END_PATHWAY-PATHWAY_START+1)/2 )), 
        PATHWAY_ABOVE_GENE, 
        'Unknown', # these are the three taxonomic annotations
        'Unknown',
        'Unknown',
        current_pathway,
        'Pathway coverages'
  ))
to_plot_coverage_pathway[,c(10:13)] <-
  as.matrix(to_plot_coverage_pathway[,c(10:13)])
to_plot_coverage_pathway_ALL <- rbind(to_plot_coverage_pathway_ALL, to_plot_coverage_pathway)
# HORIZONTAL
to_plot_hor_coverage_pathway <- to_plot_coverage_pathway # HOR
to_plot_hor_coverage_pathway[,14] <- 'Horizontal pathway coverages'
to_plot_hor_coverage_pathway[,5] <- PATHWAY_SUM_HOR # HOR
to_plot_hor_coverage_pathway_ALL <- rbind(to_plot_hor_coverage_pathway_ALL, to_plot_hor_coverage_pathway) # HOR

# add the taxonomic breakdown for each PATHWAY
PATHWAY_RANGE <- X_END_PATHWAY-PATHWAY_START+1
x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_pathway[,5])), by=list(as.vector(to_plot_coverage_gene_current_pathway[,10])), FUN=sum, na.rm=T )
y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; P <- y[pick,]
P <- rbind(P, x[x[,1]=='Unknown',] )
P$x <- (P$x*1/sum(P$x) * PATHWAY_RANGE)
P$x_pos <- cumsum(P$x)-P$x/2 + PATHWAY_START
x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_pathway[,5])), by=list(as.vector(to_plot_coverage_gene_current_pathway[,11])), FUN=sum, na.rm=T )
y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; G <- y[pick,]
G <- rbind(G, x[x[,1]=='Unknown',] )
G$x <- (G$x*1/sum(G$x) * PATHWAY_RANGE)
G$x_pos <- cumsum(G$x)-G$x/2 + PATHWAY_START 
# fix to remove species colors
x <- aggregate(as.numeric(as.vector(to_plot_coverage_gene_current_pathway[,5])), by=list(as.vector(to_plot_coverage_gene_current_pathway[,12])), FUN=sum, na.rm=T )
y <- x[x[,1] != 'Unknown',]; pick <- sort(y[, 2], dec=T, index=T)$ix; S <- y[pick,]
S <- rbind(S, x[x[,1]=='Unknown',] )
S$x <- (S$x*1/sum(S$x) * PATHWAY_RANGE)
S$x_pos <- cumsum(S$x)-S$x/2 + PATHWAY_START 
S2 <- S
if(!print_species){
  S2[S2[,1] != 'Unknown',1] <- 'ANNOTATED SPECIES'
}
l1 <- 2.5
l2 <- 1.5
l3 <- 0.5
if (sum(selected_types %in% 'Pathway species abundances') == 1 &
      sum(selected_types %in% 'Pathway genus abundances') == 1 &
      sum(selected_types %in% 'Pathway phylum abundances') == 1) {
  l1 <- 2.5
  l2 <- 1.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 0 &
      sum(selected_types %in% 'Pathway genus abundances') == 1 &
      sum(selected_types %in% 'Pathway phylum abundances') == 1) {
  l1 <- 2.5
  l2 <- 1.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 1 &
      sum(selected_types %in% 'Pathway genus abundances') == 0 &
      sum(selected_types %in% 'Pathway phylum abundances') == 1) {
  l1 <- 1.5
  l2 <- 1.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 1 &
      sum(selected_types %in% 'Pathway genus abundances') == 1 &
      sum(selected_types %in% 'Pathway phylum abundances') == 0) {
  l1 <- 1.5
  l2 <- 0.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 1 &
      sum(selected_types %in% 'Pathway genus abundances') == 0 &
      sum(selected_types %in% 'Pathway phylum abundances') == 0) {
  l1 <- 0.5
  l2 <- 0.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 0 &
      sum(selected_types %in% 'Pathway genus abundances') == 1 &
      sum(selected_types %in% 'Pathway phylum abundances') == 0) {
  l1 <- 0.5
  l2 <- 0.5
  l3 <- 0.5
}
if (sum(selected_types %in% 'Pathway species abundances') == 0 &
      sum(selected_types %in% 'Pathway genus abundances') == 0 &
      sum(selected_types %in% 'Pathway phylum abundances') == 1) {
  l1 <- 0.5
  l2 <- 0.5
  l3 <- 0.5
}
to_plot_pathway_s <- as.data.frame(
  cbind(sample,
        'not applicable', # names
        current_ko, 
        current_pathway,
        HEIGHT_KO_TAXA, # horizontal coverage
        S$x, # width
        1, # this is alpha, used for later
        S$x_pos, # X position
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR, # this is the Y position
        'Unset', # these are the three taxonomic annotations
        'Unset',
        S[,1],
        S2[,1], # this is the color
        'Pathway species abundances',
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l1-SPACE_BW_KO_COV_HOR # Y_original
  ))
to_plot_pathway_s_ALL <- rbind(to_plot_pathway_s_ALL, to_plot_pathway_s)
to_plot_pathway_g <- as.data.frame(
  cbind(sample,
        'not applicable', # names
        current_ko, 
        current_pathway,
        HEIGHT_KO_TAXA, # horizontal coverage
        G$x, # width
        1, # this is alpha, used for later
        G$x_pos, # X position
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR, # this is the Y position
        'Unset', # these are the three taxonomic annotations
        G[,1],
        'Unset',
        G[,1], # this is the color
        'Pathway genus abundances',
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l2-SPACE_BW_KO_COV_HOR # Y_original
  ))
to_plot_pathway_g_ALL <- rbind(to_plot_pathway_g_ALL, to_plot_pathway_g)    
to_plot_pathway_p <- as.data.frame(
  cbind(sample,
        'not applicable', # names
        current_ko, 
        current_pathway,
        HEIGHT_KO_TAXA, # horizontal coverage
        P$x, # width
        1, # this is alpha, used for later
        P$x_pos, # X position
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR, # this is the Y position
        P[,1], # these are the three taxonomic annotations
        'Unset',
        'Unset',
        P[,1], # this is the color
        'Pathway phylum abundances',
        PATHWAY_ABOVE_GENE-HEIGHT_KO_TAXA*l3-SPACE_BW_KO_COV_HOR # Y_original
  ))
to_plot_pathway_p_ALL <- rbind(to_plot_pathway_p_ALL, to_plot_pathway_p)
# PATHWAY ################################################



X_START <- X_START + SPACE_BW_PATHWAYS
    } # end pathways


# we need to shift the height up or down FOR GENES
colnames(to_plot_coverage_gene_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_coverage_gene_ALL$Y <- as.numeric(as.vector(to_plot_coverage_gene_ALL$Y))
to_plot_coverage_gene_ALL$height <- as.numeric(as.vector(to_plot_coverage_gene_ALL$height))
to_plot_coverage_gene_ALL$height <- to_plot_coverage_gene_ALL$height * GENE_HEIGHT / (max(to_plot_coverage_gene_ALL$height))
to_plot_coverage_gene_ALL$Y_original <- to_plot_coverage_gene_ALL$Y
to_plot_coverage_gene_ALL$Y <- to_plot_coverage_gene_ALL$Y + to_plot_coverage_gene_ALL$height/2
# we need to shift the height up or down FOR GENES HORIZONTAL
colnames(to_plot_hor_coverage_gene_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_hor_coverage_gene_ALL$Y <- as.numeric(as.vector(to_plot_hor_coverage_gene_ALL$Y))
to_plot_hor_coverage_gene_ALL$height <- as.numeric(as.vector(to_plot_hor_coverage_gene_ALL$height))
to_plot_hor_coverage_gene_ALL$height <- to_plot_hor_coverage_gene_ALL$height * GENE_HEIGHT / (max(to_plot_hor_coverage_gene_ALL$height))
to_plot_hor_coverage_gene_ALL$Y_original <- to_plot_hor_coverage_gene_ALL$Y
move_it <- sum(selected_types %in% c('Gene species abundances', 'Gene genus abundances', 'Gene phylum abundances'))
to_plot_hor_coverage_gene_ALL$Y <- to_plot_hor_coverage_gene_ALL$Y - to_plot_hor_coverage_gene_ALL$height/2 - SPACE_FOR_EACH_TAXA/3*move_it

# time to add the taxa 
l1 <- 1.2
l2 <- 2
l3 <- 6
if (sum(selected_types %in% 'Gene species abundances') == 1 &
      sum(selected_types %in% 'Gene genus abundances') == 1 &
      sum(selected_types %in% 'Gene phylum abundances') == 1) {
  l1 <- 1.2
  l2 <- 2
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 0 &
      sum(selected_types %in% 'Gene genus abundances') == 1 &
      sum(selected_types %in% 'Gene phylum abundances') == 1) {
  l1 <- 1.2
  l2 <- 2
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 1 &
      sum(selected_types %in% 'Gene genus abundances') == 0 &
      sum(selected_types %in% 'Gene phylum abundances') == 1) {
  l1 <- 2
  l2 <- 2
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 1 &
      sum(selected_types %in% 'Gene genus abundances') == 1 &
      sum(selected_types %in% 'Gene phylum abundances') == 0) {
  l1 <- 2
  l2 <- 6
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 1 &
      sum(selected_types %in% 'Gene genus abundances') == 0 &
      sum(selected_types %in% 'Gene phylum abundances') == 0) {
  l1 <- 6
  l2 <- 6
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 0 &
      sum(selected_types %in% 'Gene genus abundances') == 1 &
      sum(selected_types %in% 'Gene phylum abundances') == 0) {
  l1 <- 6
  l2 <- 6
  l3 <- 6
}
if (sum(selected_types %in% 'Gene species abundances') == 0 &
      sum(selected_types %in% 'Gene genus abundances') == 0 &
      sum(selected_types %in% 'Gene phylum abundances') == 1) {
  l1 <- 6
  l2 <- 6
  l3 <- 6
}
to_plot_phylum_ALL <- to_plot_coverage_gene_ALL
to_plot_phylum_ALL[,14] <- 'Gene phylum abundances'
to_plot_phylum_ALL$Color <- to_plot_phylum_ALL$Phylum
to_plot_phylum_ALL$height <- SPACE_FOR_EACH_TAXA/3
to_plot_phylum_ALL$Y <- to_plot_phylum_ALL$Y_original - SPACE_FOR_EACH_TAXA/l3
to_plot_genus_ALL <- to_plot_coverage_gene_ALL
to_plot_genus_ALL[,14] <- 'Gene genus abundances'
to_plot_genus_ALL$Color <- to_plot_genus_ALL$Genus
to_plot_genus_ALL$height <- SPACE_FOR_EACH_TAXA/3
to_plot_genus_ALL$Y <- to_plot_genus_ALL$Y_original - SPACE_FOR_EACH_TAXA/l2
to_plot_species_ALL <- to_plot_coverage_gene_ALL
to_plot_species_ALL[,14] <- 'Gene species abundances'
if(!print_species){
  to_plot_species_ALL[to_plot_species_ALL$Species != 'Unknown', 'Color'] <- 'ANNOTATED SPECIES' #to_plot_species_ALL$Species
}
to_plot_species_ALL[to_plot_species_ALL$Species == 'Unknown', 'Color'] <- 'Unknown'
to_plot_species_ALL$height <- SPACE_FOR_EACH_TAXA/3
to_plot_species_ALL$Y <- to_plot_species_ALL$Y_original - SPACE_FOR_EACH_TAXA/l1

# time to add the KOs
colnames(to_plot_coverage_ko_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_coverage_ko_ALL$Y <- as.numeric(as.vector(to_plot_coverage_ko_ALL$Y))
to_plot_coverage_ko_ALL$height <- as.numeric(as.vector(to_plot_coverage_ko_ALL$height))
to_plot_coverage_ko_ALL$Y_original <- to_plot_coverage_ko_ALL$Y
to_plot_coverage_ko_ALL$Y <- to_plot_coverage_ko_ALL$Y + to_plot_coverage_ko_ALL$height/2
# HOR
colnames(to_plot_hor_coverage_ko_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_hor_coverage_ko_ALL$Y <- as.numeric(as.vector(to_plot_hor_coverage_ko_ALL$Y))
to_plot_hor_coverage_ko_ALL$height <- as.numeric(as.vector(to_plot_hor_coverage_ko_ALL$height))
to_plot_hor_coverage_ko_ALL$Y_original <- to_plot_hor_coverage_ko_ALL$Y
move_it <- sum(selected_types %in% c('KO species abundances', 'KO genus abundances', 'KO phylum abundances'))
to_plot_hor_coverage_ko_ALL$Y <- to_plot_hor_coverage_ko_ALL$Y - to_plot_hor_coverage_ko_ALL$height/2 - 2*SPACE_BW_KO_COV_HOR - 3*move_it
# TAXA for KO
colnames(to_plot_ko_s_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_ko_p_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_ko_g_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
# TAXA for MOD
colnames(to_plot_mod_s_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_mod_p_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_mod_g_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
# add the MOD
colnames(to_plot_coverage_mod_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_coverage_mod_ALL$Y <- as.numeric(as.vector(to_plot_coverage_mod_ALL$Y))
to_plot_coverage_mod_ALL$height <- as.numeric(as.vector(to_plot_coverage_mod_ALL$height))
to_plot_coverage_mod_ALL$Y_original <- to_plot_coverage_mod_ALL$Y
to_plot_coverage_mod_ALL$Y <- to_plot_coverage_mod_ALL$Y + to_plot_coverage_mod_ALL$height/2
# HOR
colnames(to_plot_hor_coverage_mod_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_hor_coverage_mod_ALL$Y <- as.numeric(as.vector(to_plot_hor_coverage_mod_ALL$Y))
to_plot_hor_coverage_mod_ALL$height <- as.numeric(as.vector(to_plot_hor_coverage_mod_ALL$height))
to_plot_hor_coverage_mod_ALL$Y_original <- to_plot_hor_coverage_mod_ALL$Y
move_it <- sum(selected_types %in% c('Module species abundances', 'Module genus abundances', 'Module phylum abundances'))
to_plot_hor_coverage_mod_ALL$Y <- to_plot_hor_coverage_mod_ALL$Y - to_plot_hor_coverage_mod_ALL$height/2 - 2*SPACE_BW_KO_COV_HOR - 3*move_it

# TAXA for PATHWAY
colnames(to_plot_pathway_s_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_pathway_p_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
colnames(to_plot_pathway_g_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
# add the MOD
colnames(to_plot_coverage_pathway_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_coverage_pathway_ALL$Y <- as.numeric(as.vector(to_plot_coverage_pathway_ALL$Y))
to_plot_coverage_pathway_ALL$height <- as.numeric(as.vector(to_plot_coverage_pathway_ALL$height))
to_plot_coverage_pathway_ALL$Y_original <- to_plot_coverage_pathway_ALL$Y
to_plot_coverage_pathway_ALL$Y <- to_plot_coverage_pathway_ALL$Y + to_plot_coverage_pathway_ALL$height/2
# HOR
colnames(to_plot_hor_coverage_pathway_ALL) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type')
to_plot_hor_coverage_pathway_ALL$Y <- as.numeric(as.vector(to_plot_hor_coverage_pathway_ALL$Y))
to_plot_hor_coverage_pathway_ALL$height <- as.numeric(as.vector(to_plot_hor_coverage_pathway_ALL$height))
to_plot_hor_coverage_pathway_ALL$Y_original <- to_plot_hor_coverage_pathway_ALL$Y
move_it <- sum(selected_types %in% c('Pathway species abundances', 'Pathway genus abundances', 'Pathway phylum abundances'))
to_plot_hor_coverage_pathway_ALL$Y <- to_plot_hor_coverage_pathway_ALL$Y - to_plot_hor_coverage_pathway_ALL$height/2 - 2*SPACE_BW_KO_COV_HOR - 3*move_it






#to_plot_mod_s_ALL$Color <- factor(to_plot_mod_s_ALL$Color)

TO_PLOT <- as.data.frame(rbind(TO_PLOT, 
                               to_plot_coverage_gene_ALL, to_plot_hor_coverage_gene_ALL,
                               to_plot_phylum_ALL, to_plot_genus_ALL, to_plot_species_ALL,
                               to_plot_coverage_ko_ALL, to_plot_hor_coverage_ko_ALL,
                               to_plot_ko_p_ALL, to_plot_ko_g_ALL, to_plot_ko_s_ALL,
                               to_plot_coverage_mod_ALL, to_plot_hor_coverage_mod_ALL,
                               to_plot_coverage_pathway_ALL, to_plot_hor_coverage_pathway_ALL,
                               to_plot_mod_p_ALL, to_plot_mod_g_ALL, to_plot_mod_s_ALL,
                               to_plot_pathway_p_ALL, to_plot_pathway_g_ALL, to_plot_pathway_s_ALL
))


















X_START <- X_START + SPACE_BW_GROUPS
  } # end groups actually run
} # end groups









colnames(TO_PLOT) <- c('sample', 'gene', 'ko', 'module', 'height', 'width', 'alpha', 'X', 'Y', 'Phylum', 'Genus', 'Species', 'Color', 'Type', 'Y_original')
TO_PLOT$width <- as.numeric(as.vector(TO_PLOT$width))
TO_PLOT$height <- as.numeric(as.vector(TO_PLOT$height))
TO_PLOT$alpha <- as.numeric(as.vector(TO_PLOT$alpha))
TO_PLOT$X <- as.numeric(as.vector(TO_PLOT$X))
TO_PLOT$Y <- as.numeric(as.vector(TO_PLOT$Y))
TO_PLOT$Color <- as.vector(TO_PLOT$Color)

#print(unique(TO_PLOT$Type))
TO_PLOT <- TO_PLOT[TO_PLOT$Type %in% selected_types,]
#print(unique(TO_PLOT$Type))

# let's get the legend colors in a nice way
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}



COLORS <- 'black' # Coverage

#if (sum(grepl('NO PATHWAY',selected_pathways))==1){ # pathways
#  COLORS <- c(COLORS, ggplotColours( length(selected_pathways)-1 ), 'grey')
#} else {
  COLORS <- c(COLORS, ggplotColours( length(selected_pathways)))
#}


#if (sum(grepl('NO MODULE',names(order_module)))==1){ # modules
#  COLORS <- c(COLORS, ggplotColours( length(order_module)-1 ), 'grey')
#} else {
  COLORS <- c(COLORS, ggplotColours( length(selected_modules)))
#}

tmp <- sort(as.vector(unique(TO_PLOT$Phylum)))
pc <- tmp[tmp != 'Unknown']
tmp <- sort(as.vector(unique(TO_PLOT$Genus)))
gc <- tmp[tmp != 'Unknown']
tmp <- sort(as.vector(unique(TO_PLOT$Species)))
sc <- tmp[tmp != 'Unknown']



COLORS <- c(COLORS, ggplotColours( length(selected_kos))) # KOs
COLORS <- c(COLORS, ggplotColours( length(pc) )) # phylum
COLORS <- c(COLORS, ggplotColours( length(gc) )) # genus

if(!print_species){
  COLORS <- c(COLORS, 'dimgray') # species
}
if(print_species){
  COLORS <- c(COLORS, ggplotColours( length(sc) )) # species
}
COLORS <- c(COLORS, 'grey90') # species
if(!print_species){
  COLORS_NAMES <- c('GENE', selected_pathways, selected_modules, selected_kos,
                    as.vector(pc),
                    as.vector(gc),
                    'ANNOTATED SPECIES', 'Unknown')
}
if(print_species){
  COLORS_NAMES <- c('GENE', selected_pathways, selected_modules, selected_kos,
                    as.vector(pc),
                    as.vector(gc),
                    as.vector(sc),
                    'Unknown')
}
names(COLORS) <- COLORS_NAMES
TO_PLOT$Color <- factor(TO_PLOT$Color, levels=COLORS_NAMES[!duplicated(COLORS_NAMES)])      

#save(list=c('kegg2name', 'ko2pathway', 'SAMPLE.GROUPS', 'all_genes_horizontal', 'TAXA', 'gene_map', 'ko2module', 'our_ko2module', 'ko', 'module', 'pathway', 'Hko', 'Hmodule', 'Hpathway', 'all_genes'), file = '~/WILCOX/data.R')
#save(list=c('kegg2name', 'gene2ko', 'ko2pathway', 'SAMPLE.GROUPS', 'all_genes_horizontal', 'TAXA', 'gene_map', 'ko2module', 'our_ko2module', 'ko', 'module', 'pathway', 'Hko', 'Hmodule', 'Hpathway', 'all_genes'), file = '~/WILCOX/data.R')