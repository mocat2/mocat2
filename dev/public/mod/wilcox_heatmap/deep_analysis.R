# The idea is that using the previous scripts in the
# wilcox heatmap package we have identified a number of
# e.g. modules that we'd like to look closer at, and that's
# what this script does.Currently it's semi focused on module selection.

# specify the file to load, it should be a file
# that is the result of step two of this package.
FILE <- '/g/bork1/kultima/projects/m3/testset/maxmin10/TMP/out/2015Jun16_195233/data/session.RData'
load(FILE)

ko2pathway <- read.table('/g/bork/kultima/db/kegg/ko2pathway.txt')
ko2module <- read.table('/g/bork/kultima/db/kegg/ko2module.txt')

# Select a module
module <- 'M00080'

# because this is a module, we can check the associated pathways
# we could either go via the genes that are significant, or do it via the KEGG map
# I prefer the first
file <- SETTINGS['significant_gene', 'setting']
genes <- read.table(file, header=T, quote="", row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char="", fill = T)
selected_genes <- genes[grepl(module, genes$module),] # get genes
pathways <- sort(table(selected_genes$pathway)) # get possible
biggest_pathway <- names(pathways[length(pathways)]) # get bigegst
selected_pathways <- genes[grepl(biggest_pathway, genes$pathway),] # after finding biggest pathway, let's get the other modules of that pathway
tmp <- as.vector(unique(selected_pathways$module))
tmp <- tmp[tmp != '']
all_modules_in_top_pathway <- genes[genes$module %in% tmp,]
tmp <- unique(all_modules_in_top_pathway$module)
unique_modules <- as.vector(tmp[tmp != ''])
unique_modules <- unique(unlist(strsplit(unique_modules,',')))
tmp <- unique(all_modules_in_top_pathway$ko)
unique_kos <- as.vector(tmp[tmp != ''])
unique_kos <- unique(unlist(strsplit(unique_kos,',')))
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

ADDON <- NULL
HADDON <- NULL
ORDER <- NULL
ORDERHA <- NULL
ORDERA <- NULL
for (sample in c(g1, 'CTRL', g2)){
  
  if (sample == "CTRL") {
    ORDERA <- c(ORDERA,p1, '', '', '', '', '', '', '', '', '')
  } else{
  #sample <- 'CCIS31434951ST-20-0'
  #gene_abundance <- gene[rownames(all_modules_in_top_pathway), sample]
  #set <- cbind(all_modules_in_top_pathway, gene_abundance)
  
  tmp <- unique_kos[unique_kos %in% rownames(ko)] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  s_ko <- ko[tmp, sample, drop=F]
  colnames(s_ko) <- paste(sample, 'KO')
  tmp <- unique_modules[unique_modules %in% rownames(module)] # this might have to be improved upon later, but here we get the subset of KOs that are among the singificant ones
  s_module <- module[tmp, sample, drop=F]
  Hs_module <- Hmodule[tmp, sample, drop=F]
  colnames(s_module) <- paste(sample, 'MODULE')
  s_pathway <- pathway[biggest_pathway, sample, drop=F]
  Hs_pathway <- Hpathway[biggest_pathway, sample, drop=F]
  colnames(s_pathway) <- paste(sample, 'PATHWAY')
  pm <- cbind(as.vector(s_pathway), s_module)
  Hpm <- cbind(as.vector(Hs_pathway), Hs_module)
  #colnames(pm)[1] <- colnames(s_pathway)
  #colnames(Hpm)[1] <- colnames(s_pathway)
  addon <- NULL
  Haddon <- NULL
  for (m in rownames(pm)) {
    get_these <- rownames(ko) %in% as.vector(ko2module[ko2module[,2] %in% m,1])
    to_add <- ko[get_these,sample,drop=F] # note that here we use all KOs, and not only the significant ones...
    get_these <- rownames(to_add)
    Hto_add <- Hko[get_these,sample,drop=F] # note that here we use all KOs, and not only the significant ones...
    rownames(to_add) <- c(paste(rownames(to_add), '|', m, '|', biggest_pathway, sep=''))
    rownames(Hto_add) <- c(paste(rownames(Hto_add), '|', m, '|', biggest_pathway, sep=''))
    addon <- rbind(addon, cbind(as.vector(pm[m,1]),  as.vector(pm[m,2]) ,to_add))
    Haddon <- rbind(Haddon, cbind(as.vector(Hpm[m,1]),  as.vector(Hpm[m,2]) ,Hto_add))
  }
  colnames(addon) <- c('C-pathway', 'C-module', 'C-ko')
  rownames(addon) <- paste(sample, rownames(addon))
  colnames(Haddon) <- c('H-pathway', 'H-module', 'H-ko')
  rownames(Haddon) <- paste(sample, rownames(Haddon))
  ADDON <- rbind(ADDON, addon)
  HADDON <- rbind(HADDON, Haddon)
  p1 <- rownames(addon)
  p2 <- rownames(Haddon)
  ORDER <- c(ORDER,p1,p2)
  ORDERA <- c(ORDERA,p1, '', '', '')
  ORDERHA <- c(ORDERHA,p2)
  }
}

mHADDON <- melt(HADDON)
mADDON <- melt(ADDON)
mADDON$value <- -(mADDON$value * (max(mHADDON$value)/max(mADDON$value))) # rescale

YORDER <- rev(c('C-pathway', 'H-pathway', 'C-module', 'H-module', 'C-ko', 'H-ko'))

bm <- rbind(mADDON, mHADDON)

ggplot(bm, aes(X1, X2)) +
  geom_tile(aes(fill = value), colour = NA) +
  scale_fill_gradient2(low ="darkorange4", mid = "white", high = "steelblue4") +
  scale_x_discrete(limits=(ORDERA)) +
  scale_y_discrete(limits=unique(YORDER)) +
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
