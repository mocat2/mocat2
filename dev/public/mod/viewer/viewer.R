args <- commandArgs(trailingOnly = TRUE)

# GENE_ZIP <- '/g/bork1/kultima/projects/m3/testset/maxmin10/PROFILES/gene.profiles/all/all.gene.profiles.screened.screened.adapter.on.hg19.on.779.CC.RefGeneCatalog.padded.1-2.solexaqa.allbest.l45.p95.zip'
# gene_zip <- 'all.gene.profile.screened.screened.adapter.on.hg19.on.779.CC.RefGeneCatalog.padded.1-2.solexaqa.allbest.l45.p95.insert.scaled.gene'
# gene_hor_zip <- 'all.gene.profile.screened.screened.adapter.on.hg19.on.779.CC.RefGeneCatalog.padded.1-2.solexaqa.allbest.l45.p95.horizontal.gene'
# KEGG_ZIP <- '/g/bork1/kultima/projects/m3/testset/maxmin10/PROFILES/functional.profiles/all/KEGG/all.KEGG.profiles.screened.screened.adapter.on.hg19.on.779.CC.RefGeneCatalog.padded.1-2.solexaqa.allbest.l45.p95.zip'
# kegg_zip <- 'all.functional.profile.screened.screened.adapter.on.hg19.on.779.CC.RefGeneCatalog.padded.1-2.solexaqa.allbest.l45.p95'
# kegg_type <- 'insert.scaled'
# SAVE <- '/g/bork1/kultima/projects/m3/testset/out'
# MAP <- '/g/bork1/kultima/MOCAT/data/779.CC.RefGeneCatalog.padded.1-2.functional.map.gz'
# pME <- '/g/bork5/mocat/RESOURCES/annotations/kegg/kegg2name.txt'
# TAXONOMIC <- '/g/bork1/kultima/projects/m3/gene_catalogs/779.CC.RefGeneCatalog.faa.ko00540.AGAINST.freeze11.blastp.besthit.v4'
# PATHWAYS <- c('ko00010', 'ko00540')

GENE_ZIP <- args[1]
gene_zip <- args[2]
gene_hor_zip <- args[3]
KEGG_ZIP <- args[4]
kegg_zip <- args[5]
kegg_type <- args[6]
SAVE <- args[7]
MAP <- args[8]
pME <- args[9]
TAXONOMIC <- args[10]
PATHWAYS <- strsplit(args[11], " ")[[1]]
KEGG2NAME <- args[12]
KO2MODULE <- args[13]
KO2PATHWAY <- args[14]


require(data.table)

# The data file already exists, we load it
if(file.exists(paste(SAVE, '.all.RData', sep=''))) {
  cat('Loading existing RData file\n')
  load(paste(SAVE, '.all.RData', sep=''))
} else {
  

  cat ('Processing ko2pathway\n')
  # pathway
  #   map <- unique(cbind(as.vector(tmp$ko), as.vector(tmp$pathway)))
  #   p1 <- strsplit(map[,1], ',')
  #   p2 <- strsplit(map[,2], ',')
  #   ko2pathwayList <- list()
  #   for (list in c(1:length(p1))) {
  #     ko2pathwayList[[list]] <- expand.grid(p1[[list]], p2[[list]] )
  #   }
  #   ko2pathway <- as.data.frame(rbindlist (ko2pathwayList))
  ko2pathway <- read.table(KO2PATHWAY, sep='\t', header = F)
  colnames(ko2pathway) <- c('ko', 'pathway')
  #ko2pathway$ko <- as.vector(ko2pathway$ko)
  #ko2pathway$pathway <- as.vector(ko2pathway$pathway)
  if (KEGG2NAME !="0") {
    cat ('Processing kegg2name\n')
    kegg2name <-read.table(KEGG2NAME, sep='\t', comment.char="", quote = "")
  }
  
  cat ('Processing ko2module\n')
  ko2module <- read.table(KO2MODULE, sep='\t', header = F)
  colnames(ko2module) <- c('ko', 'module')
  
  
  # Gene coverages and horizontal gene coverages
  cat ('Loading gene coverages (this will take a few minutes)\n')
  
  # we can't do this because we need to normalize the gene counts... all_genes <- read.table(pipe(  paste( 'echo "',PATHWAYS,'" | sed "s/,/\\n/g" | sed "s/ //g" > tmp1 && gunzip -c ', MAP, ' | fgrep -f tmp1 - | cut -f 1 > tmp && unzip -p ', GENE_ZIP, ' ', gene_zip , ' | fgrep -f tmp -', sep=' ') ), sep="\t", header=T, comment.char="#", row.names=1)
  all_genes <- read.table(unz(description = GENE_ZIP, filename = gene_zip), sep="\t", header=T, comment.char="#", row.names=1, check.names = F, quote = "")
  cat ('Loading horizontal gene coverages (this will take a few minutes)\n')
  # we can't do this because we need to normalize the gene counts... all_genes_horizontal <- read.table(pipe(  paste( 'echo "',PATHWAYS,'" | sed "s/,/\\n/g" | sed "s/ //g" > tmp1 && gunzip -c ', MAP, ' | fgrep -f tmp1 - | cut -f 1 > tmp && unzip -p ', GENE_ZIP, ' ', gene_hor_zip , ' | fgrep -f tmp -', sep=' ') ), sep="\t", header=T, comment.char="#", row.names=1)
  all_genes_horizontal <- read.table(unz(description = GENE_ZIP, filename = gene_hor_zip), sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  cat ('Converting gene coverages (this will take a few minutes)\n')
  #print(summary(is.na(all_genes)))
  all_genes <- all_genes / colSums(all_genes,na.rm = T)
  all_genes <- all_genes[!(rownames(all_genes) == '-1' | rownames(all_genes) == 'mapped' | rownames(all_genes) == 'unassigned' | rownames(all_genes) == 'sum_annotated' | rownames(all_genes) == 'sum_not_annotated' | rownames(all_genes) == 'sum_not_annotated_and_annotated' | rownames(all_genes) == 'mapped_inserts_or_bases' | rownames(all_genes) == 'total_inserts_or_bases'),, drop=F]
  cat ('Converting horizontal gene coverages\n')
  all_genes_horizontal <- all_genes_horizontal[!(rownames(all_genes_horizontal) == '-1' | rownames(all_genes_horizontal) == 'mapped' | rownames(all_genes_horizontal) == 'unassigned' | rownames(all_genes_horizontal) == 'sum_annotated' | rownames(all_genes_horizontal) == 'sum_not_annotated' | rownames(all_genes_horizontal) == 'sum_not_annotated_and_annotated' | rownames(all_genes_horizontal) == 'mapped_inserts_or_bases' | rownames(all_genes_horizontal) == 'total_inserts_or_bases'),, drop=F]
  
  # KEGG data
  cat ('Loading (horizontal) ko, module and pathway coverages\n')
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, kegg_type, 'ko', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  tmp <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  ko <- tmp / colSums(tmp)
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, 'horizontal', 'ko', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  Hko <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, kegg_type, 'module', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  tmp <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  module <- tmp / colSums(tmp)
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, 'horizontal', 'module', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  Hmodule <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, kegg_type, 'pathway', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  tmp <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  pathway <- tmp / colSums(tmp)
  tmp <- read.table(unz(description = KEGG_ZIP, filename = paste(kegg_zip, 'horizontal', 'pathway', sep='.') ) , sep="\t", header=T, comment.char="#", row.names=1, check.names = F)
  Hpathway <- tmp[!(rownames(tmp) == '-1' | rownames(tmp) == 'mapped' | rownames(tmp) == 'unassigned' | rownames(tmp) == 'sum_annotated' | rownames(tmp) == 'sum_not_annotated' | rownames(tmp) == 'sum_not_annotated_and_annotated' | rownames(tmp) == 'mapped_inserts_or_bases' | rownames(tmp) == 'total_inserts_or_bases'),,drop=F]
  
  
  cat ('Processing gene map\n')
  # Because we have instances like when a gene is mapped to 2 KOs and then in turn into 2 modules, but these 2 modules only correspond to the second KO,
  # we cannot see the diff between these 2 using the gene abundance files, but actually need to sue the ko2module files...
  gene_map <- read.table(gzfile(description = MAP) , sep="\t", header=T, comment.char="", row.names=1, check.names = F)
  gene_map$ko <- as.vector(gene_map$ko)
  gene_map$module <- as.vector(gene_map$module)
  gene_map$pathway <- as.vector(gene_map$pathway)
  gene_map$ko[gene_map$ko == ''] <- 'NO KO'
  gene_map$module[gene_map$module == ''] <- 'NO MODULE'
  gene_map$pathway[gene_map$pathway == ''] <- 'NO PATHWAY'
  gene_map$ko[is.na(gene_map$ko)] <- 'NO KO'
  gene_map$module[is.na(gene_map$module)] <- 'NO MODULE'
  gene_map$pathway[is.na(gene_map$pathway)] <- 'NO PATHWAY'
  # load only genes with KOs
   tmp <- gene_map[gene_map$ko != '' & gene_map$ko != 'NO KO',]
  #   # module
  #   map <- unique(cbind(as.vector(tmp$ko), as.vector(tmp$module)))
  #   p1 <- strsplit(map[,1], ',')
  #   p2 <- strsplit(map[,2], ',')
  #   ko2moduleList <- list()
  #   for (list in c(1:length(p1))) {
  #     ko2moduleList[[list]] <- expand.grid(p1[[list]], p2[[list]] )
  #   }
  #   ko2module <- as.data.frame(rbindlist (ko2moduleList))
  
    
  # gene
  cat ('Processing gene2ko\n')
  map <- cbind(as.vector(rownames(tmp)), as.vector(tmp$ko))
  map2a <- as.data.frame(map[!grepl(',', map[,2]),])
  colnames(map2a) <- c('gene', 'ko')
  map2b <- map[grepl(',', map[,2]),]
  p1 <- strsplit(map2b[,1], ',')
  p2 <- strsplit(map2b[,2], ',')
  gene2koList <- list()
  for (list in c(1:length(p1))) {
    gene2koList[[list]] <- expand.grid(p1[[list]], p2[[list]] )
  }
  gene2ko <- as.data.frame(rbindlist (gene2koList))
  colnames(gene2ko) <- c('gene', 'ko')
  gene2ko <- rbind(gene2ko, map2a)
  
  if (TAXONOMIC !="0") {
    cat ('Processing gene2taxa\n')
    TAXA <- read.table(TAXONOMIC, header = F, sep = "\t", fill = T, row.names=1)
    TAXA[TAXA == ''] <- NA
    TAXA[TAXA == 'NA'] <- NA
  }
  
  cat ("Removing objects\n")
  rm(list=(c('gene2koList', 'p1', 'p2', 'map', 'map2a', 'map2b')))
  all_genes_rownames <- rownames(all_genes)
  cat('Saving RData file for quicker loading next time (this will take a few minutes)\n')
  ONLY_THESE_PATHWAYS <- NULL
  save.image(file = paste(SAVE, 'all.RData', sep='.'))
} # end run if RData file doesn't exist

SAVE <- args[7]
PATHWAYS <- strsplit(args[11], " ")[[1]]
if (PATHWAYS[1] != "0") {
  cat('Extracting selected pathways\n')
  # all_genes, all_genes_horiztonal, gene2ko, TAXA, all_genes_rownames, gene_map
  keep <- unique(as.vector(gene2ko[gene2ko[,2]  %in%  ko2pathway[ko2pathway[,2] %in% PATHWAYS,1],1]))
  all_genes <- all_genes[rownames(all_genes) %in% keep,,drop=F]
  all_genes_horizontal <- all_genes_horizontal[rownames(all_genes_horizontal) %in% keep,,drop=F]
  gene2ko <- gene2ko[gene2ko$gene %in% keep,]
  TAXA <- TAXA[rownames(TAXA) %in% keep,]
  all_genes_rownames <- keep
  gene_map <- gene_map[rownames(gene_map) %in% keep,]
  cat('Saving selected RData for MOCAT Viewer\n')
  ONLY_THESE_PATHWAYS <- PATHWAYS
  save.image(file = paste(SAVE, paste(sort(PATHWAYS), collapse='.'), 'RData', sep='.') )
}

cat('Saving minimal RData for MOCAT Viewer\n')
save(list=c('kegg2name', 'gene2ko', 'ko2pathway', 'TAXA', 'all_genes_rownames', 'ONLY_THESE_PATHWAYS', 'pathway'), file = paste(SAVE, 'minimal.RData', sep='.'))
cat('Initial startup completed\n')
