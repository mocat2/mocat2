library(shiny)
require(ggplot2)

load('minimal.RData')

TAXA <- as.matrix(TAXA)
TAXA[is.na(TAXA)] <- 'Unknown'
all_phylum <- sort(unique(TAXA[,3]))
all_genus <- sort(unique(TAXA[,7]))
all_species <- sort(unique(TAXA[,8]))
if (length(ONLY_THESE_PATHWAYS)>0) {
  pathway_options1 <- ONLY_THESE_PATHWAYS
} else {
  pathway_options1 <- sort(unique(as.vector(ko2pathway[ko2pathway[,1] %in% as.vector(gene2ko[gene2ko[,1] %in% all_genes_rownames,2]),2])))  
}
pathway_options2 <- as.vector(paste(
  as.vector(kegg2name[kegg2name[,1] %in% pathway_options1,1]), ' (',
  as.vector(kegg2name[kegg2name[,1] %in% pathway_options1,2]), ')',sep=''))
pathway_options3 <- as.vector(paste(
  as.vector(pathway_options1[!pathway_options1 %in% kegg2name[,1]]), '(no desc)'))
pathway_options <- NULL
if (pathway_options2[1] != ' ()') {
  pathway_options <- c(pathway_options, pathway_options2)
}
if (pathway_options3[1] != ' (no desc)') {
  pathway_options <- c(pathway_options, pathway_options3)
}

types <- c('Gene coverages',
           'Horizontal gene coverages',
           'KO coverages',
           'Horizontal KO coverages',
           'Module coverages',
           'Horizontal module coverages',
           'Pathway coverages',
           'Horizontal pathway coverages',
           'Gene species abundances',
           'Gene genus abundances',
           'Gene phylum abundances',
           'KO species abundances',
           'KO genus abundances',
           'KO phylum abundances',
           'Module species abundances',
           'Module genus abundances',
           'Module phylum abundances',
           'Pathway species abundances',
           'Pathway genus abundances',
           'Pathway phylum abundances'
           
           )


shinyUI(
  fluidPage(
    
    
    titlePanel = 'MOCAT Viewer',
    sidebarLayout(
      sidebarPanel(
        fixedRow(align="center",checkboxInput('plotPreview', 'Update preview automatically', value = FALSE)),
        fixedRow(align="center",downloadButton('downloadPlot', 'Save plot')),
        hr(),
        numericInput("width", "Image width", 2000),
        numericInput("height", "Image height", 1000),
        hr(),
        
        uiOutput("samples"),
        selectInput(
          'pathway', 'Selected pathways', choices = pathway_options, selected = '', multiple = T
        ),
        selectInput(
          'module', 'Selected Modules',
          choices = '',
          selected = '', multiple = T
        ),        
        selectInput(
          'KO', 'Selected KOs', choices = '', selected = '', multiple = T
        ),
        selectInput(
          'type', 'Plot the following', choices = types, selected = types, multiple = T
        ),
        checkboxInput('printSpecies', 'Print species names', value = FALSE),
        hr(),
        selectInput('excludePhylum', 'Exclude phyla', choices = all_phylum, multiple = T),
        selectInput('excludeGenus', 'Exclude genera', choices = all_genus, multiple = T),
        selectInput('excludeSpecies', 'Exclude species', choices = all_species, multiple = T),
        hr(),
        h4("Plotting settings"),
        numericInput('GENE_WIDTH', 'Gene width', 20, min = NA, max = NA, step = NA),
        numericInput('GENE_HEIGHT', 'Gene coverage height', 10, min = NA, max = NA, step = NA),
        numericInput('MAX_GENES', 'Genes per row', 50, min = NA, max = NA, step = NA),
        numericInput('KO_HEIGHT', 'KO coverage height', 30, min = NA, max = NA, step = NA),
        numericInput('MOD_HEIGHT', 'Module coverage height', 50, min = NA, max = NA, step = NA),
        numericInput('KO_ABOVE_GENE', 'Space b/w genes and KOs', 40, min = NA, max = NA, step = NA),
        numericInput('MOD_ABOVE_GENE', 'Space b/w genes and modules', 80, min = NA, max = NA, step = NA),
        numericInput('PATHWAY_ABOVE_GENE', 'Space b/w genes and pathways', 120, min = NA, max = NA, step = NA),
        numericInput('SPACE_BW_KO_COV_HOR', 'Space between KO and module cov. and hor. cov.', 0.2, min = NA, max = NA, step = NA),
        numericInput('SPACE_BW_KOS', 'Space b/w KOs', 10, min = NA, max = NA, step = NA),
        numericInput('SPACE_BW_MODULES', 'Space b/w modules', 1000, min = NA, max = NA, step = NA),
        numericInput('SPACE_BW_GROUPS', 'Space b/w groups', 3000, min = NA, max = NA, step = NA),
        numericInput('SPACE_FOR_EACH_TAXA', 'height gene taxa', 10, min = NA, max = NA, step = NA),
        numericInput('SPACE_BW_GENE_ROWS', 'Space b/w gene rows', 40, min = NA, max = NA, step = NA)
        #submitButton("Update View")
      ),
      mainPanel(
        plotOutput(outputId = "plot", height = "1000px", width = "2000px")
        #plotOutput(outputId = "plot", height = paste(output()$height, 'px', sep=""), width = paste(output()$width, 'px', sep=''))
      )
    )
  ))