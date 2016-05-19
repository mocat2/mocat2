library(shiny)

# Load objects saved previously
base <- as.matrix(read.table("base", header=F))[,1]
load(paste(base, '.results.RData', sep=''))

all.taxa <- sort(unique(colnames(reads[['sum_samples']][[1]]), colnames(inserts[['sum_samples']][[1]])))

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("MOCAT : Multiple Mappers"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    selectInput("dataset", "Choose a type:", choices = c("reads", "inserts")),
    numericInput("obs", "Number of top multiple mappings and samples to view:", 10),
    selectInput("sample", "Select sample:", choices = samples),
    selectInput("taxa1", "Select Taxa 1:", choices = all.taxa),
    selectInput("taxa2", "Select Taxa 2:", choices = all.taxa),
    helpText("Note that TAXA1 vs TAXA1 are not technically counted as multiple mappers, but the diagonal has the total number of reads mapped to that TAXA")
    
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  mainPanel(
    h2("Summary All Samples"),
    h4("Taxa vs Taxa"),
    tableOutput("t1t2a"),
    h4("Top Taxa vs Taxa"),
    tableOutput("top"),
    
    h4(textOutput("caption1")),
    tableOutput("ALLtaxa1VSall"),
    
    h2("Summary Top Samples"),
    h4("Taxa vs Taxa"),
    tableOutput("t1t2as"),
    
    h2(textOutput("caption")),
    h4("Taxa vs Taxa"),
    tableOutput("t1t2s"),
    
    h4(textOutput("caption2")),
    tableOutput("SAMPLEtaxa1VSall"),
        
    h4("Top Taxa vs Taxa"),
    tableOutput("view"),
    h4("Number of reads mapped"),
    tableOutput("view2")
    
    
  )
))
