library(shiny)

# Load objects saved previously
base <- as.matrix(read.table("base", header=F))[,1]
load(paste(base, '.results.RData', sep=''))

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive(function() {
    switch(input$dataset,
           "reads" = reads,
           "inserts" = inserts,
           )
  })
  
  taxa1 <- reactive(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    switch(input$taxa1, colnames(sum.list.read[[1]]))
  })
  
  taxa2 <- reactive(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    switch(input$taxa2, colnames(sum.list.read[[1]]))
  })
  
  
  output$caption <- reactiveText(function() {
    paste("Summary", input$sample)
  })
  output$caption1 <- reactiveText(function() {
    paste(input$taxa1, 'vs Others')
  })
  output$caption2 <- reactiveText(function() {
    paste(input$taxa1, 'vs Others')
  })
  
  output$t1t2a <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    
    t1 <- sum.read[input$taxa1, input$taxa1]
    t2 <- sum.read[input$taxa2, input$taxa2]
    
    sub <- sum.read[input$taxa1, input$taxa2, drop=F]
    frac <- sum.read[input$taxa1, input$taxa2] / sum.read.interactions * 100
    names(frac) <- '% of total interactions'
    ret <- cbind(sub, frac, sub/t1, sub/t2)
    colnames(ret) <- c(colnames(sub), names(frac), paste('% multiple mappers (', input$taxa1, ')', sep=''), paste('% multiple mappers (', input$taxa2, ')', sep=''))
    return(ret)
  })

  output$t1t2s <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    ret <- sum.list.read[[input$sample]][input$taxa1, input$taxa2, drop=F]
    n <- colnames(ret)
    ret <- cbind(ret, ret[1,1] / read.stats[[input$sample]][3,3] * 100)
    colnames(ret) <- c(n, '% of sample interactions')
    return(ret)
  })

  output$t1t2as <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    all <- matrix(nrow=length(names(sum.list.read)), ncol=2)
    for (i in 1:length(names(sum.list.read))) {
      all[i,1] <- sum.list.read[[i]][input$taxa1, input$taxa2]
      all[i,2] <- all[i,1] / read.stats[[i]][3,3] * 100
    }
    rownames(all) <- names(sum.list.read)
    colnames(all) <- c(paste(input$taxa1, input$taxa2, sep=':'), '% of sample interactions')
    all.order <- order(all[,1], decreasing=T)
    all <- all[all.order,,drop=F]
    if (input$obs <= dim(all)[1]) {
      max <- input$obs
    } else {
      max <- dim(all)[1]
    }
    all <- all[1:max,,drop=F]
    return (all)
  })
  
  output$top <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    ret <- as.matrix(head(top.read, n = input$obs))
    ret <- cbind(ret, ret[,1]/sum.read.interactions*100)
    colnames(ret) <- c('Multiple mappers', '% of total interactions')
    return (ret)
  })
  
  # Show the first "n" observations
  output$view <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    ret <- as.matrix(head(top.list.read[[input$sample]], n = input$obs))
    ret <- cbind(ret, ret[,1] / read.stats[[input$sample]][3,3] * 100)
    colnames(ret) <- c(paste('Multiple mappers', sep=' : '), '% of sample interactions')
    return (ret)
  })
  
  output$view2 <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    if (input$obs <= dim(read.stats[[input$sample]])[1]) {
      max <- input$obs + 5
    } else {
      max <- dim(read.stats[[input$sample]])[1]
    }
    ret <- read.stats[[1]][1:(max),3,drop=F]
    rownames(ret) <- read.stats[[1]][1:(max),2]
    colnames(ret) <- 'reads mapped to X taxa'
    return(ret)
  })
  
  output$ALLtaxa1VSall <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    all <- top.read[grep(input$taxa1, names(top.read), perl=T)]
    all <- all[order(all, decreasing=T)]
    all <- as.matrix(head(all, n = input$obs))
    all <- cbind(all, all[,1]/sum.read.interactions*100)
    colnames(all) <- c('Multiple mappers', '% of total interactions')
    return(all)
  })
    
  output$SAMPLEtaxa1VSall <- reactiveTable(function() {
    tmp <- datasetInput()
    sum.read <- tmp[['sum']]
    sum.read.interactions <- tmp[['interactions']]
    sum.list.read <- tmp[['sum_samples']]
    read.stats <- tmp[['stats']]
    top.read <- tmp[['top']]
    top.list.read <- tmp[['top_samples']]
    
    all <- top.list.read[[input$sample]][grep(input$taxa1, names(top.list.read[[input$sample]]), perl=T)]
    all <- all[order(all, decreasing=T)]
    all <- as.matrix(head(all, n = input$obs))
    all <- cbind(all, all[,1] / read.stats[[input$sample]][3,3] * 100 )
    colnames(all) <- c('Multiple mappers', '% of sample interactions')
    return(all)
  })
    
})