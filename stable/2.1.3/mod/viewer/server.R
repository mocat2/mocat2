library(shiny)
require(ggplot2)

GROUPING <- '0'

# METADATA, may change, so it's outside the RData file
cat ('Processing metadata\n')
SAMPLE.GROUPS <- list()
if (GROUPING !="0") {
  metadata <- as.matrix(read.table(METADATA, quote="", header=T, row.names=1, sep="\t", check.names=F, na.strings="NA", comment.char=""))  
  metadata <- as.data.frame(metadata[rownames(metadata) %in% colnames(ko),,drop=F])
  levels <- levels(as.factor(metadata[,GROUPING]))  
  for (i in 1:length(levels)) {
    SAMPLE.GROUPS[[levels[i]]] <- rownames(metadata)[metadata[,GROUPING] == levels[i]][!is.na(rownames(metadata)[metadata[,GROUPING] == levels[i]])]
  }
} else {
  SAMPLE.GROUPS$Samples <- colnames(ko)
}
cat ('Ready and waiting\n')

shinyServer(function(input, output, session) {
  
  
  output$samples <- renderUI({
    lapply(1:length(names(SAMPLE.GROUPS)), function(i) {
      selectInput(
        inputId = paste0("gr", i), names(SAMPLE.GROUPS)[i], choices = as.vector(SAMPLE.GROUPS[[i]]), multiple = T
      )
    })
  })
  
  
  #################################################
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('Figure.png', sep='')
    },
    content = function(file) {
      if(length(input$pathway)==0) return (NULL)
      if(length(input$module)==0) return (NULL)
      if(length(input$KO)==0) return (NULL)
      source('plot3.R', local = TRUE)
      png('Figure.png', width = as.numeric(as.vector(input$width)), height = as.numeric(as.vector(input$height)))
      p <- ggplot() +
        geom_tile(data=TO_PLOT, aes(X, Y, fill=Color, alpha=alpha, height=height, width=width), colour=NA) +
        scale_fill_manual(values = COLORS, name="", na.value="pink") +
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
        theme(legend.text = element_text(size = rel(1.5), colour = "black")) +
        guides(fill=guide_legend(ncol=2))
      print(p)
      dev.off()
    },contentType = 'image/png')  
  ####### PLOT #######
  #################################################
  
  
  #################################################
  output$plot <- renderPlot({
    if(!input$plotPreview) return (NULL)
    if(length(input$pathway)==0) return (NULL)
    if(length(input$module)==0) return (NULL)
    if(length(input$KO)==0) return (NULL)
    samples <- NULL
    for (i in c(1:length(names(SAMPLE.GROUPS)))){
      if( length(as.vector(input[[paste0("gr", i)]])  )) {
        samples <- c(samples, as.vector(input[[paste0("gr", i)]]))
      }
    }
    if(length(samples) == 0) { return (NULL) }
    
    source('plot3.R', local = TRUE)
    
    p <- ggplot() +
      geom_tile(data=TO_PLOT, aes(X, Y, fill=Color, alpha=alpha, height=height, width=width), colour=NA) +
      scale_fill_manual(values = COLORS, name="", na.value="pink") +
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
      theme(legend.text = element_text(size = rel(1.5), colour = "black")) +
      guides(fill=guide_legend(ncol=2))
    print(p)
    rm(TO_PLOT)
  })
  ####### PLOT #######
  #################################################
  
  
  outMOD = reactive({
    if(length(input$pathway)==0) return (NULL)
    # biggest pathway here is a reminant from when we selected only one pathway, now we can have many
    X <- strsplit(as.vector(input$pathway), " \\(")
    X <- sapply(X, "[[", 1)
    all_kos_in_biggest_pathway <- unique(as.vector(ko2pathway[ko2pathway[,2] %in% X,1]))
    all_modules_in_biggest_pathway <- unique(as.vector(ko2module[ko2module[,1] %in% all_kos_in_biggest_pathway,2]))
    order_module <- rev(sort(table(ko2module[ko2module[,2] %in% all_modules_in_biggest_pathway,2])))
    order_module1 <- names(order_module[order_module > 0])
    order_module2 <- as.vector(paste(
      as.vector(kegg2name[kegg2name[,1] %in% order_module1,1]), ' (',
      as.vector(kegg2name[kegg2name[,1] %in% order_module1,2]), ')',sep=''))
    order_module3 <- as.vector(paste(
      as.vector(order_module1[!order_module1 %in% kegg2name[,1]]), '(no desc)'))
    order_module <- NULL
    if (order_module2[1] != ' ()') {
    order_module <- c(order_module, order_module2)
    }
    if (order_module3[1] != ' (no desc)') {
      order_module <- c(order_module, order_module3)
    }
    return(order_module)
  })
  observe({
    updateSelectInput(session, "module", choices = outMOD(), selected = outMOD())
  })
  outKO = reactive({
    if(length(input$pathway)==0) return (NULL)
    if(length(input$module)==0) return (NULL)
    X <- strsplit(as.vector(input$module), " \\(")
    X <- sapply(X, "[[", 1)
    pw <- strsplit(as.vector(input$pathway), " \\(")
    pw <- sapply(pw, "[[", 1)
    order_ko <- rev(sort(table(gene2ko[gene2ko[,2] %in% unique(gene2ko[,2]),2]))) # here we pick only form KOs in this sample
    order_ko <- order_ko[order_ko > 0]
    order_ko <- order_ko[names(order_ko) %in% as.vector(ko2pathway[ko2pathway[,2] %in% pw,1])] # these are all the KOs in this pathway
    
    current_kos <- as.vector(ko2module[ko2module[,2] %in% X,1]) # get the KOs for these modules, including NO MODULE ones
    current_kos1 <- names(rev(sort(order_ko[names(order_ko) %in% current_kos])))
    #current_kos1 <- names(current_kos[current_kos > 0])
    current_kos2 <- as.vector(paste(
      as.vector(kegg2name[kegg2name[,1] %in% current_kos1,1]), ' (',
      as.vector(kegg2name[kegg2name[,1] %in% current_kos1,2]), ')',sep=''))
    current_kos3 <- as.vector(paste(
      as.vector(current_kos1[!current_kos1 %in% kegg2name[,1]]), '(no desc)'))
    current_kos <- NULL
    if (current_kos2[1] != ' ()') {
      current_kos <- c(current_kos, current_kos2)
    }
    if (current_kos3[1] != ' (no desc)') {
      current_kos <- c(current_kos, current_kos3)
    }
    return(current_kos)
  })
  observe({
    updateSelectInput(session, "KO", choices = outKO(), selected = outKO())
  })
})
