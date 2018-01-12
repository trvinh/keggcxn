library(rcytoscapejs)
library(DT)

shinyServer(function(input, output, session) {
  session$allowReconnect(TRUE)
  
  ##### render upload main input #####
  output$mainInput <- renderUI({
    fileInput("mainInput","Upload input file:")
  })
  
  ##### render list of pathway names and IDs #####
  pathDf <- reactive({
    pathDf <- as.data.frame(read.table("data/pathways.list", quote = "", sep='\t',header=FALSE,check.names=FALSE,comment.char="",colClasses=c(rep("factor",4))))
    colnames(pathDf) <- c("pathID","pathName","pathGroup","pathType")
    return(pathDf)
  })
  
  output$pathSelect <- renderUI({
    pathDf <- pathDf()
    pathDf$pathName <- factor(pathDf$pathName, levels = unique(pathDf$pathName))
    
    selectInput('pathName', label = "Pathway name:",
                choices = as.list(levels(pathDf$pathName)),
                selected = levels(pathDf$pathName)[11])
  })
  
  output$pathID.ui <- renderUI({
    pathDf <- pathDf()
    pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
    
    selectInput('pathID', label = "Pathway ID:",
                choices = as.list(levels(pathDf$pathID)),
                selected = levels(pathDf$pathID)[11])
  })
  
  ### Update value of selected pathID and pathName
  observe({
    if(!is.null(input$pathID)){
      pathDf <- pathDf()
      pathDf$pathName <- factor(pathDf$pathName, levels = unique(pathDf$pathName))
      selectedName <- pathDf[pathDf$pathID == input$pathID,]$pathName
      
      updateSelectInput(
        session, 'pathName', label = "Pathway Name:",
        choices = as.list(levels(pathDf$pathName)),
        selected = selectedName
      )
    }
  })
  
  observe({
    if(!is.null(input$pathName)){
      pathDf <- pathDf()
      pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
      selectedID <- pathDf[pathDf$pathName == input$pathName,]$pathID
      
      updateSelectInput(
        session, 'pathID', label = "Pathway ID:",
        choices = as.list(levels(pathDf$pathID)),
        selected = selectedID
      )
    }
  })
  
  ### show pathway description
  output$pathDesc <- renderText({
    pathDf <- pathDf()
    if(is.null(input$pathID)){paste("")} 
    else {paste(pathDf[pathDf$pathID == input$pathID,]$pathGroup,pathDf[pathDf$pathID == input$pathID,]$pathType, sep = ", ")}
    # pathDesc <- pathDf[pathDf$pathID == input$pathID,]$pathGroup
  })
  
  ##### load cxn file for selected pathway #####
  checkCxnFile <- reactive({
    cxnFile <- suppressWarnings(paste0("data/keggcxn/",input$pathID,".cxn"))
    return(file.exists(cxnFile))
  })
  
  output$msg <- renderUI({
    link <- paste0("http://www.genome.jp/kegg-bin/show_pathway?ko",input$pathID)
    if(checkCxnFile() == FALSE){
      msg <- paste0(
        "<p><em>No reaction has been found in this pathway or it does not have KGML file!! </em></p>
        <p><em>Please check:</em></p>
        <a href=\"",link,"\" target=\"_blank\">",link,"</a>"
      )
      HTML(msg)
    }
    else {
      msg <- paste0(
        "<p><em>KEGG pathway link: </em><a href=\"",link,"\" target=\"_blank\">",link,"</a></p>"
      )
      HTML(msg)
    }
  })
  
  networkDf <- reactive({
    if(is.null(input$pathID)){return ()}
    else {
      cxnFile <- suppressWarnings(paste0("data/keggcxn/",input$pathID,".cxn"))
      if(file.exists(cxnFile)){
        inputDf <- as.data.frame(read.table(cxnFile, sep='\t',header=F,check.names=FALSE,comment.char=""))
        colnames(inputDf) <- c("source","target","type")
        return(inputDf)
      } else {
        return()
      }
    }
  })
  
  ##### read input file #####
  annoDf <- reactive({
    filein <- input$mainInput
    
    if(is.null(filein)){ return()}
    annoDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE))
    colnames(annoDf)[2] <- "name"
    return (annoDf)
  })
  
  ##### get list of nodes and mark annotated nodes from input file #####
  nodeData <- reactive({
    inputDf <- networkDf()
    
    if(is.null(inputDf)){return()}
    
    ### join source and target KOs to get full list of nodes
    sourceDf <- data.frame(inputDf$source,stringsAsFactors = FALSE)
    colnames(sourceDf) <- "name"
    targetDf <- data.frame(inputDf$target,stringsAsFactors = FALSE)
    colnames(targetDf) <- "name"
    
    nodeDf <- unique(rbind(sourceDf,targetDf))
    nodeDf$id <- nodeDf$name
    
    ### map annotated nodes
    annoDf <- annoDf()

    if(!is.null(annoDf)){
      joinedDf <- merge(nodeDf,annoDf, by="name", all.x=TRUE)
      
      ### remove un-annotated nodes if requested
      if(input$filterNodes == "yes"){
        joinedDf <- joinedDf[complete.cases(joinedDf),]
      }
      
      ### change colors of nodes based on patristic distance of reference annotations
      subNodeDf <- joinedDf[,c("id","name","color")]
      if(length(unique(complete.cases(joinedDf))) > 1){
        subNodeDf[is.na(subNodeDf$color),]$color <- "#878787"
      }
      
      ### return nodeDf
      return(subNodeDf)
    } else {
      return(nodeDf)
    }
    
  })
  
  ##### get list of edges 
  edgeData <- reactive({
    inputDf <- networkDf()
    edgeDf <- inputDf[,(1:2)]
    
    return(edgeDf)
  })
  
  
  
  # setwd("/Users/trvinh/work/OLD/koAnnotation/lcaAnnotation/kegg_connectivity")
  # inputDf <- as.data.frame(read.table("00030.cxn", sep='\t',header=F,check.names=FALSE,comment.char=""))
  # colnames(inputDf) <- c("source","target","type")
  # # inputDf <- networkDf()
  # 
  # sourceDf <- data.frame(inputDf$source,stringsAsFactors = FALSE)
  # colnames(sourceDf) <- "name"
  # targetDf <- data.frame(inputDf$target,stringsAsFactors = FALSE)
  # colnames(targetDf) <- "name"
  # 
  # nodeDf <- unique(rbind(sourceDf,targetDf))
  # nodeDf$id <- nodeDf$name
  # 
  # edgeDf <- inputDf[,(1:2)]
  # 
  # annoDf <- as.data.frame(read.table("/Users/trvinh/work/R_projects/keggcxn/data/lca_micros.ko", sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE))
  # 
  # colnames(annoDf)[2] <- "name"
  # 
  # joinedDf <- merge(nodeDf,annoDf, by="name", all.x=TRUE)
  # head(joinedDf)
  # 
  # subNodeDf <- joinedDf[,c("id","name","color")]
  # subNodeDf[is.na(subNodeDf$color),]$color <- "#878787"
  # 
  # nodeData <- subNodeDf
  # edgeData <- edgeDf
  
  
  # network <- read.table("cbioportal_top1.sif", sep="\t", stringsAsFactors=FALSE)
  # 
  # colnames(network) <- c("source", "interaction", "target")
  # 
  # #maxInteractions <-  input$maxInteractions
  # maxInteractions <- 50
  # 
  # if(nrow(network) <= maxInteractions) {
  #   maxInteractions <- nrow(network)
  # } else {
  #   maxInteractions <- maxInteractions
  # }
  # 
  # network <- network[1:maxInteractions, ]
  # 
  # edgeList <- network[, c("source","target")]
  # 
  # nodes <- unique(c(edgeList$source, edgeList$target))
  # 
  # id <- nodes
  # name <- nodes
  # addLinks <- TRUE
  # 
  # if(addLinks) {
  #   href <- paste0("https://www.google.com/search?q=", nodes)
  #   tooltip <- paste0("https://www.google.com/search?q=", nodes)
  #   nodeData <- data.frame(id, name, href, tooltip, stringsAsFactors=FALSE)
  # } else {
  #   nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  # }
  # 
  # nodeData$color <- rep("#888888", nrow(nodeData))
  # nodeData$color[which(grepl("[a-z]", nodes))] <- "#FF0000"
  # 
  # nodeData$shape <- rep("ellipse", nrow(nodeData))
  # nodeData$shape[which(grepl("[a-z]", nodes))] <- "octagon"
  # 
  # edgeData <- edgeList

  # NOTE: Reactive variables used as functions networkReactive()
  networkReactive <- reactive({
    if(is.null(input$connectedNodes)) {
      return(network)
    } else {
      t1 <- which(network$source %in% input$connectedNodes)
      t2 <- which(network$target %in% input$connectedNodes)
      idx <- unique(c(t1, t2))
      return(network[idx,])
    }
  })
  
  # output$nodeDataTable <- DT::renderDataTable({
  #   tmp <- nodeData[which(id == input$clickedNode),]
  #   DT::datatable(tmp, filter='bottom', style='bootstrap', options=list(pageLength=5))
  # })
  # 
  # output$edgeDataTable <- DT::renderDataTable({
  #   DT::datatable(networkReactive(), filter='bottom', style='bootstrap', options=list(pageLength=5))
  # })
  
  output$clickedNode = renderPrint({
    input$clickedNode
  })
  
  output$connectedNodes = renderPrint({
    input$connectedNodes
  })
  
  output$plot <- renderRcytoscapejs({
    inputDf <- networkDf()
    if(is.null(inputDf)){return()}
    
    cyNetwork <- createCytoscapeJsNetwork(nodeData(), edgeData())
    rcytoscapejs(nodeEntries=cyNetwork$nodes, edgeEntries=cyNetwork$edges)
  })
  
  observeEvent(input$saveImage, {
    # NOTE: Message cannot be an empty string "", nothing will happen    
    session$sendCustomMessage(type="saveImage", message="NULL")
  })
})
