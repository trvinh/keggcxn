# if (!require("network")) {install.packages("network")}
# if (!require("igraph")) {install.packages("igraph")}
if (!require("DT")) {install.packages("DT")}
if (!require("dplyr")) {install.packages("dplyr")}

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
                selected = levels(pathDf$pathName)[93])
  })
  
  output$pathID.ui <- renderUI({
    pathDf <- pathDf()
    pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
    
    selectInput('pathID', label = "Pathway ID:",
                choices = as.list(levels(pathDf$pathID)),
                selected = levels(pathDf$pathID)[93])
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
        colnames(inputDf) <- c("from","to","interaction")
        
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
    colnames(annoDf)[2] <- "id"
    return (annoDf)
  })
  
  ##### get list of nodes and mark annotated nodes from input file #####
  nodeData <- reactive({
    inputDf <- networkDf()
    
    if(is.null(inputDf)){return()}
    
    ### join source and target KOs to get full list of nodes for reference network
    sourceDf <- data.frame(inputDf$from,stringsAsFactors = FALSE)
    colnames(sourceDf) <- "id"
    targetDf <- data.frame(inputDf$to,stringsAsFactors = FALSE)
    colnames(targetDf) <- "id"
    
    nodeDf <- unique(rbind(sourceDf,targetDf))
    nodeDf$group <- "reference"
    
    ### map annotated nodes to reference net
    annoDf <- annoDf()

    if(!is.null(annoDf)){
      joinedDf <- merge(nodeDf,annoDf, by="id", all.x=TRUE)
      
      ### remove un-annotated nodes if requested
      if(input$filterNodes == "yes"){
        joinedDf <- joinedDf[complete.cases(joinedDf),]
        joinedDf$group <- "annotated"
      }
      
      ### change group type for annotated nodes
      if(length(unique(complete.cases(joinedDf))) > 1){
        joinedDf[!is.na(joinedDf$geneID),]$group <- "annotated"
      }
      
      ### process duplicated node IDs
      dupID <- joinedDf[duplicated(joinedDf$id),]$id
      if(length(dupID) > 0){
        # join gene IDs
        joinedDf$joinedGeneID <- joinedDf$geneID
        for(i in 1:length(dupID)){
          joinedDf[joinedDf$id == dupID[i],]$joinedGeneID <- toString(joinedDf[joinedDf$id == dupID[i],]$geneID)
        }
        
        # get max FAS for this joined node
        aggrFAS <- aggregate(joinedDf[,"FAS"],list(joinedDf$id),FUN="max")
        colnames(aggrFAS) <- c("id","FAS")
        
        # remove duplicated rows
        joinedDf <- joinedDf[!duplicated(joinedDf$id),]
        joinedDf <- merge(aggrFAS, joinedDf, by="id", all.x = TRUE)
        
        # return joinedDf after removing duplicated nodes
        joinedDf <- joinedDf[,c("id","group","FAS.x","joinedGeneID","color","patristicDist")]
        colnames(joinedDf) <- c("id","group","FAS","geneID","color","patristicDist")
      }
      
      ### node color, title, label, (size)
      if(nrow(joinedDf[is.na(joinedDf$color),]) > 0){
        joinedDf[is.na(joinedDf$color),]$color <- "#878787"
      }
      
      joinedDf <- rename(joinedDf, title = geneID)
      if(nrow(joinedDf[is.na(joinedDf$title),]) > 0){
        joinedDf[is.na(joinedDf$title),]$title <- as.character(joinedDf[is.na(joinedDf$title),]$id)
      }
      
      joinedDf$label <- joinedDf$id
      
      joinedDf$value <- as.integer(100^(joinedDf$FAS)*10)
      if(nrow(joinedDf[is.na(joinedDf$FAS),]) > 0){
        joinedDf[is.na(joinedDf$FAS),]$value <- 1
      }
      
      ### return node data
      subNodeDf <- joinedDf[,c("id","label","title","group","value","color")]
      return(subNodeDf)
    } else {
      nodeDf$color <- "#878787"
      nodeDf$title <- nodeDf$id
      nodeDf$label <- nodeDf$id
      nodeDf$shape <- "box"
      return(nodeDf)
    }
    
  })
  
  ##### get list of edges #####
  edgeData <- reactive({
    inputDf <- networkDf()
    edgeDf <- inputDf[,(1:2)]
    return(edgeDf)
  })
  
  # 
  # networkReactive <- reactive({
  #   if(is.null(input$pathID)){return ()}
  #   network <- networkDf()
  #   if(is.null(input$connectedNodes)) {
  #     return(network)
  #   } else {
  #     t1 <- which(network$source %in% input$connectedNodes)
  #     t2 <- which(network$target %in% input$connectedNodes)
  #     idx <- unique(c(t1, t2))
  #     return(network[idx,])
  #   }
  # })
  
  # 
  # ##### print list of connected nodes for clicked nodes #####
  
  
  ##### RENDER NETWORK #####
  networkPlot <- reactive({
    nodeData <- nodeData()
    edgeData <- edgeData()
    if(is.null(nodeData)){return()}
    
    if(nlevels(as.factor(nodeData$group)) > 1){
      network <- visNetwork(nodeData, edgeData) %>%
        # group 2 types of nodes and show legend
        visGroups(groupname = "reference", shape = "box", color = "#878787", shadow = FALSE) %>%
        visGroups(groupname = "annotated", shape = "circle", color = "#75F514", shadow = TRUE) %>%
        visLegend(width = 0.1, position = "right") %>%
        # highlight nearest nodes
        visOptions(highlightNearest = list(enabled = T, hover = T),  nodesIdSelection = TRUE, selectedBy = "group", collapse = TRUE) %>%
        # set fix randomSeed for getting the same network everytime
        visLayout(randomSeed = 123) %>%
        # physical parameters
        visPhysics(stabilization = TRUE, maxVelocity = input$maxSpeed) %>%
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}") %>% 
        # export
        visExport(type = "pdf", name = input$pathID, float = "top", style = "warning") %>% addExport(pdf = TRUE)
    } else {
      network <- visNetwork(nodeData, edgeData) %>%
        # highlight nearest nodes
        visOptions(highlightNearest = list(enabled = T, hover = T),  nodesIdSelection = TRUE, collapse = TRUE) %>%
        # set fix randomSeed for getting the same network everytime
        visLayout(randomSeed = 123) %>%
        # physical parameters
        visPhysics(stabilization = TRUE, maxVelocity = input$maxSpeed) %>%
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}") %>% 
        # export
        visExport(type = "pdf", name = input$pathID, float = "top", style = "warning") %>% addExport(pdf = TRUE)
    }
    
    return(network)
  })
  
  output$plot <- renderVisNetwork({
    networkPlot()
  })
  
  ##### render selected node #####
  output$clickedNode = renderText({
    if (is.null(input$click)) {return()}
    
    nodeData <- nodeData()
    paste0(input$click," (desc: ",nodeData[nodeData$id == input$click,]$title,")")
  })
  
  ##### print list of connected nodes for clicked nodes #####
  observe({
    if(length(input$click) == 1){
      visNetworkProxy("plot") %>%
        visGetConnectedNodes(id = input$click, input = "connected_nodes")
    }
  })
  
  connectedNodes <- reactive({
    if (is.null(input$connected_nodes)) {return()}
    
    nodeData <- nodeData()
    
    # get node tooltip
    connectedNodes <- data.frame(input$connected_nodes)
    colnames(connectedNodes)[1] <- "id"
    connectedNodesOut <- merge(connectedNodes, nodeData, by = "id", all.x = TRUE)
    connectedNodesOut <- connectedNodesOut[!duplicated(connectedNodesOut),]
    # return list
    connectedNodesOut <- connectedNodesOut[,c("id","title")]
    colnames(connectedNodesOut) <- c("Node","Desc")

    return(connectedNodesOut)
  })

  output$connectedNodes <- DT::renderDataTable({
    DT::datatable(connectedNodes(), style='bootstrap', options=list(pageLength=5, sDom  = '<"top">lrt<"bottom">ip'))
  })
  
})
