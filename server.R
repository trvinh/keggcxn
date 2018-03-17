if (!require("DT")) {install.packages("DT")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("reshape")) {install.packages("reshape")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("visNetwork")) {install.packages("visNetwork")}
if (!require("igraph")) {install.packages("igraph")}
if (!require("NetIndices")) {install.packages("NetIndices")}

shinyServer(function(input, output, session) {
  session$allowReconnect(TRUE)
  
  ##### render upload main input #####
  # output$mainInput <- renderUI({
  #   fileInput("mainInput","Upload input file:")
  # })
  
  output$fileUploaded <- reactive({
    return(!is.null(input$mainInput))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  ##### render removeNode.ui #####
  output$removeNode.ui <- renderUI({
    if(input$netType == "KO"){
      radioButtons(
        inputId="filterNodes",
        label="Remove un-annotated node(s):",
        choices=list("no","yes"),
        selected="no",
        inline=T
      )
    } else {
      radioButtons(
        inputId="filterNodes",
        label="Remove un-connected node(s):",
        choices=list("no","yes"),
        selected="no",
        inline=T
      )
    }
  })
  
  ##### render list of pathway names and IDs #####
  pathDf <- reactive({
    pathDf <- as.data.frame(read.table("data/pathways.list", quote = "", sep='\t',header=FALSE,check.names=FALSE,comment.char="",colClasses=c(rep("factor",4))))
    colnames(pathDf) <- c("pathID","pathName","pathGroup","pathType")
    return(pathDf)
  })
  
  output$pathID.ui <- renderUI({
    pathDf <- pathDf()
    pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
    selectInput('pathID', label = "Pathway ID:",
                choices = as.list(levels(pathDf$pathID)),
                selected = levels(pathDf$pathID)[93])
  })
  
  output$pathSelect <- renderUI({
    pathDf <- pathDf()
    pathDf$pathName <- factor(pathDf$pathName, levels = unique(pathDf$pathName))
    
    selectInput('pathName', label = "Pathway name:",
                choices = as.list(levels(pathDf$pathName)),
                selected = levels(pathDf$pathName)[93])
  })
  
  output$multiPathID.ui <- renderUI({
    pathDf <- pathDf()
    pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
    selectedIDs <- c(levels(pathDf$pathID)[9],levels(pathDf$pathID)[10],levels(pathDf$pathID)[11])
    
    selectInput('multiPathID', label = "Pathway ID:",
                choices = as.list(levels(pathDf$pathID)),
                selected = selectedIDs,
                multiple = TRUE)
  })
  
  output$multiPathSelect <- renderUI({
    pathDf <- pathDf()
    pathDf$pathName <- factor(pathDf$pathName, levels = unique(pathDf$pathName))
    selectedNames <- c(levels(pathDf$pathName)[9],levels(pathDf$pathName)[10],levels(pathDf$pathName)[11])

    selectInput('multiPathName', label = "Pathway name:",
                choices = as.list(levels(pathDf$pathName)),
                selected = selectedNames,
                selectize=FALSE,
                multiple = TRUE)
  })

  ### Update value(s) of selected pathID/multiPathID and pathName/multiPathName
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
  
  observe({
    if(!is.null(input$multiPathID)){
      pathDf <- pathDf()
      pathDf$pathName <- factor(pathDf$pathName, levels = unique(pathDf$pathName))
      selectedNames <- pathDf[pathDf$pathID %in% input$multiPathID,]$pathName
      
      updateSelectInput(
        session, 'multiPathName', label = "Pathway Name:",
        choices = as.list(levels(pathDf$pathName)),
        selected = selectedNames
      )
    }
  })
  
  observe({
    if(!is.null(input$multiPathName)){
      pathDf <- pathDf()
      pathDf$pathID <- factor(pathDf$pathID, levels = unique(pathDf$pathID))
      selectedIDs <- pathDf[pathDf$pathName %in% input$multiPathName,]$pathID
      
      updateSelectInput(
        session, 'multiPathID', label = "Pathway ID:",
        choices = as.list(levels(pathDf$pathID)),
        selected = selectedIDs
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
    if(input$netType == "CPD"){
      cxnFile <- suppressWarnings(paste0("data/keggcxncpd/",input$pathID,".cpd"))
    }
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
  
  networkDf <- function(pathID){
    # if(is.null(input$pathID)){return ()}
    # else {
      cxnFile <- suppressWarnings(paste0("data/keggcxn/",pathID,".cxn"))
      if(input$netType == "CPD"){
        cxnFile <- suppressWarnings(paste0("data/keggcxncpd/",pathID,".cpd"))
      }
      if(file.exists(cxnFile)){
        inputDf <- as.data.frame(read.table(cxnFile, sep='\t',header=F,check.names=FALSE,comment.char=""))
        if(input$netType == "KO"){
          colnames(inputDf) <- c("from","to","interaction")
        } else {
          colnames(inputDf) <- c("from","to","reaction","id")
        }
        
        return(inputDf)
      } else {
        return()
      }
    # }
  }
  
  ##### read input file #####
  annoDf <- reactive({
    filein <- input$mainInput
    
    if(is.null(filein)){ return()}
    annoDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE))
    colnames(annoDf)[2] <- "id"
    return (annoDf)
  })
  
  ##### render list of source species #####
  output$sourceList.ui <- renderUI({
    annoDf <- annoDf()
    if(!is.null(annoDf)){
      source <- c("all",unique(annoDf$source))
      selectInput('sourceList', label = "Choose source taxon:", choices = source, selected = "all", width = "80%")
    }
  })
  
  ##### get list of edges #####
  edgeData <- reactive({
    inputDf <- networkDf(input$pathID)
    
    if(input$netType == "CPD"){
      ### filter edges by annotated proteins
      annoDf <- annoDf()
      if(!is.null(annoDf)){
        inputDf$color <- "grey"
        inputDf$shadow <- FALSE
        inputDf$dashes <- TRUE
        
        ### filter anno data based on selected source species 
        if(!is.null(input$sourceList)){
          if(input$sourceList != "all"){
            annoDf <- annoDf[annoDf$source == input$sourceList,]
          }
        }
        
        ### remove reactions that have no enzymes (no annotated proteins)
        inputDf$color[inputDf$id %in% annoDf$id] <- "red"
        inputDf$shadow[inputDf$id %in% annoDf$id] <- TRUE
        inputDf$dashes[inputDf$id %in% annoDf$id] <- FALSE
        
        ### return
        edgeDf <- unique(inputDf[,c("from","to","color","shadow","dashes")])
        return(edgeDf)
      } else {
        edgeDf <- inputDf[,(1:2)]
        return(edgeDf)
      }
    } else {
      edgeDf <- inputDf[,(1:2)]
      return(edgeDf)
    }
  })
  
  ##### get list of nodes and mark annotated nodes from input file #####
  nodeDataPre <- function(pathID){
    inputDf <- networkDf(pathID)
    if(is.null(inputDf)){return()}
    
    ### join source and target KOs to get full list of nodes for reference network
    sourceDf <- data.frame("id" = inputDf$from,stringsAsFactors = FALSE)
    targetDf <- data.frame("id" = inputDf$to,stringsAsFactors = FALSE)
    nodeDf <- unique(rbind(sourceDf,targetDf))
    
    annoDf <- annoDf()
    if(!is.null(annoDf)){
      ### filter data based on selected source species
      if(!is.null(input$sourceList)){
        if(input$sourceList != "all"){
          annoDf <- annoDf[annoDf$source == input$sourceList,]
        }
      }
      
      if(input$netType == "KO"){
        ### map annotated nodes to reference net
        joinedDf <- merge(nodeDf,annoDf, by="id", all.x=TRUE)
        joinedDf$group <- joinedDf$source
        
        ### remove un-annotated nodes if requested
        if(input$filterNodes == "yes"){
          # joinedDf <- joinedDf[complete.cases(joinedDf),]
          joinedDf <- joinedDf[!is.na(joinedDf$geneID),]
        }
        
        return(joinedDf)
      } else {
        ### remove un-connected nodes if requested
        if(input$filterNodes == "yes"){
          inputDf$id <- as.character(inputDf$id)
          inputDf <- inputDf[inputDf$id %in% annoDf$id,]
          meltInputDf <- melt(inputDf[,c("from","to","id")], id="id")
          colnames(meltInputDf) <- c("ko","variable","id")
          
          outDf <- unique(data.frame("id" = meltInputDf$id,stringsAsFactors = FALSE))
          return(outDf)
        } else {
          return(nodeDf)
        }
      }
    } else {
      return(nodeDf)
    }
  }
  
  nodeData <- reactive({
    inputDf <- networkDf(input$pathID)
    
    if(is.null(inputDf)){return()}
    
    nodeDf <- nodeDataPre(input$pathID)
    
    ### read annotated data
    annoDf <- annoDf()
    if(!is.null(annoDf)){
      if(input$netType == "KO"){
        joinedDf <- nodeDf
        
        if(nrow(joinedDf) == 0){
          nodeEmpty <- data.frame("id"="1","label"="Selected source taxon does not have any annotated proteins","title"="EMPTY","color"="white")
          return(nodeEmpty)
        } else {
          ### change group type for reference nodes
          # if(length(unique(complete.cases(joinedDf))) > 1){
          if(nrow(joinedDf[is.na(joinedDf$geneID),]) > 0){
            joinedDf[is.na(joinedDf$geneID),]$group <- "reference"
          }
          
          ### process duplicated node IDs
          dupID <- joinedDf[duplicated(joinedDf$id),]$id
          if(length(dupID) > 0){
            # join gene IDs
            joinedDf$joinedGeneID <- joinedDf$geneID
            for(i in 1:length(dupID)){
              joinedDf[joinedDf$id == dupID[i],]$joinedGeneID <- toString(unique(sort(joinedDf[joinedDf$id == dupID[i],]$geneID)))
            }
            
            # join groups (source species)
            joinedDf$joinedGroup <- joinedDf$group
            for(i in 1:length(dupID)){
              joinedDf[joinedDf$id == dupID[i],]$joinedGroup <- toString(unique(sort(joinedDf[joinedDf$id == dupID[i],]$source)))
            }
            
            # get max FAS and min Dist for joined nodes
            aggrDist <- aggregate(joinedDf[,"patristicDist"],list(joinedDf$id),FUN="min")
            colnames(aggrDist) <- c("id","patristicDist")
            
            aggrFAS <- aggregate(joinedDf[,"FAS"],list(joinedDf$id),FUN="max")
            colnames(aggrFAS) <- c("id","FAS")
            
            aggrDf <- merge(aggrFAS, aggrDist, by="id")
            # for(i in 1:nrow(aggrFAS)){
            #   aggrDf$joinedColor[i] <- joinedDf[joinedDf$id %in% aggrDf$id[i] & joinedDf$FAS %in% aggrDf$FAS[i],]$color[1]
            # }
            
            # remove duplicated IDs from joinedDf and merge with aggrDf
            joinedDf <- joinedDf[!duplicated(joinedDf$id),]
            joinedDf <- merge(aggrDf, joinedDf, by="id", all.x = TRUE)
            
            # select needed columns of joinedDf
            joinedDf <- joinedDf[,c("id","joinedGroup","FAS.x","joinedGeneID","patristicDist.x")]
            colnames(joinedDf) <- c("id","group","FAS","geneID","patristicDist")
          } 
          # else {
          #   # select needed columns of joinedDf
          #   joinedDf <- joinedDf[,c("id","joinedGroup","FAS.x","joinedGeneID","patristicDist.x","color")]
          #   colnames(joinedDf) <- c("id","group","FAS","geneID","patristicDist","color")
          # }
          
          ### create node color based on patristic distance
          colorCodeDf = data.frame("colorLv" = c(0,1,2,3,4,5,6,7,8,9,10), 
                                   "color" = c('#AE250A','#A8390B','#A24E0C','#9C630D','#97780E','#918D0F','#8BA110','#86B611','#80CB12','#7AE013','#75F514'))
          joinedDf$colorLv <- 10 - as.integer(as.numeric(joinedDf$patristicDist)*10)
          joinedDf <- merge(joinedDf,colorCodeDf, by = "colorLv", all.x = TRUE)
          drop <- c("colorLv")
          joinedDf <- joinedDf[ , !(names(joinedDf) %in% drop)]
          
          joinedDf$color <- as.character(joinedDf$color)
          if(nrow(joinedDf[is.na(joinedDf$color),]) > 0){
            joinedDf[is.na(joinedDf$color),]$color <- "#878787"
          }
          
          ### node title, label, size(~FAS)
          joinedDf <- dplyr::rename(joinedDf, title = geneID)
          if(nrow(joinedDf[is.na(joinedDf$title),]) > 0){
            joinedDf[is.na(joinedDf$title),]$title <- as.character(joinedDf[is.na(joinedDf$title),]$id)
          }
          
          joinedDf$label <- joinedDf$id
          
          joinedDf$value <- as.integer(100^(joinedDf$FAS)*10)
          if(nrow(joinedDf[is.na(joinedDf$FAS),]) > 0){
            joinedDf[is.na(joinedDf$FAS),]$value <- 1
          }
          
          ### add node shapes
          shapeList <- c("circle", "triangle", "square", "dot", "star", "ellipse", "database", "text", "diamond")
          groupList <- sort(unique(joinedDf$group[joinedDf$group != "reference"]))
          n = 1
          for(i in 1:length(groupList)){
            if(n > 9){n = 1}
            joinedDf$shape[joinedDf$group == groupList[i]] <- as.character(shapeList[n])
            n=n+1
          }
          joinedDf$shape[joinedDf$group == "reference"] <- "box"
          
          ### add shadow
          joinedDf$shadow <- FALSE#TRUE
          joinedDf$shadow[joinedDf$group == "reference"] <- FALSE
          
          ### return node data
          subNodeDf <- joinedDf[,c("id","label","title","group","value","color","shape","shadow")]
          return(subNodeDf)
        }
      } else {
        nodeDf$color <- "#7aa8e7"
        nodeDf$title <- nodeDf$id
        nodeDf$label <- nodeDf$id
        nodeDf$shape <- "circle"
        nodeDf$shadow <- FALSE
        return(nodeDf)
      }
    } else {
      nodeDf$color <- "#878787"
      nodeDf$title <- nodeDf$id
      nodeDf$label <- nodeDf$id
      nodeDf$shape <- "box"
      nodeDf$shadow <- FALSE
      return(nodeDf)
    }
  })
  
  ##### RENDER NETWORK #####
  ### reset visOptions
  observeEvent(input$resetVisOption, {
    shinyjs::reset("performance")
    shinyjs::reset("maxSpeed")
    shinyjs::reset("layout")
  })
  
  ### create network
  networkPlot <- reactive({
    nodeData <- nodeData()
    edgeData <- edgeData()
    # print(nrow(nodeData))
    # print(nrow(edgeData))
    
    if(is.null(nodeData)){return()}
    
    network <- visNetwork(nodeData, edgeData, main = input$pathName) %>%
      # highlight nearest nodes
      visOptions(highlightNearest = list(enabled = TRUE, hover = FALSE),  nodesIdSelection = TRUE, collapse = TRUE) %>%
      # set fix randomSeed for getting the same network everytime
      visLayout(randomSeed = 123) %>%
      # physical parameters
      visPhysics(stabilization = TRUE, maxVelocity = input$maxSpeed) %>%
      visEvents(click = "function(nodes){
                Shiny.onInputChange('click', nodes.nodes[0]);
                ;}") %>% 
      # export
      visExport(type = "png", name = input$pathID, float = "top", style = "warning") %>% addExport(pdf = TRUE)
    
    return(network)
})
  
  ### create legend for mapped network
  networkPlotLegend <- reactive({
    nodeData <- nodeData()
    network <- networkPlot()

    if(nlevels(as.factor(nodeData$group)) > 1){
      # create legend info
      lnodes <- nodeData[,c("group","color","shape")]
      lnodes <- lnodes %>% dplyr::rename(label=group)
      # lnodes <- nodeData %>% select(group,color,shape) %>% dplyr::rename(label=group)
      lnodes <- lnodes[!duplicated(lnodes$label),]
      groupList <- as.character(levels(as.factor(lnodes$label)))
      for(i in 1:length(groupList)){
        lnodes$label[lnodes$label == groupList[i]] <- paste0("group_",i)
      }
      lnodes %>% distinct %>% mutate(title=label)
      lnodes <- lnodes[order(lnodes$label),]

      # create mapped network with legend
      networkLegend <- network %>% visLegend(position = "right",addNodes = lnodes,useGroups=F) %>% visOptions(selectedBy = "group", highlightNearest = list(enabled = TRUE, hover = FALSE),  nodesIdSelection = TRUE, collapse = TRUE)
      return(networkLegend)
      # return(network)
    } else {
      return(network)
    }
  })
  
  ### increase performance of the visualization
  networkPlotFinal <- reactive({
    nodeData <- nodeData()
    networkPlot <- networkPlotLegend()
    
    # disable visEdges smooth
    maxNode <- (1-input$performance)*100 + 30
    if(length(nodeData$id) >= maxNode){
      networkPlot <- networkPlot %>% visEdges(smooth = FALSE)# %>% visPhysics(stabilization = FALSE)
    }
    
    # use iGraph layout
    if(length(nodeData$id) > 500){
      networkPlot <- networkPlot %>% visIgraphLayout(layout = input$layout)# %>% visPhysics(stabilization = FALSE)
    }
    
    return(networkPlot)
  })
  
  ### output network
  output$plot <- renderVisNetwork({
    networkPlotFinal()
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
  
  ##### generate node property table #####
  igraphObj <- reactive({
    # get node and edge data
    nodeDf <- nodeData()
    edgeDf <- edgeData()
    edgeDf <- edgeDf[edgeDf$from %in% nodeDf$id & edgeDf$to %in% nodeDf$id,]
    
    # create igraph object
    graph <- graph_from_data_frame(edgeDf, directed = FALSE, vertices = nodeDf)
  })
  
  netProp <- reactive({
    graph <- igraphObj()
    # use GenInd function from NetIndices package to output network properties
    graph.adj<-get.adjacency(graph,sparse=FALSE)
    graph.properties<-GenInd(graph.adj)
    
    # get degree for all nodes
    all.deg.graph <- as.data.frame(degree(graph,v=V(graph),mode="all"))
    
    # network properties
    netProp <- data.frame(
      "Nodes" = graph.properties$N,  # number of nodes
      "Edges" = graph.properties$Ltot/2, # number of links
      "Avg_degree" = graph.properties$LD, # same as mean(all.deg.graph[,1]); link density (average # of links per node)
      "Max_degree" = max(all.deg.graph[,1]),
      "Avg_path_length" = average.path.length(graph, unconnected=TRUE),
      "Diameter" = diameter(graph)
    )
    return(netProp)
  })
  
  netPropMapped <- function(pathID){
    if(is.null(input$mainInput)){ return()}
    
    inputDf <- networkDf(pathID)
    pathDf <- pathDf()
    
    print(head(inputDf))
    
    networkProperty <- data.frame()
    if(input$netType == "KO"){
      ### read annotated data
      annoDf <- annoDf()
      if(!is.null(annoDf)){
        joinedDf <- nodeDataPre(pathID)
        
        if(nrow(joinedDf) == 0){
          nodeEmpty <- data.frame("id"="1","label"="Selected source taxon does not have any annotated proteins","title"="EMPTY","color"="white")
          return(nodeEmpty)
        } else {
          ### change group type for reference nodes
          # if(length(unique(complete.cases(joinedDf))) > 1){
          if(nrow(joinedDf[is.na(joinedDf$geneID),]) > 0){
            joinedDf[is.na(joinedDf$geneID),]$group <- "reference"
          }
          joinedDf <- joinedDf[!is.na(joinedDf$source),]
          for(source in levels(as.factor(joinedDf$source))){
            # node data
            nodeDf <- data.frame("id" = as.character(unique(joinedDf[joinedDf$source == source,"id"])))
            ### edge data
            edgeDf <- inputDf[,(1:2)]
            edgeDf <- edgeDf[edgeDf$from %in% nodeDf$id & edgeDf$to %in% nodeDf$id,]
            
            # create igraph object
            graph <- graph_from_data_frame(edgeDf, directed = FALSE, vertices = nodeDf)
            
            # use GenInd function from NetIndices package to output network properties
            graph.adj<-get.adjacency(graph,sparse=FALSE)
            graph.properties<-GenInd(graph.adj)
            
            # get degree for all nodes
            all.deg.graph <- as.data.frame(degree(graph,v=V(graph),mode="all"))
            
            # network properties
            netProp <- data.frame(
              "Nodes" = graph.properties$N,  # number of nodes
              "Edges" = graph.properties$Ltot/2, # number of links
              "Avg_Degree" = round(graph.properties$LD,2), # same as mean(all.deg.graph[,1]); link density (average # of links per node)
              "Max_Degree" = max(all.deg.graph[,1]),
              "Avg_Path_Len" = round(average.path.length(graph, unconnected=TRUE),2),
              "Diameter" = diameter(graph)
            )
            
            netProp$Pathway <- pathDf[pathDf$pathID == pathID,]$pathName
            netProp$Source <- source
            networkProperty <- rbind(networkProperty,netProp)
          }
          # return(networkProperty)
        }
      }
      networkProperty <- networkProperty[,c("Pathway","Source","Nodes","Edges","Avg_Degree","Max_Degree","Avg_Path_Len","Diameter")]
    } else {
      annoDf <- annoDf()
      if(!is.null(annoDf)){
        for(source in levels(as.factor(annoDf$source))){
          annoDfsub <- annoDf[annoDf$source == source,]
          inputDfsub <- inputDf[inputDf$id %in% annoDfsub$id,]
          
          if(nrow(inputDfsub) > 0){
            edgeDf <- unique(inputDfsub[,c(1:2)])
            
            sourceDf <- data.frame("id" = inputDfsub$from,stringsAsFactors = FALSE)
            targetDf <- data.frame("id" = inputDfsub$to,stringsAsFactors = FALSE)
            nodeDf <- unique(rbind(sourceDf,targetDf))
            
            # create igraph object
            graph <- graph_from_data_frame(edgeDf, directed = FALSE, vertices = nodeDf)
            
            # use GenInd function from NetIndices package to output network properties
            graph.adj<-get.adjacency(graph,sparse=FALSE)
            graph.properties<-GenInd(graph.adj)
            
            # get degree for all nodes
            all.deg.graph <- as.data.frame(degree(graph,v=V(graph),mode="all"))
            
            # network properties
            netProp <- data.frame(
              "Nodes" = graph.properties$N,  # number of nodes
              "Edges" = graph.properties$Ltot/2, # number of links
              "Avg_Degree" = round(graph.properties$LD,2), # same as mean(all.deg.graph[,1]); link density (average # of links per node)
              "Max_Degree" = max(all.deg.graph[,1]),
              "Avg_Path_Len" = round(average.path.length(graph, unconnected=TRUE),2),
              "Diameter" = diameter(graph)
            )
            
            netProp$Pathway <- pathDf[pathDf$pathID == pathID,]$pathName
            netProp$Source <- source
            networkProperty <- rbind(networkProperty,netProp)
          } else {
            netProp <- data.frame("Nodes" = 0, "Edges" = 0, "Avg_Degree" = 0, "Max_Degree" = 0, "Avg_Path_Len" = 0, "Diameter" = 0)
            netProp$Pathway <- pathDf[pathDf$pathID == pathID,]$pathName
            netProp$Source <- source
            networkProperty <- rbind(networkProperty,netProp)
          }
        }
        networkProperty <- networkProperty[,c("Pathway","Source","Nodes","Edges","Avg_Degree","Max_Degree","Avg_Path_Len","Diameter")]
      }
    }
    
    return(networkProperty)
  }
  
  multiNetProp <- reactive({
    pathIDs <- input$multiPathID
    multiNetProperty <- data.frame()
    for(pathID in pathIDs){
      print(pathID)
      netProp <- netPropMapped(pathID)
      print(netProp)
      multiNetProperty <- rbind(multiNetProperty,netProp)
    }
    
    # multiNetProperty <- netPropMapped(input$pathID)
    return(multiNetProperty)
  })
  
  nodeProp <- reactive({
    graph <- igraphObj()
    
    # use GenInd function from NetIndices package to output network properties
    graph.adj<-get.adjacency(graph,sparse=FALSE)
    graph.properties<-GenInd(graph.adj)
    
    # get degree for all nodes
    nodeDegree <- as.data.frame(degree(graph,v=V(graph),mode="all"))
    colnames(nodeDegree) <- c("Degree")
    nodeDegree$id = rownames(nodeDegree)
    nodeDegree <- nodeDegree[,c("id","Degree")]
    
    # get node information from input file
    annoDf <- annoDf()
    
    if(is.null(annoDf)){
      nodeDegree$href <- paste0("http://www.genome.jp/dbget-bin/www_bget?",nodeDegree$id)
      colnames(nodeDegree) <- c("Node","Degree","Link")
      return(nodeDegree)
    } else {
      mergeDf <- merge(nodeDegree,annoDf, by="id", all.x = TRUE)
      mergeDf$href <- paste0("http://www.genome.jp/dbget-bin/www_bget?",mergeDf$id)
      outDf <- mergeDf[,c("id","Degree","geneID","FAS","patristicDist","href")]
      colnames(outDf) <- c("Node","Degree","Gene ID","FAS","Patristic distance","Link")
      return(outDf)
    }
  })
  
  ##### network statistics #####
  ### network statistic tab
  output$stat.network.table <- renderTable({
    netProp()
  })
  
  output$stat.network.table.mapped <- renderDataTable({
    datatable(
      netPropMapped(input$pathID), 
      style='bootstrap', 
      options=list(pageLength=5), 
      rownames = FALSE)
  })
  
  output$stat.node.table <- renderDataTable({
    datatable(
      nodeProp(), 
      style='bootstrap', 
      options=list(pageLength=5), 
      rownames = FALSE)
  })
  
  output$download_stat.node.table <- downloadHandler(
    filename = function(){c("node_properties.txt")},
    content = function(file){
      write.table(nodeProp(),file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  degreeDistPlot <- reactive({
    if (is.null(nodeProp())){return ()}
    
    nodeProp <- nodeProp()
    labelMean <- paste0("mean = ",round(mean(nodeProp$Degree),2))
    
    if(input$distPlotType == "Density"){
      p <- ggplot(nodeProp, aes(x=Degree)) +
        geom_histogram(aes(y=..density..), binwidth=.5, alpha = .8, position="identity") +
        geom_density(alpha = .2, fill='#7ea4d6', color="#7ea4d6") +
        geom_vline(data=nodeProp, aes(xintercept=mean(nodeProp$Degree),colour="mean"),
                   linetype="solid", size=.5) +
        scale_color_manual(name = "", values = c(mean = "red"), label = labelMean) +
        theme_minimal()
      p <- p + theme(legend.title = element_blank(), legend.text = element_text(size=input$textLegend),
                     axis.text = element_text(size=input$textSize), axis.title = element_text(size=input$textSize))
      return(p)
    } else {
      p <- ggplot(nodeProp, aes(x=Degree)) +
        geom_histogram(binwidth=.5, alpha = .8, position="identity") +
        geom_vline(data=nodeProp, aes(xintercept=mean(nodeProp$Degree),colour="mean"),
                   linetype="solid", size=.5) +
        scale_color_manual(name = "", values = c(mean = "red"), label = labelMean) +
        theme_minimal()
      p <- p + theme(legend.title = element_blank(), legend.text = element_text(size=input$textLegend),
                     axis.text = element_text(size=input$textSize), axis.title = element_text(size=input$textSize))
      return(p)
    }
  })
  
  output$degreeDistPlot <- renderPlot(width = 512, height = 356,{
    degreeDistPlot()
  })
  
  output$degreePlot.ui <- renderUI({
    plotOutput("degreeDistPlot")
  })
  
  output$plotDownload_dist <- downloadHandler(
    filename = function() {paste0("distributionPlot.pdf")},
    content = function(file) {
      ggsave(file, plot = degreeDistPlot(), dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ### multi network analysis tab
  output$stat.multinetwork.table <- renderDataTable({
    datatable(
      multiNetProp(), 
      style='bootstrap', 
      options=list(pageLength=5), 
      rownames = FALSE)
  })
  
  ### properties statistics
  networkPropStat <- reactive({
    networkProperty <- multiNetProp()

    if(nrow(networkProperty) == 0){return()}
    
    meltedNetworkProp <- melt(networkProperty, id.vars = c("Pathway","Source"))
    colnames(meltedNetworkProp) <- c("Pathway","Source","Property","Value")
    
    meltedNetworkProp$Property <- as.character(meltedNetworkProp$Property)
    meltedNetworkProp$Property[meltedNetworkProp$Property == "Avg_Degree"] <- "Avg. degree"
    meltedNetworkProp$Property[meltedNetworkProp$Property == "Avg_Path_Len"] <- "Avg. path length"

    return(meltedNetworkProp)
  })
    
  ### for nodes and edges
  nodes_edges_plot <- reactive({
    meltedNetworkProp <- networkPropStat()
    if (is.null(meltedNetworkProp)){return ()}
    
    meltedNetworkPropSub2 <- meltedNetworkProp[meltedNetworkProp$Property %in% c("Nodes","Edges"),]
    meltedNetworkPropSub2$Property <- factor(meltedNetworkPropSub2$Property, levels = c("Nodes","Edges"))
    
    p <- ggplot(meltedNetworkPropSub2, aes(x=Pathway, y=Value,fill=Source)) +
      facet_wrap( ~ Property) +
      geom_point(aes(col=Source)) +
      coord_flip() +
      labs(y="Count") +
      theme(axis.title.x = element_text(size=input$textSizeNodes), axis.text.x = element_text(size=input$textSizeNodes),
            axis.title.y = element_text(size=input$textSizeNodes), axis.text.y = element_text(size=input$textSizeNodes),
            axis.ticks.x = element_blank(),
            strip.text.x = element_text(size = input$textSizeNodes),
            legend.position = "top", legend.text = element_text(size=input$textLegendNodes),legend.title = element_text(size=input$textLegendNodes)) +
      scale_fill_brewer(palette = "Set2")
    p
    
    return(p)
  })
  
  output$nodes_edges_plot <- renderPlot({
    nodes_edges_plot()
  })
  
  output$nodes_edges.ui <- renderUI({
    plotOutput("nodes_edges_plot",width = input$widthNodes, height = input$heightNodes)
  })
  
  output$plotDownload_nodes_edges <- downloadHandler(
    filename = function() {paste0("network_node_edge.pdf")},
    content = function(file) {
      ggsave(file, plot = nodes_edges_plot(), width = input$widthNodes*0.056458333, height = input$heightNodes*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ### for avg. degree, avg. path length and max path length (density)
  networkStatPlot <- reactive({
    meltedNetworkProp <- networkPropStat()
    if (is.null(meltedNetworkProp)){return ()}
    
    meltedNetworkPropSub <- meltedNetworkProp[!(meltedNetworkProp$Property %in% c("Nodes","Edges","Max_Degree")),]
    
    meltedNetworkPropSubSummary <- meltedNetworkPropSub %>%
      group_by(Property,Source) %>%
      summarize(mean = mean(Value))
    
    g <- ggplot(meltedNetworkPropSub, aes(x=Source, y=Value, fill=Source)) +
      facet_wrap( ~ Property) +
      geom_violin() + 
      geom_point(data = meltedNetworkPropSubSummary, aes(y = mean), color = "black", size = 2) +
      theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=input$textSize),
            axis.ticks.x = element_blank(),
            strip.text.x = element_text(size = input$textSize),
            legend.position = "top", legend.text = element_text(size=input$textLegend),legend.title = element_blank()) +
      scale_fill_brewer(palette = "Set2")
    
    return(g)
  })
  
  output$networkStatPlot <- renderPlot({
    networkStatPlot()
  })
  
  output$networkStat.ui <- renderUI({
    plotOutput("networkStatPlot",width = input$width, height = input$height)
  })
  
  output$plotDownload_networkStat <- downloadHandler(
    filename = function() {paste0("network_stat.pdf")},
    content = function(file) {
      ggsave(file, plot = networkStatPlot(), width = input$width*0.056458333, height = input$height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
})
