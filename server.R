if (!require("DT")) {install.packages("DT")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("visNetwork")) {install.packages("visNetwork")}

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
  
  ##### render list of source species #####
  output$sourceList.ui <- renderUI({
    annoDf <- annoDf()
    if(!is.null(annoDf)){
      source <- c("all",unique(annoDf$source))
      selectInput('sourceList', label = "Choose source taxon:", choices = source, selected = "all", width = "80%")
    }
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
    # nodeDf$group <- "reference"
    
    ### read annotated data
    annoDf <- annoDf()
    if(!is.null(annoDf)){
      ### filter data based on selected source species
      if(!is.null(input$sourceList)){
        if(input$sourceList != "all"){
          annoDf <- annoDf[annoDf$source == input$sourceList,]
        }
      }
      
      ### map annotated nodes to reference net
      joinedDf <- merge(nodeDf,annoDf, by="id", all.x=TRUE)
      joinedDf$group <- joinedDf$source
      
      ### remove un-annotated nodes if requested
      if(input$filterNodes == "yes"){
        # joinedDf <- joinedDf[complete.cases(joinedDf),]
        joinedDf <- joinedDf[!is.na(joinedDf$geneID),]
      }
      
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
        joinedDf$shadow <- TRUE
        joinedDf$shadow[joinedDf$group == "reference"] <- FALSE
        
        ### return node data
        subNodeDf <- joinedDf[,c("id","label","title","group","value","color","shape","shadow")]
        return(subNodeDf)
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
  
  ##### get list of edges #####
  edgeData <- reactive({
    inputDf <- networkDf()
    edgeDf <- inputDf[,(1:2)]
    return(edgeDf)
  })
  
  ##### RENDER NETWORK #####
  ### reset visOptions
  observeEvent(input$resetVisOption, {
    shinyjs::reset("performance")
    shinyjs::reset("maxSpeed")
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
      visExport(type = "pdf", name = input$pathID, float = "top", style = "warning") %>% addExport(pdf = TRUE)
    
    return(network)
  })
  
  ### create legend for mapped network
  networkPlotLegend <- reactive({
    nodeData <- nodeData()
    network <- networkPlot()
    
    if(nlevels(as.factor(nodeData$group)) > 1){
      # create legend info
      lnodes <- nodeData %>% select(group,color,shape) %>% dplyr::rename(label=group)
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
      networkPlot <- networkPlot %>% visIgraphLayout()# %>% visPhysics(stabilization = FALSE)
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
  nodeProp <- reactive({
    # get node and edge data
    nodeDf <- nodeData()
    edgeDf <- edgeData()
    edgeDf <- edgeDf[edgeDf$from %in% nodeDf$id & edgeDf$to %in% nodeDf$id,]

    # calculate node degrees
    nodeDegree <- data.frame("id" = nodeDf$id)
    
    allConnect <- c(as.character(edgeDf$from), as.character(edgeDf$to))
    if(length(allConnect) > 0){
      degree <- as.data.frame(table(allConnect))
      colnames(degree) <- c("id","Degree")
      nodeDegree <- merge(nodeDegree, degree, by="id", all.x = TRUE)
    } else {
      nodeDegree$Degree <- 0
    }

    if(nrow(nodeDegree[is.na(nodeDegree$Degree),]) > 0){
      nodeDegree[is.na(nodeDegree$Degree),]$Degree <- 0
    }
    
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
  output$stat.table <- renderDataTable({
    datatable(nodeProp(), style='bootstrap', options=list(pageLength=5, sDom  = '<"top">lrt<"bottom">ip'), filter = "bottom")
  })
  
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
      p <- p + theme(legend.title = element_blank(), legend.text = element_text(size=input$distTextSize),
                     axis.text = element_text(size=input$distTextSize), axis.title = element_text(size=input$distTextSize))
      return(p)
    } else {
      p <- ggplot(nodeProp, aes(x=Degree)) +
        geom_histogram(binwidth=.5, alpha = .8, position="identity") +
        geom_vline(data=nodeProp, aes(xintercept=mean(nodeProp$Degree),colour="mean"),
                   linetype="solid", size=.5) +
        scale_color_manual(name = "", values = c(mean = "red"), label = labelMean) +
        theme_minimal()
      p <- p + theme(legend.title = element_blank(), legend.text = element_text(size=input$distTextSize),
                     axis.text = element_text(size=input$distTextSize), axis.title = element_text(size=input$distTextSize))
      return(p)
    }
  })
  
  output$degreeDistPlot <- renderPlot(width = 512, height = 356,{
    degreeDistPlot()
  })
  
  output$stat.ui <- renderUI({
    plotOutput("degreeDistPlot")
  })
  
  output$plotDownload_dist <- downloadHandler(
    filename = function() {paste0("distributionPlot.pdf")},
    content = function(file) {
      ggsave(file, plot = degreeDistPlot(), dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
})
