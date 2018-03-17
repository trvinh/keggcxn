if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("shinyjs")) {install.packages("shinyjs")}
if (!require("visNetwork")) {install.packages("visNetwork")}
if (!require("shinycssloaders")) {
  if("devtools" %in% installed.packages() == FALSE){
    install.packages("devtools")
  }
  devtools::install_github('andrewsali/shinycssloaders')
}

shinyUI(
  fluidPage(
    tags$style(type="text/css", "body {padding-top: 80px;}"),
    tags$head(tags$script(src="cyjs.js")),
    
    # Application title
    titlePanel(""),
    
    useShinyjs(),
    
    ################### main narvarpage tabs ##########################
    navbarPage(
      em(strong("KEGGcxn")),
      id ="tabs",
      collapsible = TRUE,
      inverse = TRUE,
      fluid = TRUE,
      position = "fixed-top",
      
      ########## NETWORK TAB ###########
      tabPanel(
        "Network",
        wellPanel(
          fluidRow(
            column(
              2,
              # uiOutput("mainInput")
              fileInput("mainInput","Annotated file:"),
              radioButtons(
                inputId="netType",
                label="Network type:",
                choices=list("KO","CPD"),
                selected="KO",
                inline=T
              )
            ),
            column(
              2,
              uiOutput("pathID.ui")
            ),
            column(
              3,
              uiOutput("pathSelect"),
              textOutput("pathDesc")
            ),
            column(
              3,
              uiOutput("sourceList.ui"),
              uiOutput("removeNode.ui")
            ),
            column(
              2,
              bsButton("visOption","visOptions"),
              hr(),
              bsButton("toPhyloProfile","Link to PhyloProfile")
            )
          )
        ),
        
        column(
          9,
          uiOutput("msg"),
          withSpinner(visNetworkOutput("plot",height = "600px")),
          conditionalPanel(
            condition = "output.fileUploaded == true",
            h5("Patristic distance to reference orthologs"),
            img(src='FAS_color_scale.png', width="60%")
          )
        ),
        
        column(
          3,
          h4("Selected Node"),
          verbatimTextOutput("clickedNode"),
          
          h4("Connected Nodes"),
          dataTableOutput("connectedNodes")
        )
      ),
      
      ########## STATISTIC TAB ###########
      tabPanel(
        "Network statistic",
        h4("NETWORK PROPERTIES"),
        withSpinner(tableOutput("stat.network.table"), proxy.height="50px", type=7, size = 0.5),
        hr(),
        conditionalPanel(
          condition = "output.fileUploaded == true",
          h4("MAPPED NETWORK PROPERTIES"),
          withSpinner(dataTableOutput("stat.network.table.mapped"), proxy.height="50px", type=7, size = 0.5),
          hr()
        ),
        h4("NODE PROPERTIES"),
        withSpinner(dataTableOutput("stat.node.table"), proxy.height="50px", type=7, size = 0.5),
        downloadButton('download_stat.node.table','Export table'),
        hr(),
        h4("NODE DEGREE DISTRIBUTION"),
        column(
          4,
          radioButtons(inputId="distPlotType", label="Choose type of distribution:", choices=list("Absolute count","Density"), selected="Absolute count", inline=T),
          numericInput("distTextSize","Label size",min=2,max=99,step=1,value=12,width=100),
          hr(),
          tags$head(
            tags$style(HTML('#plotDownload_dist{background-color:#A9E2F3}'))
          ),
          downloadButton('plotDownload_dist','Download plot')
        ),
        column(
          8,
          withSpinner(uiOutput("degreePlot.ui"), proxy.height="50px", type=7, size = 0.5)
        )
      ),
      
      ########## MULTIPLE NETWORKS ANALYSIS ###########
      tabPanel(
        "Multi network analysis",
        wellPanel(
        fluidRow(
            column(
              6,
              uiOutput("multiPathID.ui")
            ),
            column(
              6,
              uiOutput("multiPathSelect")
            )
          )
        ),
        conditionalPanel(
          condition = "output.fileUploaded == true",
          column(
            12,
            h4("NETWOK PROPERTY PLOTS")
          ),
          ### nodes vs edges plot
          column(
            2,numericInput("widthNodes","Width (px)",min=300,max=3200,step=50,value=800,width=100)
          ),
          column(
            2,numericInput("heightNodes","Height (px)",min=300,max=3200,step=50,value=300,width=100)
          ),
          column(
            2,numericInput("textSizeNodes","Text plot",min=2,max=99,step=1,value=12,width=100)
          ),
          column(
            2,numericInput("textLegendNodes","Text legend",min=2,max=99,step=1,value=12,width=100)
          ),
          column(
            4,
            strong("Download"),
            br(),
            tags$head(
              tags$style(HTML('#plotDownload_nodes_edges{background-color:#A9E2F3}'))
            ),
            downloadButton('plotDownload_nodes_edges','Download plot')
          ),
          column(
            12,
            withSpinner(uiOutput("nodes_edges.ui"), proxy.height="50px", type=7, size = 0.5),
            hr()
          ),
          
          ### degree and path length plot
          column(
            2,numericInput("width","Width (px)",min=300,max=3200,step=50,value=800,width=100)
          ),
          column(
            2,numericInput("height","Height (px)",min=300,max=3200,step=50,value=300,width=100)
          ),
          column(
            2,numericInput("textSize","Text plot",min=2,max=99,step=1,value=12,width=100)
          ),
          column(
            2,numericInput("textLegend","Text legend",min=2,max=99,step=1,value=12,width=100)
          ),
          column(
            4,
            strong("Download"),
            br(),
            tags$head(
              tags$style(HTML('#plotDownload_networkStat{background-color:#A9E2F3}'))
            ),
            downloadButton('plotDownload_networkStat','Download plot')
          ),
          column(
            12,
            withSpinner(uiOutput("networkStat.ui"), proxy.height="50px", type=7, size = 0.5),
            hr()
          ),
          
          ### network property table
          column(
            12,
            h4("NETWORK PROPERTIES"),
            withSpinner(dataTableOutput("stat.multinetwork.table"), proxy.height="50px", type=7, size = 0.5),
            hr()
          )
        )
      )
    ),
    
    ####### popup to confirm parsing data from input file
    bsModal(
      "visOptionBs", "visOptions", "visOption", size = "small",
      sliderInput("performance",
                  "Performance:", min = 0, max = 1, step = 0.1, value = 1, width = 200),
      sliderInput("maxSpeed",
                  "Max velocity:", min = 0, max = 1000, step = 10, value = 10, width = 200),
      selectInput('layout', label = "iGraph layout:",
                  choices = c("layout_nicely","layout.fruchterman.reingold","layout.kamada.kawai","layout.lgl","layout.mds","layout.reingold.tilford","layout.sphere","layout.star"),
                  selected = "layout_nicely"),
      hr(),
      bsButton("resetVisOption","Default")
    )
  )
)