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
              fileInput("mainInput","Annotated file:")
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
              radioButtons(
                inputId="filterNodes",
                label="Remove un-annotated node(s):",
                choices=list("no","yes"),
                selected="no",
                inline=T
              )
            ),
            column(
              2,
              bsButton("visOption","visOptions")
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
        h4("NODE PROPERTIES"),
        dataTableOutput("stat.table"),
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
          uiOutput("stat.ui")
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
      hr(),
      bsButton("resetVisOption","Default")
    )
  )
)