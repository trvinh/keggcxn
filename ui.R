if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("shinyjs")) {install.packages("shinyjs")}

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
              uiOutput("mainInput")
            ),
            column(
              2,
              uiOutput("pathID.ui")
            ),
            column(
              4,
              uiOutput("pathSelect"),
              textOutput("pathDesc")
            ),
            column(
              4,
              radioButtons(
                inputId="filterNodes",
                label="Remove un-annotated node(s):",
                choices=list("no","yes"),
                selected="no",
                inline=T
              ),
              sliderInput("maxSpeed",
                          "Max velocity:", min = 0, max = 500, step = 5, value = 10, width = 200)
            )
          )
        ),
        
        column(
          8,
          uiOutput("msg"),
          visNetworkOutput("plot",height = "600px")
        ),
        
        column(
          4,
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
    )
  )
)