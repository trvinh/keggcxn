if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("shinyjs")) {install.packages("shinyjs")}
if (!require("rcytoscapejs")) {
  library("devtools");
  devtools::install_github("cytoscape/r-cytoscape.js");
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
      
      ########## INPUT TAB ###########
      # tabPanel(
      #   "Input & settings"
      # ),
      
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
              )
            )
          )
        ),
        
        column(
          8,
          uiOutput("msg"),
          rcytoscapejsOutput("plot", height="600px")
        ),
        
        column(4,
               h4("Download Network"),
               actionLink("saveImage", "Download as PNG"),
               h4("Clicked Node"),
               verbatimTextOutput("clickedNode"),
               h4("Connected Nodes"),
               verbatimTextOutput("connectedNodes"),
               hr(),
               h4("Clicked Node Data"),
               dataTableOutput("nodeDataTable"),
               hr(),
               h4("Data for Edges between Connected Nodes"),
               dataTableOutput("edgeDataTable")
        )
      )
    )
  )
)