library(shiny)
library(OmicSelector)
library(magick)
library("tools")
library(plotly)


options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("OmicSelector: Tool templete."),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The file should be prepared as for the analysis in OmicSelector. Features of interest should start with `hsa` prefix. The file should contain `Batch` and `Class` variables. See the exemplary file in the documentation."),
      fileInput("file2", "Upload data (*.csv or *.xlsx):",accept = c(".csv",".xlsx")),
      hr(),

    
      
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      p("Input data preview:"),
      dataTableOutput("inputprev"),
      hr(),
      
      
      
      
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  observe({
    hideTab("tabs", "Input data")
    hideTab("tabs", "Results (output)")
    
    inFile2 <- input$file2
    if (!is.null(inFile2)) { 
      showTab("tabs", "Input data")
      showTab("tabs", "Results (output)")
      file_path = inFile2$datapath
      # output$filepath = renderPrint(file_path)
      if(file_ext(file_path) == "xlsx") {
        library(xlsx)
        dane = xlsx::read.xlsx(file_path, 1)
      } else {
        library(data.table)
        dane = fread(file_path)
      }
      output$inputprev = renderDataTable(dane)
      danex = dplyr::select(dane, starts_with("hsa"))
      metadane = dplyr::select(dane, -starts_with("hsa"))
      
      
      
      
      
      
      
      
      
      
    }
  })}

# Run the application 

shinyApp(ui = ui, server = server)