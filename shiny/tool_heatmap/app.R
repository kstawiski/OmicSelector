library(shiny)
library(OmicSelector)
library(magick)
library("tools")
library(plotly)


options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("OmicSelector: Heatmap generator."),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The file should be prepared as for the analysis in OmicSelector. Features of interest should start with `hsa` prefix. The file should contain `Batch` and `Class` variables. See the exemplary file in the documentation."),
      fileInput("file2", "Upload data (*.csv or *.xlsx):",accept = c(".csv",".xlsx")),
      hr(),
      selectInput("labels", "Label matadata:",
                  c("Class" = "class",
                    "Batch" = "batch",
                    "Class+Batch" = "class+batch")),
      selectInput("zscored", "Z-score before plotting:",
                  c("No" = FALSE,
                    "Yes" = TRUE)),
      textInput("expression_name", "Expression metric:", "log10(TPM)"),
      selectInput("legend_pos", "Position of legend:",
                  c("Top right" = "topright",
                    "Top left" = "topleft",
                    "Bottom right" = "bottomright",
                    "Bottom left" = "bottomleft",
                    "Center" = "center")),
      sliderInput("legend_cex", "Legend size:",
                  min = 0.1, max = 5, value = 0.8),
      textInput("centered_on", "Centered on: (leave empty for median)", ""),
      textInput("trim_min", "Trim minimum: (leave empty for actual minimum)", ""),
      textInput("trim_max", "Trim maximum: (leave empty for actual maxium)", ""),
      sliderInput("height", "Plot height:",
                  min = 100, max = 2000, value = 600),
      sliderInput("margin", "Plot margin:",
                  min = 1, max = 100, value = 10),
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # p("Input data preview:"),
      # dataTableOutput("inputprev"),
      # hr(),
      plotOutput("heatmap")

      
      
      
      
      
      
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
      # output$inputprev = renderDataTable(dane)
      danex = dplyr::select(dane, starts_with("hsa"))
      metadane = dplyr::select(dane, -starts_with("hsa"))
      
      rlab = data.frame()
      if(input$labels == "class") { rlab = data.frame(Class = dane$Class) }
      if(input$labels == "batch") { rlab = data.frame(Batch = dane$Batch) }
      if(input$labels == "class+batch") { rlab = data.frame(Batch = dane$Batch, Class = dane$Class) }
      
      trim_min = NULL
      trim_max = NULL
      centered_on = NULL
      if(grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", input$trim_min)) { trim_min = as.numeric(input$trim_min) }
      if(grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", input$trim_max)) { trim_max = as.numeric(input$trim_max) }
      if(grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", input$centered_on)) { centered_on = as.numeric(input$centered_on) }
      
      
      output$heatmap = renderPlot({
        OmicSelector_heatmap(
          x = danex,
          rlab = rlab,
          zscore = input$zscored,
          margins = c(input$margin, input$margin),
          expression_name = input$expression_name,
          trim_min = trim_min,
          trim_max = trim_max,
          centered_on = centered_on,
          legend_pos = input$legend_pos,
          legend_cex = input$legend_cex
        )
      }, height = input$height)
      
      
      
      
      
      
      
      
    }
  })}

# Run the application 

shinyApp(ui = ui, server = server)