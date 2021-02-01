library(shiny)
library(OmicSelector)
library(magick)
library("tools")
library(plotly)


options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("OmicSelector: Batch-effect correction."),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The file should be prepared as for the analysis in OmicSelector. Features of interest should start with `hsa` prefix. The file should contain `Batch` and `Class` variables. See the exemplary file in the documentation."),
      fileInput("file2", "Upload data (*.csv or *.xlsx):",accept = c(".csv",".xlsx")),
      hr(),
      selectInput("model", "Batch correction mode (covariates):",
                  c("~ Batch" = "~ Batch",
                    "~ Batch + Class" = "~ Batch + Class")),
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      p("Input data preview:"),
      dataTableOutput("inputprev"),
      hr(),
      p("Initial PCA:"),
      plotlyOutput("init_pca3d", height = 500),
      plotOutput("init_pca1", height = 500),
      plotOutput("init_pca2", height = 500),
      hr(),
      p("Corrected dataset:"),
      dataTableOutput("corrected"),
      downloadButton("downloadData", "Download batch-corrected data"),
      hr(),
      p("Corrected PCA:"),
      plotlyOutput("corrected_pca3d", height = 500),
      plotOutput("corrected_pca1", height = 500),
      plotOutput("corrected_pca2", height = 500),

      
      
      
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
      dane$meta = paste0(dane$Batch, " - ", dane$Class)
      output$init_pca3d = renderPlotly({ OmicSelector_PCA_3D(ttpm_features = danex, meta = dane$meta) })
      output$init_pca1 = renderPlot({ OmicSelector_PCA(ttpm_features = danex, meta = dane$Class) })
      output$init_pca2 = renderPlot({ OmicSelector_PCA(ttpm_features = danex, meta = dane$Batch) })
      
      corrected = OmicSelector_combat(
        danex,
        metadane = metadane,
        model = input$model
      )
      output$corrected = renderDataTable(corrected)
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          data.table::fwrite(corrected, file)
        }
      )
      
      danex2 = dplyr::select(corrected, starts_with("hsa"))
      metadane2 = dplyr::select(corrected, -starts_with("hsa"))
      dane$meta = paste0(dane$Batch, " - ", dane$Class)
      output$corrected_pca3d = renderPlotly({ OmicSelector_PCA_3D(ttpm_features = danex2, meta = dane$meta) })
      output$corrected_pca1 = renderPlot({ OmicSelector_PCA(ttpm_features = danex2, meta = dane$Class) })
      output$corrected_pca2 = renderPlot({ OmicSelector_PCA(ttpm_features = danex2, meta = dane$Batch) })
      
      
      
      
      
      
      
    }
})}

# Run the application 

shinyApp(ui = ui, server = server)
