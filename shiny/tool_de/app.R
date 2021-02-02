library(shiny)
library(OmicSelector)
library(magick)
library("tools")
library(plotly)
library(summarytools)
library(DT)

suppressMessages(library(naniar))
suppressMessages(library(VIM))
library(shinyjqui)

options(shiny.maxRequestSize = 30*1024^2)

library(waiter)
waiting_screen <- tagList(
  spin_3(),
  h4("OmicSelector is working...")
) 

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  use_waiter(),
  # Application title
  titlePanel("OmicSelector: Differential analysis using corrected t-test."),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The file should be prepared as for the analysis in OmicSelector. Features of interest should start with `hsa` prefix. The file should contain `Class` variable. See the exemplary file in the documentation."),
      fileInput("file2", "Upload data (*.csv or *.xlsx):",accept = c(".csv",".xlsx")),
      hr(),
      textInput("caselabel","Label for cases (in `Class` variable):", "Cancer"),
      selectInput("mode","Type of values:",
                  c("Log-transformed (e.g. log(TPM))" = "logtpm",
                    "Crude (e.g. deltaCt)" = "deltact")),
      selectInput("adjust_p","Use adjusted p-value in volcano plot:",
                  c("No" = FALSE,
                    "Yes (adjusted with BH)" = TRUE)),
      uiOutput("only_label_ui")
      
      
      
      
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # p("Input data preview (max 100 rows and 8 columns):"),
      # dataTableOutput("inputprev"),
      # hr(),
      p("Principal Component Analysis:"),
      jqui_resizable(plotOutput("misspca1")),
      plotlyOutput("misspca1_3d", height = 500),
      hr(),
      p("Volcano plot:"),
      jqui_resizable(plotOutput("volc1")),
      hr(),
      p("Differential analysis results:"),
      dataTableOutput("DE"),
      downloadButton("downloadData",label = "Download DE results"),
      # hr(),
      # p("Descriptive analysis:"),
      # uiOutput("descriptiveanalysis")

      
      
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  
  
  
  observe({
    waiter_show(html = waiting_screen, color = "black")
    
    inFile2 <- input$file2
    if (!is.null(inFile2)) { 
      
      # CODE:
      
      
      file_path = inFile2$datapath
      # output$filepath = renderPrint(file_path)
      if(file_ext(file_path) == "xlsx") {
        library(xlsx)
        dane = xlsx::read.xlsx(file_path, 1)
      } else {
        library(data.table)
        dane = fread(file_path)
      }
      
      dane$Class = ifelse(dane$Class == input$caselabel, "Case", "Control")
      
      max_col = 8
      if(ncol(dane) < 8) { max_col = ncol(dane) }
      max_row = 100
      if(nrow(dane) < 100) { max_row = nrow(dane) }
      output$inputprev = renderDataTable(dane[1:max_row,1:max_col])
      danex = dplyr::select(dane, starts_with("hsa"))
      metadane = dplyr::select(dane, -starts_with("hsa"))
      
  
      output$misspca1 = renderPlot({ OmicSelector_PCA(danex, metadane$Class) })
      output$misspca1_3d = renderPlotly({ OmicSelector_PCA_3D(ttpm_features = danex, meta = metadane$Class) })
      
      # Descriptive
      # dfSum = dfSummary(dane)
      # output$descriptiveanalysis <- renderUI({
      #  print(dfSummary(dane), method = "render")
      # })
      
      # DE
      DE = OmicSelector_differential_expression_ttest(danex, classes = metadane$Class, mode = input$mode)
      output$DE = renderDataTable({ DE }, options = list(scrollX = TRUE))
      
      # Volcano
      if(is.null(input$only_label)) { output$only_label_ui = renderUI({ selectizeInput('only_label', 'Label features in plot:', choices = DE$miR, multiple = TRUE, selected = DE$miR[DE$`p-value` < 0.05]) }) }

      output$volc1 = renderPlot({ OmicSelector_vulcano_plot(selected_miRNAs = DE$miR, DE = DE, only_label = input$only_label, take_adjusted_p = input$adjust_p) }) 
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("de-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          data.table::fwrite(DE, file)
        }
      )
      
      
      
      
    }
    waiter_hide()
  })}

# Run the application 

shinyApp(ui = ui, server = server)