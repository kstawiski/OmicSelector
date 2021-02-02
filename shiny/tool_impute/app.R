library(shiny)
library(OmicSelector)
library(magick)
library("tools")
library(plotly)

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
  titlePanel("OmicSelector: Missing data analysis and imputation."),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The file should be prepared as for the analysis in OmicSelector. Features of interest should start with `hsa` prefix. The file should contain `Class` variable. See the exemplary file in the documentation."),
      fileInput("file2", "Upload data (*.csv or *.xlsx):",accept = c(".csv",".xlsx")),
      hr(),
      selectInput("method","Choose imputing method:",
                  c("Predictive mean matching" = "pmm",
                    "Replace with mean (for each variable, only for numeric variables)" = "mean"))

      
      
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      p("Input data preview (max 8 columns):"),
      dataTableOutput("inputprev"),
      hr(),
      p("Missing values in dataset: (plots are resizeable)"),
      jqui_resizable(plotOutput("missplot1")),
      jqui_resizable(plotOutput("missplot2")),
      jqui_resizable(plotOutput("missplot3")),
      p("PCA for complete cases before imputing:"),
      jqui_resizable(plotOutput("misspca1")),
      plotlyOutput("misspca1_3d"),
      hr(),
      p("Imputed dataset (without missing values, preview up to 8 columns):"),
      dataTableOutput("completedprev"),
      downloadButton("downloadData","Download imputed dataset (without missing values) for futher analysis"),
      hr(),
      p("PCA for all cases after imputation:"),
      jqui_resizable(plotOutput("misspca2")),
      plotlyOutput("misspca2_3d")
      
      
      
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
      max_col = 8
      if(ncol(dane) < 8) { max_col = ncol(dane) }
      output$inputprev = renderDataTable(dane[,1:max_col])
      danex = dplyr::select(dane, starts_with("hsa"))
      metadane = dplyr::select(dane, -starts_with("hsa"))
      
      
      # Missing plots
      library(ggplot2)
      output$missplot1 = renderPlot({ vis_miss(danex, warn_large_data = FALSE, cluster = T, sort_miss = T) + theme(plot.margin = margin(2, 2, 2, 2, "cm")) }) 
      output$missplot2 = renderPlot({ gg_miss_var(danex) })
      output$missplot3 = renderPlot({ if (n_var_miss(danex) > 2) { gg_miss_upset(danex, nsets = ncol(danex), nintersects = 1000) } })
      output$misspca1 = renderPlot({ OmicSelector_PCA(danex[complete.cases(danex),], metadane$Class[complete.cases(danex)]) })
      output$misspca1_3d = renderPlotly({ OmicSelector_PCA_3D(ttpm_features = danex[complete.cases(danex),], meta = metadane$Class[complete.cases(danex)]) })
      
      # Imputing
      tempdane = as.data.frame(cbind(`Class` = metadane$Class, danex))
      if(input$method == "pmm"){
        #tempdane = dataset
        suppressMessages(library(mice))
        temp1 = mice(tempdane, m=1)
        temp2 = temp1$data
        temp3 = mice::complete(temp1)
        completed = temp3
      } else if(input$method == "mean") {
        for(i in 2:ncol(tempdane)){
          tempdane[is.na(tempdane[,i]), i] <- mean(tempdane[,i], na.rm = TRUE)
        }
        completed = tempdane
      }
      
      
      
      
      # Assess results
      completedx = dplyr::select(completed, starts_with("hsa"))
      metadane2 = dplyr::select(metadane, -Class)
      completed_with_metadata = cbind(completed, metadane2)
      output$completedprev = renderDataTable(completed[,1:max_col])
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          data.table::fwrite(completed_with_metadata, file)
        }
      )
      output$misspca2 = renderPlot({ OmicSelector_PCA(completedx, completed$Class) })
      output$misspca2_3d = renderPlotly({ OmicSelector_PCA_3D(ttpm_features = completedx, meta = completed$Class) })
      
      
      
      
    }
    waiter_hide()
  })}

# Run the application 

shinyApp(ui = ui, server = server)