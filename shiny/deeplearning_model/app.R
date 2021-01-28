library(shiny)
library(OmicSelector)
library(magick)
OmicSelector::OmicSelector_load_extension("deeplearning")
options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("OmicSelector: Deep learning model viewer."),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            textInput("analysisid", "Analysis ID:", ""),
            textInput("modelid", "Model name:", ""),
            fileInput("file2", "or Upload model file (*.zip):",accept = c("application/zip",".zip")),
            hr(),
            p("Case inputs:"),
            uiOutput("caseinputs"),
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Model details", 
                                 p("Hyperparameters and results:"),
                                 tableOutput("wynikicsv"),
                                 p("ROC curve and cutoff:"),
                                 #uiOutput("ROCcurve"),
                                 plotOutput("ROCcurve", height = "500"),
                                 p("Model training log:"),
                                 verbatimTextOutput("summary"),
                                 p("Training and testing induction curve:"),
                                 plotOutput("trainingcurve", height = "500"),
                                 downloadButton("downloadModel", "Download model file"),
                            p(" ")),
                        tabPanel("Predict (single case)", p("Provided predictors:"),
                                 tableOutput("singlecasepredictors"),
                                 p("Prediction:"),
                                 tableOutput("singlecasepredictiontable"),
                                 p("Prediction details:"),
                                 verbatimTextOutput("singlecaseprediction"),
                                 
                                 ),
                        tabPanel("Predict (csv file)", 
                                 
                                 fileInput("file1", "Choose CSV File",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 p("Please upload new file .csv file. The file must contain variables consistent with predictors used in the network. If it contrains the binary variable 'Class' with 'Case' and 'Control' levels, additional performance metrics are calculated."),
                                 hr(),
                                 p("Data with predictions:"),
                                 downloadButton("downloadDataWithPredictions", "Download new data with predictions"),
                                 p("Performance and details:"),
                                 verbatimTextOutput("bulkpredictions"),
                                 )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observe({
        query <- parseQueryString(session$clientData$url_search)
        if (!is.null(query[['analysisid']])) {
            updateTextInput(session, "analysisid", value = query[['analysisid']])
        }
        if (!is.null(query[['modelid']])) {
            updateTextInput(session, "modelid", value = query[['modelid']])
        }
        
        
        
        file_path = paste0("/OmicSelector/", input$analysisid, "/models/deeplearning/", input$modelid, ".zip")
        # if(file.exists(file_path)) { toggle("file2") }
        manual = F
        inFile2 <- input$file2
        if (!is.null(inFile2)) { 
            file_path = inFile2$datapath
            manual = T }
        
        # file_path = "/OmicSelector/doktorat1/models/deeplearning/deeplearning_189-1611103892.zip"
        
        if(file.exists(file_path) || manual == T) {
            
            # Load model log
            file_in_zip = dplyr::filter(unzip(file_path, list = T), grepl("training.log",Name))[1,"Name"]
            unzip(file_path, file_in_zip, exdir = tempdir())
            file_path_unzipped = paste0(tempdir(), "/", file_in_zip)
            training_log = readLines(file_path_unzipped, warn = F)
            output$summary = renderText({ paste0(training_log, sep = "\n") })
            
            # Allow extracting model
            output$downloadModel <- downloadHandler(
                filename = function() {
                    paste(input$analysisid, "__", input$modelid, ".zip", sep = "")
                },
                content = function(file) {
                    file.copy(file_path, file)
                },
                contentType = "application/zip"
            )
            
            # Load model hyperparameter
            file_in_zip = dplyr::filter(unzip(file_path, list = T), grepl("wyniki.csv",Name))[1,"Name"]
            unzip(file_path, file_in_zip, exdir = tempdir())
            file_path_unzipped = paste0(tempdir(), "/", file_in_zip)
            hyperparams = read.csv(file_path_unzipped)
            hyperparams2 = data.table(`Parameter` = colnames(hyperparams), `Value` = as.character(hyperparams[1,]))
            output$wynikicsv = renderTable(hyperparams2)
            
            # Load ROC curve
            output$ROCcurve <- renderImage({
                file_in_zip = dplyr::filter(unzip(file_path, list = T), grepl("cutoff.png",Name))[1,"Name"]
                unzip(file_path, file_in_zip, exdir = tempdir())
                file_path_unzipped = paste0(tempdir(), "/", file_in_zip)
                width  <- session$clientData$output_ROCcurve_width
                if(width>500) { width = 500 }
                
                img = image_read(file_path_unzipped)
                img2 = image_scale(img, width)
                outfile <- tempfile(paste0(input$modelid, "roc"), fileext='.png')
                image_write(img2, outfile)
                list(src = outfile)
            }, deleteFile = T)
            
            # Load training log
            output$trainingcurve <- renderImage({
                file_in_zip = dplyr::filter(unzip(file_path, list = T), grepl("training.png",Name))[1,"Name"]
                unzip(file_path, file_in_zip, exdir = tempdir())
                file_path_unzipped = paste0(tempdir(), "/", file_in_zip)
                width  <- session$clientData$output_trainingcurve_width
                if(width>500) { width = 500 }
                
                img = image_read(file_path_unzipped)
                img2 = image_scale(img, width)
                outfile <- tempfile(paste0(input$modelid, "training"), fileext='.png')
                image_write(img2, outfile)
                list(src = outfile)
            }, deleteFile = T)
            
            # Single case inputs
            output$caseinputs <- renderUI({ })
            predict_csv = data.frame()
            x <- reactiveValuesToList(input)
            miRNAs = strsplit(hyperparams$formula, " \\+ ")[[1]]
            miRNAs = str_replace_all(miRNAs, "\n    ","")
            miRNAs_rev = rev(miRNAs)
            if(miRNAs[1] %in% names(x))
            {
                # tworzymy liste predyktorow
                
                for (i in 1:length(miRNAs)){
                    predict_csv[1,miRNAs[i]] = as.numeric(x[[miRNAs[i]]])
                }
                
                # predykcja
                output$singlecasepredictors = renderTable({ predict_csv })
                
                # predict_csv_path = tempfile("predict_csv", fileext = ".csv")
                # data.table::fwrite(predict_csv, predict_csv_path)
                singlecase_prediction = OmicSelector_deep_learning_predict(model_path = file_path, new_dataset = predict_csv, blinded = T, new_scaling = F)
                output$singlecaseprediction = renderPrint({
                    singlecase_prediction
                })
                output$singlecasepredictiontable = renderTable({ singlecase_prediction$predictions })
                
                
                
            } else {
                # jesli nie ma elementow z sieci
                for(i in 1:length(miRNAs_rev))
            {
                if(i == 1)
                {
                    insertUI(
                        selector = paste0("#caseinputs"),
                        where = "afterEnd",
                        ui = textInput(miRNAs_rev[i], miRNAs_rev[i], "0"))
                } else {
                    insertUI(
                        selector = paste0("#caseinputs"),
                        where = "afterEnd",
                        ui = textInput(miRNAs_rev[i], miRNAs_rev[i], "0"))
                }
                
            }
            }
            
            
            # Bulk prediction
            inFile <- input$file1
            
            if (is.null(inFile)) { return(NULL) } else {
                new_data = data.table::fread(inFile$datapath)
                if("Class" %in% colnames(new_data)) { 
                    blinded = F } else { blinded = T }
                
                if(all(miRNAs %in% colnames(new_data))) {
                bulk_prediction = OmicSelector_deep_learning_predict(model_path = file_path, new_dataset = new_data, blinded = blinded, new_scaling = F)
                with_predictions = cbind(bulk_prediction$predictions, new_data)
                # output$newdata = renderDataTable({ with_predictions })
                output$bulkpredictions = renderPrint({ bulk_prediction })
                output$downloadDataWithPredictions <- downloadHandler(
                    filename = function() {
                        paste("predictions.csv", sep = "")
                    },
                    content = function(file) {
                        data.table::fwrite(with_predictions, file)
                    }
                )} else {
                    output$bulkpredictions = renderPrint({ "Not all predictiors are present in the data." })
                }
            
                
            }
            
            
            
            
            
            
            
            
            
            

            # OmicSelector_deep_learning_predict(model_path = file_path,)
            
        } else {
            output$summary = renderText({ paste0("Model file located at ", file_path, " does not exist. Choose correct model name and analysis id.") })
        }
        
        
    })
    
    
    
}

# Run the application 

shinyApp(ui = ui, server = server)
