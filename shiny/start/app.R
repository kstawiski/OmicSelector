##### SETUP
library(shiny)
library(shinyjs)
library(OmicSelector)
library(dplyr)
library(tools)
library(shinyfullscreen)
suppressMessages(suppressMessages(library(data.table)))
suppressMessages(suppressMessages(library(dplyr)))
suppressMessages(suppressMessages(library(ggplot2)))
suppressMessages(suppressMessages(library(OmicSelector)))
oldwd = getwd()
options(shiny.maxRequestSize = 30*1024^2)
dummy_df = as.data.frame(matrix(rnorm(100),nrow = 10))

library(waiter)
waiting_screen <- tagList(
    spin_3(),
    h4("OmicSelector is working...")
)

##### INIT VARS
if_analysis_id_duplicated = TRUE
while(if_analysis_id_duplicated == TRUE){
    init_random_analysis_id = stringi::stri_rand_strings(1, 15)
    if_analysis_id_duplicated = dir.exists(paste0("/OmicSelector/",init_random_analysis_id))
}

##### UI
ui <- fluidPage(
    #### STYLE
    useShinyjs(),
    tags$head(
        tags$style(HTML("
      #error_msg {
        color: red;
      }
    "))
        
    ),
    
    ##### MAIN
    htmlTemplate("ui.html",
                 init_random_analysis_id = init_random_analysis_id)
)

##### SERVER
server <- function(input,output,session){
    observeEvent(input$button_startover, {
        shinyjs::toggle(id= "panelA")
        shinyjs::toggle(id= "panelB")
    })
  
  observeEvent(input$select_plot, { 
      output$wykres = renderImage({
        filename <- input$select_plot
        list(src = filename, width = "100%", height = "100%")
      })
    })
  
  
  observeEvent(input$select_df, { 
    output$download_df = downloadHandler(filename = input$select_df, content = function(file) { file.copy(input$select_df, file) })
  }) 
  
    
  observeEvent(input$button_start_analysis, {

        # VALIDATION:
        walidacja = ""
        if(dir.exists(paste0("/OmicSelector/",input$analysis_id)) == TRUE) { walidacja = "The analysis with the same analysis ID already exists. Choose a unique one." }
        if(is.null(input$file2)) { walidacja = "You have to provide xlsx or csv file for the analysis." }
        try({ 
          plik = input$file2$datapath
          if(tools::file_ext(plik) == "xlsx") { dane = xlsx::read.xlsx(plik,sheetIndex=1) } else {
          dane = data.table::fread(plik)
        } })
        if(is.null(dane)) { walidacja = "You have to provide valid xlsx or csv file for the analysis." }
        
        if(walidacja != "") { showNotification(walidacja, duration = 10, type = "error") }
        validate(need(walidacja == "", "Form validation error", "button_start_analysis"))
        
        
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Validating submitted data & performing initial analysis", value = 0)
        
        # UNIQUEID:
        uniquid = gsub("\\.","",paste0(as.numeric(Sys.time())))
        dir.create(paste0(oldwd, "/", uniquid))
        setwd(paste0(oldwd, "/", uniquid))
        OmicSelector_log(paste0("Staring pre-check of OmicSelector analysis ID: ", getwd()))
        
        # ACTION:
        n_steps <- 13
        task_status = "OK"
      
        if(file.exists("task.log")) { file.remove("task.log") }

        OmicSelector_log("Welcome! OmicSelector id loaded.");
        options(warn=-1)
        progress$inc(1/n_steps) # PROGRESS BAR 1
        
        writeLines("ERROR", "var_initcheck.txt", sep="")
        
        writeLines(input$class_interest, "var_class_interest.txt", sep="")
        writeLines(input$split, "var_split.txt", sep="")
        writeLines(input$type, "var_type.txt", sep="")
        
        colnames(dane) = make.names(colnames(dane), unique = T)
        dane[dane == ""] = NA
        
        error = FALSE
        if("Class" %in% colnames(dane)) { OmicSelector_log("✓ The file contains Class variable. ") } else 
        { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); OmicSelector_log("☒ The file contains DOES NOT Class variable. "); task_status = "ERROR"; showNotification(type = "error", ui = "☒ The file contains DOES NOT Class variable. ") }
        
        
        klasy = unique(dane$Class)
        ref = c("Case","Control")
        if(!dplyr::setequal(ref, klasy)) {
          class_interest = readLines("var_class_interest.txt", warn = F)
          dane$ClassOrginal = dane$Class
          dane$Class = ifelse(dane$Class == class_interest, "Case", "Control")
          OmicSelector_log("\n✓ The Class variable was converted successfully. The orginial values are saved in ClassOrginal variable.")
        } else {
          OmicSelector_log("\n✓ The Class variable has only case and control cases.")
        }
        progress$inc(1/n_steps)
        
        dane$Class = factor(dane$Class, levels = c("Control","Case"))
        if(table(dane$Class)[1] == 0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); OmicSelector_log("☒ There are no control cases."); task_status = "ERROR"; showNotification(type = "error", ui = "☒ There are no control cases.") }
        if(table(dane$Class)[2] == 0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); OmicSelector_log("☒ There are no cases of interest."); task_status = "ERROR"; showNotification(type = "error", ui = "☒ There are no cases of interest.") }
        OmicSelector_log(paste0("\n✓ The data contains ", table(dane$Class)[1], " `Control` cases and ", table(dane$Class)[2], " `Case` cases (cases of interest)."))
        progress$inc(1/n_steps)
        
        
        try({ 
        temp = dplyr::select(dane, starts_with("hsa"))
        if(ncol(temp)==0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); OmicSelector_log("☒ The data does not contain any features (e.g. miRNAs) for feautre selection. Remember that feature names should start from hsa..."); task_status = "ERROR"; showNotification(type = "error", ui = "☒ The data does not contain any features (e.g. miRNAs) for feautre selection. Remember that feature names should start from hsa...") }
        OmicSelector_log(paste0("\n✓ The data contains ", ncol(temp), " features (e.g. miRNAs) for selection."))
        if(ncol(temp) > 0) { go_next = T } else { go_next = F }
        })
        progress$inc(1/n_steps)
        
        
        try({
        czy_numeryczne = sapply(temp, is.numeric)
        if(sum(czy_numeryczne) != ncol(temp)) { writeLines(as.character("FAIL", "var_initcheck.txt", sep = "")); OmicSelector_log("☒ Some of the features are not numeric. Please remove them. Not numeric: ", paste0(colnames(temp)[czy_numeryczne == F], collapse = ", ")); task_status = "ERROR"; showNotification(type = "error", ui = "☒ Some of the features are not numeric. Please remove them. Not numeric: ", paste0(colnames(temp)[czy_numeryczne == F], collapse = ", "))}
        OmicSelector_log(paste0("\n✓ All features are numeric."))
        })
        progress$inc(1/n_steps)
        
        try({
        missing = FALSE
        czy_brakna = sapply(temp, is.na)
        if(sum(colSums(czy_brakna)) != 0) 
        { OmicSelector_log("☒ Some of the features contain missing data.\nWith missing values:\n  -", paste0(colnames(temp)[colSums(czy_brakna) > 0], collapse = "\n  - ")); task_status = "ERROR"; showNotification(type = "error", ui = "☒ Some of the features contain missing data.\nWith missing values:\n  -", paste0(colnames(temp)[colSums(czy_brakna) > 0], collapse = "\n  - "))
          missing = T } else 
          { OmicSelector_log(paste0("\n✓ There are no missing data in features.")) }
        })
        progress$inc(1/n_steps)
        
        try({
        batch = F
        if("Batch" %in% colnames(dane)) { OmicSelector_log("\n✓ The file contains `Batch` variable that can be used for batch-effect correction. You can apply combat-based batch correction via OmicSelector, but not via GUI. Please see the package manual. Please also check if the following contingency table is correct:\n ")
          print(table(dane$Class, dane$Batch))
          batch = T } else { OmicSelector_log("\n✓ The file does not contain `Batch` variable that can be used for batch-effect correction.") }
        })
        progress$inc(1/n_steps)
        
        try({
        x = dplyr::select(dane, starts_with("hsa"))
        like_counts = sapply(x, function(x2) (sum(unlist(na.omit(x2))%%1 == 0) + sum(unlist(na.omit(x2)) >= 0))/(2*length(unlist(na.omit(x2)))) == 1)
        positive = F
        if(sum(like_counts)/ncol(x)) { OmicSelector_log("\n✓ Feature values are positive integers. The file could represent read counts. Please note that the feature selection pipeline requires normalized data. You can use OmicSelector to convert counts to tpm, but not via GUI. Please refer to package manual."); positive = T; } else {
          OmicSelector_log("\n✓ Feature values are not positive integers. This is ok if your input data is normalized (e.g. deltaCt values or tpm-normalized counts). Features not looking like counts: ");
          OmicSelector_log(paste0(colnames(x)[which(like_counts == FALSE)], collapse = ", "))
        }
        writeLines(as.character(positive), "var_seemslikecounts.txt", sep="")
        })
        progress$inc(1/n_steps)
        
        suppressMessages(library(OmicSelector))
        try({
        split = readLines("var_split.txt", warn = F)
        if(split == "yes") {
          if("mix" %in% colnames(dane)) {
            if(sum(which(dane$mix == 'train')) > 0) { 
              OmicSelector_log(paste0("\n✓ Samples in training set: ", sum(dane$mix == 'train')));
            } else {
              OmicSelector_log(paste0("\n✓ Samples in training set: ", sum(which(dane$mix == 'train'))));
              OmicSelector_log("There are no samples in training set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups."); task_status = "ERROR"; showNotification(type = "error", ui = "There are no samples in training set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
            }
            if(sum(which(dane$mix == 'test')) > 0) { 
              OmicSelector_log(paste0("\n✓ Samples in test set: ", sum(dane$mix == 'test')));
            } else {
              OmicSelector_log(paste0("\n✓ Samples in test set: ", sum(dane$mix == 'test')));
              OmicSelector_log("There are no samples in test set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups."); task_status = "ERROR"; showNotification(type = "error", ui = "There are no samples in test set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
            }
            if(sum(which(dane$mix == 'valid')) > 0) { 
              OmicSelector_log(paste0("\n✓ Samples in validation set: ", sum(dane$mix == 'valid')));
            } else {
              OmicSelector_log(paste0("\n✓ Samples in validation set: ", sum(dane$mix == 'valid')));
              OmicSelector_log("There are no samples in validation set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups."); task_status = "ERROR"; showNotification(type = "error", ui = "There are no samples in validation set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
            }
            train = dplyr::filter(dane, mix == "train")
            fwrite(train, "mixed_train.csv")
            test = dplyr::filter(dane, mix == "test")
            fwrite(test, "mixed_test.csv")
            valid = dplyr::filter(dane, mix == "valid")
            fwrite(valid, "mixed_valid.csv")
            mixed = rbind(train,test,valid)
            fwrite(mixed, "data_start.csv")
            fwrite(mixed, "mixed.csv")
            
            if(sum(which(dane$mix == 'train_balanced')) > 0) {
              train_balanced = dplyr::filter(dane, mix == "train_balanced")
              fwrite(train_balanced, "mixed_train_balanced.csv")
              OmicSelector_log(paste0("\n✓ Balanced training set file was retored."))
              
            } 
            
            dane = OmicSelector_load_datamix(use_smote_not_rose = T)  # load mixed_*.csv files
            train = dane[[1]]
            test = dane[[2]]
            valid = dane[[3]]
            train_smoted = dane[[4]]
            train$mix = "train"
            test$mix = "test"
            valid$mix = "valid"
            train_smoted$mix = "train_balanced"
            merged = rbind(train,test,valid, train_smoted)
            fwrite(merged, "merged.csv")
            
          } else {
            OmicSelector_log("\n✓ The data is not splitted, i.e. doesn't have 'train', 'test' and 'valid' in 'mix' variable. We will perform data splitting (60% train, 20% test, 20% valid). ");
            metadane = dplyr::select(dane, -starts_with("hsa"))
            ttpm_features = dplyr::select(dane, starts_with("hsa"))
            mixed = OmicSelector_prepare_split(metadane = metadane, ttpm = ttpm_features, train_proc = 0.6)
            fwrite(mixed, "data_start.csv")
            fwrite(mixed, "mixed.csv")
          }} else {
            fwrite(dane, "data_start.csv")
            train = dane
            train$mix = "train"
            fwrite(train, "mixed_train.csv")
            test = dane
            test$mix = "test"
            fwrite(test, "mixed_test.csv")
            valid = dane
            valid$mix = "valid"
            fwrite(valid, "mixed_valid.csv")
            mixed = rbind(train,test,valid)
            fwrite(mixed, "mixed.csv")
          }
        })
        progress$inc(1/n_steps)
        
        try({
        OmicSelector_log("Checking if data can be loaded. Balancing dataset... (this may take a while)");
        dane = OmicSelector_load_datamix(use_smote_not_rose = T)  # load mixed_*.csv files
        train = dane[[1]]
        test = dane[[2]]
        valid = dane[[3]]
        train_smoted = dane[[4]]
        trainx = dane[[5]]
        trainx_smoted = dane[[6]]  # get the objects from list to make the code more readable.
        OmicSelector_log("\n✓ All datasets can be loaded. SMOTE-based balanced dataset can be loaded. ");
        })
        progress$inc(1/n_steps)
        
        try({
        OmicSelector_log("Performing DE analysis...");
        type = readLines("var_type.txt", warn = F)
        x_mix = dplyr::select(mixed, starts_with("hsa"))
        DE_mix = OmicSelector_differential_expression_ttest(x_mix, mixed$Class, mode = type)
        fwrite(DE_mix, "DE_mixed.csv")
        DE_train = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
        fwrite(DE_train, "DE_train.csv")
        OmicSelector_log("\n✓ DE was performed for whole dataset (mixed) and for training datasets. ");
        })
        progress$inc(1/n_steps)
        
        OmicSelector_log("Drawing exploratory plots...");
        try({ png("exploratory_pca.png", width = 1170, height = 658)
          a = OmicSelector_PCA(trainx, train$Class)
          a
          ggsave("exploratory_pca.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_pca_all_batch.png", width = 1170, height = 658)
          a = OmicSelector_PCA(dplyr::select(mixed, starts_with('hsa')), mixed$Batch)
          a
          ggsave("exploratory_pca_all_batch.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_pca_all_mix.png", width = 1170, height = 658)
          a = OmicSelector_PCA(dplyr::select(mixed, starts_with('hsa')), mixed$mix)
          a
          ggsave("exploratory_pca_all_mix.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_pca_all_class.png", width = 1170, height = 658)
          a = OmicSelector_PCA(dplyr::select(mixed, starts_with('hsa')), mixed$ClassOrginal)
          a
          ggsave("exploratory_pca_all_class.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_vulcano.png", width = 1170, height = 658)
          a = OmicSelector_vulcano_plot(DE_train$miR, DE = DE_train, only_label = DE_train$miR[DE_train$`p-value` < 0.05])
          a
          ggsave("exploratory_vulcano.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_vulcano2.png", width = 1170, height = 658)
          a = OmicSelector_vulcano_plot(DE_train$miR, DE = DE_train, only_label = DE_train$miR[1:10]) # top 10
          a
          ggsave("exploratory_vulcano2.png", a)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmap.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = trainx, rlab = data.frame(Class = train$Class), expression_name = "Expression")
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmapsig.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = trainx[,DE_train$miR[DE_train$`p-value` < 0.05]], rlab = data.frame(Class = train$Class), expression_name = "Expression")
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmapz.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = trainx, rlab = data.frame(Class = train$Class), expression_name = "Expression", zscore = T)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmapzsig.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = trainx[,DE_train$miR[DE_train$`p-value` < 0.05]], rlab = data.frame(Class = train$Class), expression_name = "Expression", zscore = T)
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmapmix.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix), expression_name = "Expression")
          suppressMessages(graphics.off()) })
        try({ png("exploratory_heatmapmixz.png", width = 1170, height = 658)
          OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix), expression_name = "Expression", zscore = T)
          suppressMessages(graphics.off()) })
        if(batch){
          try({ png("exploratory_pcab.png", width = 1170, height = 658)
            meta = paste0(mixed$Class, " - ", mixed$Batch)    
            OmicSelector_PCA(dplyr::select(mixed, starts_with("hsa")), meta)
            suppressMessages(graphics.off())  })
          try({ png("exploratory_heatmapmixb.png", width = 1170, height = 658)
            OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix, Batch = mixed$Batch), expression_name = "Expression")
            suppressMessages(graphics.off())  })
          try({ png("exploratory_heatmapmixzb.png", width = 1170, height = 658)
            OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix, Batch = mixed$Batch), expression_name = "Expression", zscore = T)
            suppressMessages(graphics.off())  })
        }
        OmicSelector_log("\n✓ PCA, vulcano plot and heatmaps were generated. ");
        progress$inc(1/n_steps)
        
        try({
        if(!file.exists("data_start.csv")) { fwrite(dane, "data_start.csv") }
        # Out: czy missing,czy batch
        writeLines(as.character(batch), "var_batch.txt", sep="")
        writeLines(as.character(missing), "var_missing.txt", sep="")
        writeLines(task_status, "var_initcheck.txt", sep="")
        OmicSelector_log("\n✓ The files are ready for further analysis. ");
        })
        
        OmicSelector_log("[OmicSelector: TASK COMPLETED]")
        
        output$task_log = shiny::renderPrint(cat(paste0(readLines("temp.log", warn =F), collapse = "\n")))
        output$status_log = shiny::renderPrint(cat(paste0(task_status, collapse = "\n")))
        
        wykresy = list.files(".", pattern = "*.png")
        updateSelectInput(session, "select_plot", choices = as.list(wykresy))
        
        dfs = list.files(".", pattern = "*.csv")
        updateSelectInput(session, "select_df", choices = as.list(dfs))
        
        
        
        progress$inc(1/n_steps)
        
        
        
        
        
        
        # END:
        shinyjs::toggle(id= "panelA")
        shinyjs::toggle(id= "panelB")
    })
  
  observeEvent(input$button_process, { 
    if(readLines("var_initcheck.txt", warn = F) != "OK") {
      showNotification("Cannot proceed due to file validation errors. Please correct the file and start over.", type = "error")
    } else {
      showNotification(paste0("Please wait... Creating analysis framework for ID: ", input$analysis_id))
      file.rename("temp.log","initial_check.txt")
      dir.create(paste0("/OmicSelector/",input$analysis_id))
      file.copy(list.files("."), paste0("/OmicSelector/",input$analysis_id), recursive = T)
      updateQueryString(paste0("/analysis.php?id=",input$analysis_id), mode = c("replace","push"),
                    session = getDefaultReactiveDomain())
      shinyjs::runjs(paste0('window.location.href = "/analysis.php?id=', input$analysis_id ,'";'))
    }
  })
  
  
}

shinyApp(ui=ui,server=server)