options(warn=-1)
suppressMessages(suppressMessages(library(data.table)))
suppressMessages(suppressMessages(library(dplyr)))
writeLines("ERROR", "var_initcheck.txt", sep="")
if(file.exists("data.csv")) { dane = fread("data.csv") } else {
    if(file.exists("data.xlsx")) { dane = xlsx::read.xlsx("data.xlsx",sheetIndex=1) }
}
colnames(dane) = make.names(colnames(dane), unique = T)
dane[dane == ""] = NA

error = FALSE
if("Class" %in% colnames(dane)) { cat("✓ The file contains Class variable. ") } else 
{ writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); stop("☒ The file contains DOES NOT Class variable. ") }

dane$Class = factor(dane$Class, levels = c("Control","Cancer"))
if(table(dane$Class)[1] == 0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); stop("☒ There are no control cases.") }
if(table(dane$Class)[2] == 0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); stop("☒ There are no cancer cases.") }
cat(paste0("\n✓ The data contains ", table(dane$Class)[1], " `Control` cases and ", table(dane$Class)[2], " `Cancer` cases."))

temp = dplyr::select(dane, starts_with("hsa"))
if(ncol(temp)==0) { writeLines(as.character("FAIL"), "var_initcheck.txt", sep=""); stop("☒ The data does not contain any features (e.g. miRNAs) for feautre selection. Remember that feature names should start from hsa...") }
cat(paste0("\n✓ The data contains ", ncol(temp), " features (e.g. miRNAs) for selection."))

czy_numeryczne = sapply(temp, is.numeric)
if(sum(czy_numeryczne) != ncol(temp)) { writeLines(as.character("FAIL", "var_initcheck.txt", sep = "")); stop("☒ Some of the features are not numeric. Please remove them. Not numeric: ", paste0(colnames(temp)[czy_numeryczne == F], collapse = ", "))}
cat(paste0("\n✓ All features are numeric."))

missing = FALSE
czy_brakna = sapply(temp, is.na)
if(sum(colSums(czy_brakna)) != 0) 
{ stop("☒ Some of the features contain missing data.\nWith missing values:\n  -", paste0(colnames(temp)[colSums(czy_brakna) > 0], collapse = "\n  - "))
missing = T } else 
{ cat(paste0("\n✓ There are no missing data in features.")) }

batch = F
if("Batch" %in% colnames(dane)) { cat("\n✓ The file contains `Batch` variable that can be used for batch-effect correction. You can apply combat-based batch correction via OmicSelector, but not via GUI. Please see the package manual. Please also check if the following contingency table is correct:\n ")
print(table(dane$Class, dane$Batch))
batch = T } else { cat("\n✓ The file does not contain `Batch` variable that can be used for batch-effect correction.") }

x = dplyr::select(dane, starts_with("hsa"))
like_counts = sapply(x, function(x2) (sum(unlist(na.omit(x2))%%1 == 0) + sum(unlist(na.omit(x2)) >= 0))/(2*length(unlist(na.omit(x2)))) == 1)
positive = F
if(sum(like_counts)/ncol(x)) { cat("\n✓ Feature values are positive integers. The file could represent read counts. Please note that the feature selection pipeline requires normalized data. You can use OmicSelector to convert counts to tpm, but not via GUI. Please refer to package manual."); positive = T; } else {
    cat("\n✓ Feature values are not positive integers. This is ok if your input data is normalized (e.g. deltaCt values or tpm-normalized counts). Features not looking like counts: ");
    cat(paste0(colnames(x)[which(like_counts == FALSE)], collapse = ", "))
}
writeLines(as.character(positive), "var_seemslikecounts.txt", sep="")

suppressMessages(library(OmicSelector))
split = readLines("var_split.txt", warn = F)
if(split == "yes") {
if("mix" %in% colnames(dane)) {
    if(sum(which(dane$mix == 'train')) > 0) { 
    cat(paste0("\n✓ Samples in training set: ", sum(dane$mix == 'train')));
    } else {
        cat(paste0("\n✓ Samples in training set: ", sum(which(dane$mix == 'train'))));
        stop("There are no samples in training set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
    }
    if(sum(which(dane$mix == 'test')) > 0) { 
    cat(paste0("\n✓ Samples in test set: ", sum(dane$mix == 'test')));
    } else {
        cat(paste0("\n✓ Samples in test set: ", sum(dane$mix == 'test')));
        stop("There are no samples in test set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
    }
    if(sum(which(dane$mix == 'valid')) > 0) { 
    cat(paste0("\n✓ Samples in validation set: ", sum(dane$mix == 'valid')));
    } else {
        cat(paste0("\n✓ Samples in validation set: ", sum(dane$mix == 'valid')));
        stop("There are no samples in validation set! Please fix the 'mix' variable. This variable should contain the assignment to 'train', 'test' and 'valid' groups.")
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
} else {
    cat("\n✓ The data is not splitted, i.e. doesn't have 'train', 'test' and 'valid' in 'mix' variable. We will perform data splitting (60% train, 20% test, 20% valid). ");
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

dane = OmicSelector_load_datamix(use_smote_not_rose = T)  # load mixed_*.csv files
train = dane[[1]]
test = dane[[2]]
valid = dane[[3]]
train_smoted = dane[[4]]
trainx = dane[[5]]
trainx_smoted = dane[[6]]  # get the objects from list to make the code more readable.
cat("\n✓ All datasets can be loaded. SMOTE-based balancing can be performed. ");

type = readLines("var_type.txt", warn = F)
x_mix = dplyr::select(mixed, starts_with("hsa"))
DE_mix = OmicSelector_differential_expression_ttest(x_mix, mixed$Class, mode = type)
fwrite(DE_mix, "DE_mixed.csv")
DE_train = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
fwrite(DE_train, "DE_train.csv")
cat("\n✓ DE was performed for whole dataset (mixed) and for training datasets. ");

png("exploratory_pca.png", width = 1170, height = 658)
OmicSelector_PCA(trainx, train$Class)
suppressMessages(graphics.off())
png("exploratory_vulcano.png", width = 1170, height = 658)
OmicSelector_vulcano_plot(DE_train$miR, DE = DE_train, only_label = DE_train$miR[DE_train$`p-value` < 0.05])
suppressMessages(graphics.off())
png("exploratory_heatmap.png", width = 1170, height = 658)
OmicSelector_heatmap(x = trainx, rlab = data.frame(Class = train$Class), expression_name = "Expression")
suppressMessages(graphics.off())
png("exploratory_heatmapsig.png", width = 1170, height = 658)
OmicSelector_heatmap(x = trainx[,DE_train$miR[DE_train$`p-value` < 0.05]], rlab = data.frame(Class = train$Class), expression_name = "Expression")
suppressMessages(graphics.off())
png("exploratory_heatmapz.png", width = 1170, height = 658)
OmicSelector_heatmap(x = trainx, rlab = data.frame(Class = train$Class), expression_name = "Expression", zscore = T)
suppressMessages(graphics.off())
png("exploratory_heatmapzsig.png", width = 1170, height = 658)
OmicSelector_heatmap(x = trainx[,DE_train$miR[DE_train$`p-value` < 0.05]], rlab = data.frame(Class = train$Class), expression_name = "Expression", zscore = T)
suppressMessages(graphics.off())
png("exploratory_heatmapmix.png", width = 1170, height = 658)
OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix), expression_name = "Expression")
suppressMessages(graphics.off())
png("exploratory_heatmapmixz.png", width = 1170, height = 658)
OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix), expression_name = "Expression", zscore = T)
suppressMessages(graphics.off())
if(batch){
png("exploratory_pcab.png", width = 1170, height = 658)
meta = paste0(mixed$Class, " - ", mixed$Batch)    
OmicSelector_PCA(dplyr::select(mixed, starts_with("hsa")), meta)
suppressMessages(graphics.off())
png("exploratory_heatmapmixb.png", width = 1170, height = 658)
OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix, Batch = mixed$Batch), expression_name = "Expression")
suppressMessages(graphics.off())
png("exploratory_heatmapmixzb.png", width = 1170, height = 658)
OmicSelector_heatmap(x = dplyr::select(mixed, DE_train$miR[DE_train$`p-value` < 0.05]), rlab = data.frame(Class = mixed$Class, Mix = mixed$mix, Batch = mixed$Batch), expression_name = "Expression", zscore = T)
suppressMessages(graphics.off())
}
cat("\n✓ PCA, vulcano plot and heatmaps were generated. ");

if(!file.exists("data_start.csv")) { fwrite(dane, "data_start.csv") }
# Out: czy missing,czy batch
writeLines(as.character(batch), "var_batch.txt", sep="")
writeLines(as.character(missing), "var_missing.txt", sep="")
writeLines("OK", "var_initcheck.txt", sep="")
cat("\n✓ The files are ready for further analysis. ");