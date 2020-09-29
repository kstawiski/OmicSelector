## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, message=FALSE, warning=FALSE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=12, fig.height=8)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=150),tidy=TRUE)
options(rgl.useNULL = TRUE)
options(warn=-1)
suppressMessages(library(dplyr))
set.seed(1)
options(knitr.table.format = "html")
library(OmicSelector)

## -----------------------------------------------------------------------------
readLines("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/vignettes/setup.R") %>% paste0(collapse="\n") %>% cat

## ----setup--------------------------------------------------------------------
library(OmicSelector)

## ---- eval = F----------------------------------------------------------------
#  OmicSelector_download_tissue_miRNA_data_from_TCGA()
#  OmicSelector_process_tissue_miRNA_TCGA(remove_miRNAs_with_null_var = T)

## -----------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(knitr)))
data("orginal_TCGA_data")
OmicSelector_table(table(orginal_TCGA_data$primary_site, orginal_TCGA_data$sample_type))

## -----------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(dplyr)))

cancer_cases = filter(orginal_TCGA_data, primary_site == "Pancreas" & sample_type == "PrimaryTumor")
control_cases = filter(orginal_TCGA_data, sample_type == "SolidTissueNormal")

## -----------------------------------------------------------------------------
cancer_cases$Class = "Cancer"
control_cases$Class = "Control"

dataset = rbind(cancer_cases, control_cases)

OmicSelector_table(table(dataset$Class), col.names = c("Class","Number of cases"))

## -----------------------------------------------------------------------------
boxplot(dataset$age_at_diagnosis ~ dataset$Class)
t.test(dataset$age_at_diagnosis ~ dataset$Class)
OmicSelector_table(table(dataset$gender.x, dataset$Class))
chisq.test(dataset$gender.x, dataset$Class)

## -----------------------------------------------------------------------------
old_dataset = dataset # backup
dataset = dataset[grepl("Adenocarcinomas", dataset$disease_type),]
match_by = c("age_at_diagnosis","gender.x")
tempdane = dplyr::select(dataset, match_by)
tempdane$Class = ifelse(dataset$Class == "Cancer", TRUE, FALSE)
suppressMessages(library(mice))
suppressMessages(library(MatchIt))
temp1 = mice(tempdane, m=1)
temp2 = temp1$data
temp3 = mice::complete(temp1)
temp3 = temp3[complete.cases(temp3),]
tempform = OmicSelector_create_formula(match_by)
mod_match <- matchit(tempform, data = temp3)
newdata = match.data(mod_match)
dataset = dataset[as.numeric(rownames(newdata)),]


## -----------------------------------------------------------------------------
boxplot(dataset$age_at_diagnosis ~ dataset$Class)
t.test(dataset$age_at_diagnosis ~ dataset$Class)
OmicSelector_table(table(dataset$gender.x, dataset$Class))
chisq.test(dataset$gender.x, dataset$Class)
fwrite(dataset, "balanced_dataset.csv.gz")
OmicSelector_tutorial_balanced_dataset = dataset # can be used by data("OmicSelector_tutorial_balanced_dataset")

## -----------------------------------------------------------------------------
dataset = OmicSelector_correct_miRNA_names(dataset) # Correct miRNA names based on the aliases. Useful when analyzing old datasets - to keep the results coherent with current knowledge.
danex = dplyr::select(dataset, starts_with("hsa")) # Create data.frame or matrix with miRNA counts with miRNAs in columns and cases in rows.
metadane = dplyr::select(dataset, -starts_with("hsa")) # Metadata with 'Class' variables.
OmicSelector_table(table(metadane$Class)) # Let's be sure that 'Class' variable is correct and contains only 'Cancer' and 'Control' cases.
ttpm = OmicSelector_counts_to_log10tpm(danex, metadane, ids = metadane$sample,
                                 filtr = T, filtr_minimalcounts = 100, filtr_howmany = 1/3) # We will leave only the miRNAs which apeared with at least 100 counts in 1/3 of cases.

## -----------------------------------------------------------------------------
mixed = OmicSelector_prepare_split(metadane = metadane, ttpm = ttpm, train_proc = 0.6)
OmicSelector_tutorial_balanced_mixed = mixed # can be used by data("OmicSelector_tutorial_balanced_mixed")

## -----------------------------------------------------------------------------
mixed = fread("mixed.csv")
OmicSelector_table(table(mixed$Class, mixed$mix))
OmicSelector_table(cbind(mixed[1:10,c(100:105)], Class = mixed[1:10,"Class"]))

## ----warning=FALSE------------------------------------------------------------
dane = OmicSelector_load_datamix(use_smote_not_rose = T) # load mixed_*.csv files
train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]] # get the objects from list to make the code more readable.

## -----------------------------------------------------------------------------
pca = OmicSelector_PCA(trainx, train$Class)
pca

## -----------------------------------------------------------------------------
if(is.null(sessionInfo()$loadedOnly$IRdisplay)) { # if not in the Jupyter, if you run OmicSelector_PCA_3D in learning/editing Jupyter enviorment it may cause: *** caught segfault *** address 0x1, cause 'memory not mapped'
pca3d = OmicSelector_PCA_3D(trainx, train$Class)
pca3d }

## -----------------------------------------------------------------------------
de = OmicSelector_differential_expression_ttest(trainx, train$Class)
sig_de = de %>% dplyr::filter(`p-value BH` <= 0.05) %>% dplyr::arrange(`p-value BH`) # leave only significant after Benjamini-Hochberg procedure and sort by ascending p-value
OmicSelector_table(sig_de) 

## -----------------------------------------------------------------------------
OmicSelector_heatmap(x = dplyr::select(trainx, sig_de$miR),
           rlab = data.frame(Class = train$Class),
           zscore = F, margins = c(10,10))

## -----------------------------------------------------------------------------
OmicSelector_heatmap(x = dplyr::select(trainx, sig_de$miR),
           rlab = data.frame(Class = train$Class),
           zscore = T, margins = c(10,10))

## -----------------------------------------------------------------------------
OmicSelector_vulcano_plot(selected_miRNAs = de$miR, DE = de, only_label = sig_de$miR[1:10])

## -----------------------------------------------------------------------------
de_test = OmicSelector_differential_expression_ttest(dplyr::select(test, starts_with("hsa")), test$Class)
de_valid = OmicSelector_differential_expression_ttest(dplyr::select(valid, starts_with("hsa")), valid$Class)
OmicSelector_correlation_plot(de$log2FC, de_test$log2FC, "log2FC on training set", "log2FC on test set", "", yx = T)
OmicSelector_correlation_plot(de$log2FC, de_valid$log2FC, "log2FC on training set", "log2FC on validation set", "", yx = T)
OmicSelector_correlation_plot(de_test$log2FC, de_valid$log2FC, "log2FC on test set", "log2FC on validation set", "", yx = T)

## -----------------------------------------------------------------------------
library(OmicSelector)
selected_features = OmicSelector_OmicSelector(wd = getwd(), m = c(1:3,51), max_iterations = 1, stamp = "tutorial") # For the sake of this tutorial and vignette building we will use only few fastest methods. The m parameter defines what methods will be tested. See more details below.

## -----------------------------------------------------------------------------
readLines("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/vignettes/Tutorial_OmicSelector.R") %>% paste0(collapse="\n") %>% cat

## -----------------------------------------------------------------------------
shiny::includeHTML("methods.html")

## -----------------------------------------------------------------------------
selected_sets_of_miRNAs = OmicSelector_merge_formulas(max_miRNAs = 11) # we filter out sets with more than 11 miRNAs.
selected_sets_of_miRNAs_with_own = OmicSelector_merge_formulas(max_miRNAs = 11, 
                                                     add = list("my_own_signature" = c("hsa.miR.192.5p","hsa.let.7g.5p","hsa.let.7a.5p","hsa.let.7d.5p","hsa.miR.194.5p", "hsa.miR.98.5p", "hsa.let.7f.5p", "hsa.miR.26b.5p"))) # you can also add your own signature (for example selected from literature)

## -----------------------------------------------------------------------------
all_sets = readRDS("featureselection_formulas_all.RDS")
length(all_sets) # How many feature selection methods completed in time?
final_sets = readRDS("featureselection_formulas_final.RDS")
length(final_sets) # How many feature selection methods completed in time and fulfilled max_miRNA criteria? (remember about fcsig and cfs_sig)
featureselection_formulas_final = fread("featureselection_formulas_final.csv")
OmicSelector_table(featureselection_formulas_final) # show information about selected formulas

## -----------------------------------------------------------------------------
hist(featureselection_formulas_final$ile_miRNA[-which(featureselection_formulas_final$ile_miRNA == 0)], 
     breaks = ncol(train),
     main = "Number of selected microRNAs distribution",
     xlab = "Number of selected microRNAs"
     ) # Histogram showing how many miRNAs were selected in final set.
psych::describe(featureselection_formulas_final$ile_miRNA[-which(featureselection_formulas_final$ile_miRNA == 0)]) # Descriptive statistics of how many features where selected in the final set.

## ---- eval=F------------------------------------------------------------------
#  readLines("Tutorial_benchmark.R") %>% paste0(collapse="\n") %>% cat

## -----------------------------------------------------------------------------
library(OmicSelector)
OmicSelector_tutorial_balanced_benchmark = OmicSelector_benchmark(search_iters = 5, # 5 random hyperparameter sets will be checked; 5 is set here for speed purposes.. for real projects use more, like 5000...
            algorithms = c("ctree"), # just add ctree, note that logistic regression (glm) is always included
            output_file = paste0("benchmark.csv")) # the main output
# exemplary benchmark data can be loaded using data('OmicSelector_tutorial_balanced_benchmark')

## -----------------------------------------------------------------------------
OmicSelector_table(fread("benchmark.csv"))

## -----------------------------------------------------------------------------
metody = OmicSelector_get_benchmark_methods("benchmark.csv") # gets the methods used in benchmark
par(mfrow = c(2,2))
for(i in 1:length(metody)){
    temp = OmicSelector_get_benchmark("benchmark.csv") # loads benchmark
    temp2 = dplyr::select(temp, starts_with(paste0(metody[i],"_")))
    boxplot(temp[,paste0(metody[i],"_train_Accuracy")], temp[,paste0(metody[i],"_test_Accuracy")], temp[,paste0(metody[i],"_valid_Accuracy")],
            main = paste0("Method: ", metody[i]), names = c("Training","Testing","Validation"), ylab = "Accuracy", ylim = c(0.5,1))
    tempids = c(match(paste0(metody[i],"_train_Accuracy"), colnames(temp)), match(paste0(metody[i],"_test_Accuracy"), colnames(temp)), match(paste0(metody[i],"_valid_Accuracy"), colnames(temp)))
  }
par(mfrow = c(1,1))

## -----------------------------------------------------------------------------
acc1 = OmicSelector_best_signiture_proposals(benchmark_csv = "benchmark.csv", without_train = F) # generates the benchmark sorted by metaindex
best_signatures = acc1[1:3,] # get top 3 methods
OmicSelector_table(best_signatures[,c("metaindex","method","miRy")])

## -----------------------------------------------------------------------------
acc1 = OmicSelector_best_signiture_proposals(benchmark_csv = "benchmark.csv", without_train = T) # generates the benchmark sorted by metaindex
best_signatures = acc1[1:3,] # get top 3 methods
OmicSelector_table(best_signatures[,c("metaindex","method","miRy")])

## -----------------------------------------------------------------------------
acc = OmicSelector_best_signiture_proposals_meta11(benchmark_csv = "benchmark.csv") # generates the benchmark sorted by metaindex
best_signatures = acc[1:3,] # get top 3 methods
OmicSelector_table(best_signatures[,c("metaindex","method","miRy")])

## ----fig.height=9, fig.width=16-----------------------------------------------
for(i in 1:length(metody))
  {
suppressMessages(library(PairedData))
suppressMessages(library(profileR))
pd = paired(as.numeric(acc[1:3,paste0(metody[i],"_train_Accuracy")]),as.numeric(acc[1:3,paste0(metody[i],"_test_Accuracy")]))
colnames(pd) = c("Train Accuracy","Test Accuracy")
plot2 = OmicSelector_profileplot(pd, Method.id = acc$method[1:3], standardize = F)
pd = paired(as.numeric(acc[1:3,paste0(metody[i],"_train_Accuracy")]),as.numeric(acc[1:3,paste0(metody[i],"_valid_Accuracy")]))
colnames(pd) = c("Train Accuracy","Valid Accuracy")
plot3 = OmicSelector_profileplot(pd, Method.id = acc$method[1:3], standardize = F)
pd = paired(as.numeric(acc[1:3,paste0(metody[i],"_test_Accuracy")]),as.numeric(acc[1:3,paste0(metody[i],"_valid_Accuracy")]))
colnames(pd) = c("Test Accuracy","Valid Accuracy")
plot4 = OmicSelector_profileplot(pd, Method.id = acc$method[1:3], standardize = F)



require(gridExtra)
grid.arrange(arrangeGrob(plot2, plot3, ncol=2, nrow = 1, top=metody[i]))
grid.arrange(arrangeGrob(plot4, ncol=1, nrow = 1, top=metody[i]))
}

## -----------------------------------------------------------------------------
acc2 = acc[1:6,] # get top 6 methods
accmelt = melt(acc2, id.vars = "method") %>% filter(variable != "metaindex") %>% filter(variable != "miRy")
accmelt = cbind(accmelt, strsplit2(accmelt$variable, "_"))
acctest = accmelt$value[accmelt$`2` == "test"]
accvalid = accmelt$value[accmelt$`2` == "valid"]
accmeth = accmelt$method[accmelt$`2` == "test"]
unique(accmeth)
plot5 = ggplot(, aes(x = as.numeric(acctest), y = as.numeric(accvalid), shape = accmeth)) +
  geom_point() + scale_x_continuous(name="Accuracy on test set", limits=c(0.5, 1)) +
  scale_y_continuous(name="Accuracy on validation set", limits=c(0.5, 1)) +
  theme_bw()
grid.arrange(arrangeGrob(plot5, ncol=1, nrow = 1))

## -----------------------------------------------------------------------------
OmicSelector_table(best_signatures[1:3,2:4])

## -----------------------------------------------------------------------------
selected_miRNAs = OmicSelector_get_features_from_benchmark(benchmark_csv = "benchmark.csv", best_signatures$method[1]) # for the best performing signiture
gsub("\\.", "-", selected_miRNAs) # R doesn't like hyphens, but we can introduce them easly

## -----------------------------------------------------------------------------
best_de = OmicSelector_best_signiture_de(selected_miRNAs)
OmicSelector_table(best_de)

## ----fig.height=9, fig.width=16-----------------------------------------------
for(i in 1:3){
  cat(paste0("\n\n## ", acc$method[i],"\n\n"))
  par(mfrow = c(1,2))
  acc = OmicSelector_best_signiture_proposals_meta11("benchmark.csv")
  metody = OmicSelector_get_benchmark_methods("benchmark.csv")
  ktory_set = match(acc$method[i], OmicSelector_get_benchmark("benchmark.csv")$method)
  #do_ktorej_kolumny = which(colnames(acc) == "metaindex")
  #barplot(as.numeric(acc[i,1:do_ktorej_kolumny]))
  for(ii in 1:length(metody)) {
    
    temp = OmicSelector_get_benchmark("benchmark.csv") %>% 
      dplyr::select(starts_with(paste0(metody[ii],"_t")),starts_with(paste0(metody[ii],"_v")))
    
    ROCtext = paste0("Training AUC ROC: ", round(temp[ktory_set,1],2), " (95%CI: ", round(temp[ktory_set,2],2), "-", round(temp[ktory_set,3],2), ")")
    
    temp = temp[,-c(1:3)]
    temp2 = as.numeric(temp[ktory_set,])
    temp3 = matrix(temp2, nrow = 3, byrow = T)
    colnames(temp3) = c("Accuracy","Sensitivity","Specificity")
    rownames(temp3) = c("Training","Testing","Validation")
    temp3 = t(temp3)
    
    plot1 = barplot(temp3, beside=T, ylim = c(0,1), xlab = paste0(ROCtext,"\nBlack - accuracy, blue - sensitivity, green - specificity"), width = 0.85, col=c("black", "blue", "green"), legend = F,  args.legend = list(x="topright", bty = "n", inset=c(0, -0.25)), cex.lab=0.7, main = paste0(acc$method[i], " - ", metody[ii]), font.lab=2)
    ## Add text at top of bars
    text(x = plot1, y = as.numeric(temp3), label = paste0(round(as.numeric(temp[ktory_set,])*100,1),"%"), pos = 3, cex = 0.6, col = "red")
  }
  par(mfrow = c(1,1))

}

## -----------------------------------------------------------------------------
overlap = OmicSelector_signiture_overlap(acc$method[1:3], "benchmark.csv")

## -----------------------------------------------------------------------------
attr(overlap,"intersections")

## ----warning=FALSE------------------------------------------------------------
OmicSelector_vulcano_plot(selected_miRNAs = de$miR, DE = de, only_label = selected_miRNAs)

## -----------------------------------------------------------------------------
OmicSelector_heatmap(x = dplyr::select(mixed, gsub("\\.", "-", selected_miRNAs)),
           rlab = data.frame(Class = mixed$Class, Mix = mixed$mix),
           zscore = F, margins = c(10,10))

## -----------------------------------------------------------------------------
OmicSelector_heatmap(x = dplyr::select(mixed, gsub("\\.", "-", selected_miRNAs)),
           rlab = data.frame(Class = mixed$Class, Mix = mixed$mix),
           zscore = T, margins = c(10,10))

## -----------------------------------------------------------------------------
cat(paste0(gsub("\\.", "-", selected_miRNAs), collapse = ", "))

## -----------------------------------------------------------------------------
session_info()

## -----------------------------------------------------------------------------
packageDescription("OmicSelector")

## ---- eval = FALSE------------------------------------------------------------
#  render("Tutorial.Rmd", output_file = "Tutorial.html", output_dir = "../inst/doc/")

## -----------------------------------------------------------------------------
OmicSelector_table(as.data.frame(installed.packages()))

## -----------------------------------------------------------------------------
unlink("temp", recursive=TRUE)
unlink("models", recursive=TRUE)

