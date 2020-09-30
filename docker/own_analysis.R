library(OmicSelector)
dane = OmicSelector_load_datamix(use_smote_not_rose = T)  # load mixed_*.csv files
train = dane[[1]]
test = dane[[2]]
valid = dane[[3]]
train_smoted = dane[[4]]
trainx = dane[[5]]
trainx_smoted = dane[[6]]  # get the objects from list to make the code more readable.