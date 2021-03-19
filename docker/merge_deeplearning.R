library(OmicSelector)
lista_plikow = list.files(".", pattern = "^deeplearning.*.csv$")
library(plyr)
wyniki = data.frame()
for(i in 1:length(lista_plikow)) { temp = data.table::fread(lista_plikow[i]); wyniki = rbind.fill(wyniki, temp); }


temp = dplyr::select(wyniki, training_Accuracy, test_Accuracy, valid_Accuracy)
temp = t(temp)
wyniki$metaindex = psych::harmonic.mean(temp)

# for(i in 1:nrow(wyniki))
# {
#   wyniki[i,"metaindex"] = psych::harmonic.mean(c(wyniki[i,"DORtrain"], wyniki[i,"DORtest"], wyniki[i,"DORvalid"]))
# }

#hist(wyniki$metaindex)
#summary(wyniki$metaindex)

wyniki$metaindex2 = (wyniki$training_Accuracy + wyniki$test_Accuracy + wyniki$valid_Accuracy) / 3

data.table::fwrite(wyniki, "merged_deeplearning.csv")
library(dplyr)
wynikitop = wyniki %>% arrange(desc(metaindex)) %>% filter(worth_saving == TRUE)
if(nrow(wynikitop)<1000) { max_top = nrow(wynikitop) } else { max_top = 1000}
wynikitop = wynikitop[1:max_top,]
data.table::fwrite(wynikitop, "merged_deeplearning_top.csv")
data.table::fwrite(as.data.frame(wynikitop$name), "merged_deeplearning_names.csv")

