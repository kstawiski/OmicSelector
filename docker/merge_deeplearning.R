library(OmicSelector)
lista_plikow = list.files(".", pattern = "^deeplearning.*.csv$")
library(plyr)
wyniki = data.frame()
for(i in 1:length(lista_plikow)) { temp = data.table::fread(lista_plikow[i]); wyniki = rbind.fill(wyniki, temp); }
data.table::fwrite(wyniki, "merged_deeplearning.csv")