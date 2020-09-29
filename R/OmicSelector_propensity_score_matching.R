#' OmicSelector_propensity_score_matching
#'
#' Propensity score matching.
#'
#' @param dataset Original dataset.
#' @param match_by Vector describing by which variables should the dataset be matched.
#' @param method Passed to `matchit()`.
#' @param distance Passed to `matchit()`.
#'
#' @export
#'
OmicSelector_propensity_score_matching = function(dataset, match_by = c("age_at_diagnosis","gender.x"), method = "nearest", distance = "logit"){
  suppressMessages(library(MatchIt))
  suppressMessages(library(mice))

  #tempdane = dataset
  tempdane = dplyr::select(dataset, match_by)
  tempdane$Class = ifelse(dataset$Class == "Cancer", TRUE, FALSE)
  suppressMessages(library(mice))
  temp1 = mice(tempdane, m=1)
  temp2 = temp1$data
  temp3 = mice::complete(temp1)


  temp3 = temp3[complete.cases(temp3),]

  tempform = OmicSelector_create_formula(match_by)
  mod_match <- MatchIt::matchit(tempform,
                       method = method, distance = distance,
                       data = temp3)

  newdata = match.data(mod_match)
  newdatafull = dataset[as.numeric(rownames(newdata)),]
  return(newdatafull)
}
