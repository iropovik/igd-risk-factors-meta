#+ setup, include=FALSE
library(readxl)
library(tidyverse)
library(irr)

rm(list = ls())
firstCoder <- read_delim("rater1.csv", delim = ";")
secondCoder <- read_delim("rater2.csv", delim = ";")

ratersData <- list(firstCoder, secondCoder)

variableNames <- c("focal", "design", "predictedDirection", "nTotal", "randomization",
                   "nMale", "nFemale", "sampleRestricted", "meanAge",
                   "Game_complexity","possibleCOI" ,"Risk of bias coder", "published",
                   "mExp", "sdExp",
                   "mCtrl", "sdCtrl", "nExp", "nCtrl", "controlsActive",
                   "Fstat", "tstat", "rstat",
                   "df1", "df2", "interventionDuration(m)", "posttestDelay(day)",
                   "Att_type_exp_imp", "inLabAdministration", "controlledGameplay",
                   "persuasive_mechanic", "Commercial")

ratersData <- lapply(ratersData, function(x){
  x %>% select(paperID, studyID, effectID, starts_with(variableNames)) %>%
    mutate(paperID = as.numeric(paperID),
           studyID = as.numeric(studyID),
           effectID = as.numeric(effectID))})

datCut <- as_tibble(left_join(ratersData[[1]], ratersData[[2]], by = c("paperID", "studyID", "effectID"), suffix = c("_first", "_second"), keep = TRUE))
datCut <- datCut[,order(colnames(datCut),decreasing=TRUE)] %>% select(paperID_second, paperID_first, studyID_second, studyID_first, effectID_second, effectID_first, everything())
#View(datCut)

datCut %>% select(starts_with(variableNames, ignore.case = F)) %>% apply(., 2, table)

kappas <- list()
for(i in 1:length(variableNames)){
  kappas[[i]] <- datCut %>% select(starts_with(variableNames[i], ignore.case = F)) %>% kappa2(.)
}
names(kappas) <- variableNames

agreement <- list()
for(i in 1:length(variableNames)){
  agreement[[i]] <- datCut %>% select(starts_with(variableNames[i])) %>% agree(.)
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Cohen's kappa
kappas
#'# Inter-rater agreement
agreement

