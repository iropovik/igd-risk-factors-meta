#' ---
#' title: "Risk factors for Gaming Disorder: A meta-analysis"
#' author: "Ivan Ropovik, Marcel Martončik, Gabriel Baník, Lenka Vargová, Peter Babinčák, & Matúš Adamkovič"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

#+ setup, include = FALSE
# NOTE: Please note that to run the script, you need the development versions of metafor and dmetar packages from github.

rm(list = ls())
startTime <- Sys.time()

# Settings ----------------------------------------------------------------
# What is the minimum number of effects to be synthesized
kThreshold <- 10

# Should bias-correction methods be applied to meta-analytic models?
biasOn <- FALSE

# Should the meta-analytic models exclude outlying, excessively influential effects?
outlierSensitivity <- FALSE

# Should the meta-analytic models exclude effects based on selection inference approaches?
selInfSensitivity <- FALSE

# Controls for the multiple-parameter selection models 
# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 3 # Set to 5 just to make code checking/running fast. For the final analysis, it should be set to 1000
nIterationsPcurve <- 2
nIterationVWsensitivity <- 2 # 200 Number of iterations for the Vevea & Woods (2005) step function model sensitivity analysis 
# Number of chains and iterations for Robust Bayesian model-averaging approach
runRobMA <- TRUE
robmaChains <- 2
robmaSamples <- 100

# Whether to apply a 4- or 3-parameter selection model. If fallback == TRUE, the procedure falls back to the 3-parameter selection model. 
# This should be selected when too few effects in the opposite side make the estimate unstable.
fallback <- TRUE

# Even when fallback == FALSE, the 4-parameter selection model still falls back to 3 parameters for the given iteration if,
# (1) it fails to converge or (2) the number of p-values in each of the step intervals gets smaller than minPvalues.
minPvalues <- 4

# Steps and delta parameters for Vevea & Woods selection models 
# Can be adjusted if a different selection process is assumed. 
# Please note that steps vector represents one-tailed p-values.
stepsDelta <- data.frame(
  steps =     c(.0025, .005, .0125, .025, .05, .10, .25, .50, 1),
  moderateSelection = c(1, 0.99, 0.97, 0.95, 0.80, 0.60, 0.50, 0.50, 0.50),
  severeSelection =   c(1, 0.99, 0.97, 0.95, 0.65, 0.40, 0.25, 0.25, 0.25),
  extremeSelection =  c(1, 0.98, 0.95, 0.90, 0.50, 0.20, 0.10, 0.10, 0.10))

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
# side <- "right"
# For this application, we derived the side argument empirically from the direction of the majority of effect sizes, establishing risk and protective factors
# Inside note: code edit for pre-defined side argument: puni_star, petPeese, powerEst

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# Controls for PET-PEESE
condEst <- FALSE

# Assumed constant sampling correlation
rho <- 0.5

# Correct for indirect selection bias. If FALSE, only direct uni/bivariate selection will be accounted for.
indirectSel <- FALSE

# Should the plots be displayed? If FALSE, plots will only be stored in respective list objects
displayPlots <- FALSE

# Sourcing and reading data -----------------------------------------------------------------
#+ include = FALSE
source("functions.R")
source("pcurvePlotOption.R")
source("esConversion.R")
if(outlierSensitivity == TRUE){source("maDiag.R")}
statcheck <- read.csv("statcheck.csv")
#dat <- dat[-c(129, 131, 132, 321),]

knitr::opts_chunk$set(echo = FALSE, warning = FALSE)

#'#### Appendix E: Fully analytic output
#'
#' **This is the supplementary analytic output for the paper Risk factors for Gaming Disorder: A meta-analysis.**
#' 
#' **It reports detailed results for all models reported in the paper. The analytic R script by which this html report was generated can be found on the project's OSF page at: [LINK].**
#' **Preprint of the paper is available here, Appendices A-D are available here: **
#' 
#' ------------------------------------
#' 
#' **Brief information about the methods used in the analysis:**
#' 
#' **RMA results with model-based SEs**
#' k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#' **RVE SEs with Satterthwaite small-sample correction**
#' Estimate based on a multilevel RE model with constant sampling correlation model (CHE - correlated hierarchical effects - working model) (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/). 
#' Interpretation of naive-meta-analysis should be based on these estimates.
#'
#' **Prediction interval**
#' Shows the expected range of true effects in similar studies.
#' As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between PI LB and PI UB.
#' Note that these are non-adjusted estimates. An unbiased newly conducted study will more likely fall in an interval centered around bias-adjusted estimate with a wider CI width.
#'
#' **Heterogeneity**
#' Tau can be interpreted as the total amount of heterogeneity in the true effects. 
#' I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates. Estimates calculated by two approaches are reported.
#' This is followed by separate estimates of between- and within-cluster heterogeneity and estimated intra-class correlation of underlying true effects.
#' 
#' **Proportion of significant results**
#' What proportion of effects were statistically at the alpha level of .05.
#' 
#' **ES-precision correlation**
#' Kendalls's correlation between the ES and precision.
#' 
#' **4/3PSM**
#' Applies a permutation-based, step-function 4-parameter selection model (one-tailed p-value steps = c(.025, .5, 1)). 
#' Falls back to 3-parameter selection model if at least one of the three p-value intervals contains less than 5 p-values.
#' For this meta-analysis, we applied 3-parameter selection model by default as there were only 11 independent effects in the opposite direction overall (6%), causing the estimates to be unstable across iterations.
#' pvalue = p-value testing H0 that the effect is zero. ciLB and ciUB are lower and upper bound of the CI. k = number of studies. steps = 3 means that the 4PSM was applied, 2 means that the 3PSM was applied.
#' We also ran two sensitivity analyses of the selection model, the Vevea & Woods (2005) step function model with a priori defined selection weights and the Robust Bayesian Meta-analysis model employing the model-averaging approach (Bartoš & Maier, 2020).
#' 
#' **PET-PEESE**
#' Estimated effect size of an infinitely precise study. Using 4/3PSM as the conditional estimator instead of PET (can be changed to PET). If the PET-PEESE estimate is in the opposite direction, the effect can be regarded nil. 
#' By default (can be changed to PET), the function employs a modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/). 
#' It also uses the same RVE sandwich-type based estimator in a CHE (correlated hierarchical effects) working model with the identical random effects structure as the primary (naive) meta-analytic model. 
#' 
#' We report results for both, PET and PEESE, with the first reported one being the primary (based on the conditional estimator).
#' 
#' **WAAP-WLS**
#' The combined WAAP-WLS estimator (weighted average of the adequately powered - weighted least squares) tries to identify studies that are adequately powered to detect the meta-analytic effect. 
#' If there is less than two such studies, the method falls back to the WLS estimator (Stanley & Doucouliagos, 2015). If there are at least two adequately powered studies, WAAP returns a WLS estimate based on effects from only those studies.
#' 
#' type = 1: WAAP estimate, 2: WLS estimate. kAdequate = number of adequately powered studies
#' 
#' **p-uniform**
#' P-uniform* is a selection model conceptually similar to p-curve. It makes use of the fact that p-values follow a uniform distribution at the true effect size while it includes also nonsignificant effect sizes.
#' Permutation-based version of p-uniform method, the so-called p-uniform* (van Aert, van Assen, 2021).
#' 
#' **p-curve**
#' Permutation-based p-curve method. Output should be self-explanatory. For more info see p-curve.com
#' 
#' **Power for detecting SESOI and bias-corrected parameter estimates**
#' Estimates of the statistical power for detecting a smallest effect sizes of interest equal to .20, .50, and .70 in SD units (Cohen's d). 
#' A sort of a thought experiment, we also assumed that population true values equal the bias-corrected estimates (4/3PSM or PET-PEESE) and computed power for those.
#' 
#' **Handling of dependencies in bias-correction methods**
#' To handle dependencies among the effects, the 4PSM, p-curve, p-uniform are implemented using a permutation-based procedure, randomly selecting only one focal effect (i.e., excluding those which were not coded as being focal) from a single study and iterating nIterations times.
#' Lastly, the procedure selects the result with the median value of the ES estimate (4PSM, p-uniform) or median z-score of the full p-curve (p-curve).
#' 
#' ***

# Descriptives ------------------------------------------------------------

#'# Descriptives
#'
#'## Publication year
c("from" = min(dat$pubYear, na.rm = T), "to" = max(dat$pubYear, na.rm = T))

#'## Sample sizes
#'
#'### N of effects
dat %>% filter(!is.na(yi)) %>% nrow()

#'### N of studies
c("Studies overall" = dat %>% filter(is.na(reasonNotUsed)) %$% length(unique(.$study)),
  "Studies for which ES data were available" = dat %>% filter(!is.na(yi) & is.na(reasonNotUsed)) %$% length(unique(.$study)))

#'###  N of papers
c("Papers overall" = dat %>% filter(is.na(reasonNotUsed)) %$% length(unique(.$paperID)),
  "Papers for which ES data were available" = dat %>% filter(!is.na(yi) & is.na(reasonNotUsed)) %$% length(unique(.$paperID)))

#'### Median N across all the ES eligible for meta-analysis
median(dat$ni, na.rm = T)

#'### Total meta-analytic N
out <- list(NA)
for(i in unique(dat[is.na(dat$reasonNotUsed),]$study)){
  out[i] <- dat %>% filter(study == i & is.na(reasonNotUsed)) %>% select(ni) %>% max()
}
sum(unlist(out), na.rm = T)

#'## Proportion of effects in published papers
prop.table(table(dat$published))*100

#'## Mean gender ratio (percent female)
out <- list(NA)
for(i in unique(dat[is.na(dat$reasonNotUsed),]$study)){
  out[i] <- dat %>% filter(study == i & is.na(reasonNotUsed)) %>% select(percFemale) %>% unlist() %>% median()
}
c("Mean" = mean(unlist(out), na.rm = T), "SD" = sd(unlist(out), na.rm = T))

#'## Weighted mean age of included samples
dat[is.na(dat$reasonNotUsed) & !is.na(dat$ni),] %$% weighted.mean(x = meanAge, w = ni, na.rm = T)

#'## Proportion of sample types
prop.table(table(dat$sampleType))*100

#'## Proportion for GD criteria usage
prop.table(table(dat$gdCriteria))*100

#'## GD measure used
sort(prop.table(table(dat$gdMeasure))*100, decreasing = T)

#'## Correlate type proportions
prop.table(table(dat$correlateType))*100

#'## Reason not used proportions
prop.table(table(dat$reasonNotUsed))*100

#'## Design used proportions
prop.table(table(dat$design))*100

#'## Possible CoI proportions
prop.table(table(dat$possibleCOI))*100

# Meta-analysis -----------------------------------------------------------
#'# Meta-analysis results for individual correlates
#'
#' Number of iterations run equal to `r nIterationsPcurve` for p-curve and `r nIterations` for all other bias correction functions.
#' For sensitivity analyses, we ran `r nIterationVWsensitivity` iterations for the Vevea & Woods (2005) step function model.
#' For supplementary robust bayesian model-averaging approach, we employed `r robmaChains` MCMC chains with `r robmaSamples` sampling iterations.
tab <- sort(table(dat$correlate), decreasing = T)[sort(table(dat$correlate), decreasing = T) > kThreshold]
corVect <- names(tab)
rmaObjects <- rmaResults <- metaResultsPcurve <- vector(mode = "list", length(corVect))
names(rmaObjects) <- names(rmaResults) <- names(metaResultsPcurve) <- corVect
for(i in 1:length(corVect)){
  data <- dat %>% filter(correlate == corVect[i])
  if(outlierSensitivity == TRUE){
    if(!is.null(infResults)){ 
      data <- data %>% filter(!result %in% as.numeric(infResults[[i]]$rowname))
    }
  }
  if(selInfSensitivity == TRUE){
    data <- data %>% filter(!selectiveInference == 1)
  }
  rmaObject <- data %>% rmaCustom()
  rmaResults[[i]] <- data %>% maResults(., rmaObject = rmaObject, bias = biasOn)
  if(biasOn == TRUE){metaResultsPcurve[[i]] <- if(is.list(metaResultsPcurve)){metaResultsPcurve} else {NA}}
  rmaObjects[[i]] <- rmaObject[[1]]
}
(rmaTable <- lapply(rmaResults, function(x){maResultsTable(x, bias = biasOn)}) %$% as.data.frame(do.call(rbind, .)))

#'## Meta-analysis plots (forest, funnel, p-curve plots)
forestPlots <- funnelPlots <- pcurvePlots <- list(NA)
for(i in 1:length(corVect)){
  xlab <- eval(substitute(corVect[i]))
  forest(rmaObjects[[i]], order = order(rmaObjects[[i]]$vi.f, decreasing = T), addpred = T, header="Paper/Study/Effect", xlab = xlab, mlab="", col="gray40")
  forestPlots[[i]] <- recordPlot()
  funnel(rmaObjects[[i]], level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0, pch = 20, yaxis = "sei", digits = c(1, 2), xlab = xlab)
  funnelPlots[[i]] <- recordPlot()
  tryCatch(quiet(pcurveMod(metaResultsPcurve[[i]], effect.estimation = FALSE, plot = TRUE)), error = function(e) NA)
  if(!is.null(metaResultsPcurve[[i]])){title(xlab, cex.main = 1)} else {next}
  pcurvePlots[[i]] <- recordPlot()
  if(displayPlots == FALSE) dev.off()
}
names(forestPlots) <- names(funnelPlots) <- corVect
if(biasOn == TRUE){names(pcurvePlots) <- corVect}

#'# Meta-analysis results for aggregated correlate types

# Empirical direction of individual correlate types
# Needed to scale the risk and protective factor in the same direction on aggregate
tabED <- sort(table(dat$correlate), decreasing = T)[sort(table(dat$correlate), decreasing = T) > 1]
corVectED <- names(tabED)
rmaResultsED <- vector(mode = "list", length(corVectED))
names(rmaResultsED) <- corVectED
for(i in 1:length(corVectED)){
  rmaResultsED[[i]] <- rma(yi, vi, data = dat %>% filter(correlate == corVectED[i]))
}

empiricalDirection <- cbind("correlate" = names(rmaResultsED), 
                            lapply(rmaResultsED, function(x){x[[1]]/abs(x[[1]])}) %$% as.data.frame(do.call(rbind, .)) %>% `colnames<-`("empiricalDirection")) %>% `rownames<-` (NULL)
dat <- left_join(dat, empiricalDirection, by = "correlate")

#'## Comparison of correlate types
#'
#'### Model without covariates
viMatrix <- dat %$% impute_covariance_matrix(vi, cluster = study, r = rho)
rmaObjectCorrType <- dat %>% mutate(yi = yi*empiricalDirection) %$% rma.mv(yi ~ 0 + correlateType, V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(RVEmodelCorrType <- dat %$% list("k" = table(correlateType),
                                  "test" = coef_test(rmaObjectCorrType, vcov = "CR2", test = "z", cluster = study), 
                                  "CIs" = conf_int(rmaObjectCorrType, vcov = "CR2", test = "z", cluster = study),
                                  "RVE Wald test" = Wald_test(rmaObjectCorrType, constraints = constrain_equal(1:5), vcov = "CR2")))

#'### Model with covariates
#' 
#' Controlling for design-related factors that are prognostic w.r.t. the effect sizes (i.e., might vary across moderator categories).
rmaObjectCorrTypeCov <- dat %>% mutate(yi = yi*empiricalDirection) %$% rma.mv(yi ~ 0 + correlateType + meanAge + percFemale + gamingStyle + sampleType, V = viMatrix, method = "REML", random = ~ 1|study/result, sparse = TRUE)
(RVEmodelCorrTypeCov <- dat %$% list("test" = coef_test(rmaObjectCorrTypeCov, vcov = "CR2", test = "z", cluster = study), 
                                     "CIs" = conf_int(rmaObjectCorrTypeCov, vcov = "CR2", test = "z", cluster = study),
                                     "RVE Wald test" = Wald_test(rmaObjectCorrTypeCov, constraints = constrain_equal(1:5), vcov = "CR2")))

#'### Full results for correlate type subsets
tabT <- sort(table(dat$correlateType), decreasing = T)[sort(table(dat$correlateType), decreasing = T) > kThreshold]
corVectT <- names(tabT)
rmaObjectsT <- rmaResultsT <- metaResultsPcurveT <- vector(mode = "list", length(corVectT))
names(rmaObjectsT) <- names(rmaResultsT) <- names(metaResultsPcurveT) <- corVectT
for(i in 1:length(corVectT)){
  data <- dat %>% filter(correlateType == corVectT[i]) %>% mutate(yi = yi*empiricalDirection)
  if(outlierSensitivity == TRUE){
    if(!is.null(infResults)){
      data <- data %>% filter(!result %in% as.numeric(infResults[[i]]$rowname))
    }
  }
  if(selInfSensitivity == TRUE){
    data <- data %>% filter(!selectiveInference == 1)
  }
  rmaObjectT <- data %>% rmaCustom()
  rmaResultsT[[i]] <- data %>% maResults(., rmaObject = rmaObjectT, bias = biasOn)
  if(biasOn == TRUE){metaResultsPcurveT[[i]] <- metaResultPcurve}
  rmaObjectsT[[i]] <- rmaObjectT[[1]]
}
(rmaTableT <- lapply(rmaResultsT, function(x){maResultsTable(x, bias = biasOn)}) %$% as.data.frame(do.call(rbind, .)))

# #'### Meta-analysis plots for aggregated correlate types (forest, funnel, p-curve plots)
# forestPlotsT <- pcurvePlotsT <- list(NA)
# for(i in 1:length(corVectT)){
#   xlab <- eval(substitute(corVectT[i]))
#   forest(rmaObjectsT[[i]], order = order(rmaObjectsT[[i]]$vi.f, decreasing = T), addpred = T, header="Paper/Study/Effect", xlab = xlab, mlab="", col="gray40")
#   forestPlotsT[[i]] <- recordPlot()
#   tryCatch(quiet(pcurveMod(metaResultsPcurveT[[i]], effect.estimation = FALSE, plot = TRUE)), error = function(e) NULL)
#   if(!is.null(metaResultsPcurveT[[i]])){title(xlab, cex.main = 1)} else {next}
#   pcurvePlotsT[[i]] <- recordPlot()
#   if(displayPlots == FALSE) dev.off()
# }
# names(forestPlotsT) <- corVectT
# if(biasOn == TRUE){names(pcurvePlotsT) <- corVectT}

#'## Detailed results for aggregated correlate types
rmaResultsT

#'## Summary forest plot
correlatesES <- lapply(c(rmaResults, rmaResultsT), function(x){cbind(x[[1]]$k.eff,
                                                     x[[2]]$test$beta,
                                                     x[[2]]$CIs$CI_L,
                                                     x[[2]]$CIs$CI_U)}) %$% as.data.frame(do.call(rbind, .))
rownames(correlatesES) <- c(rownames(rmaTable), rownames(rmaTableT))
correlatesES <- rownames_to_column(correlatesES)
colnames(correlatesES) <- c("Correlate", "k", "yi", "ciLB", "ciUB")

par(mar=c(4,4,1,2))
correlatesES %$% 
  forest(yi, ci.lb = ciLB, ci.ub = ciUB, slab = Correlate, ilab = k, xlim = c(-1, 1), alim = c(-1.3,1.3),
         ilab.xpos = .65, ilab.pos = 4, xlab="Meta-analytic effect size (Pearson's r)", header = "Risk or protective factor",
         at=c(-.5, -.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7), cex = .8,
         rows=c(38:10, 7:1),
         ylim=c(1, 43)
         )
text(.707, 39.5, "k", font = 2, cex = .8)
text(-1.299, c(8.5,39.5), font = 2, cex = .8, pos = 4, c("Aggregate correlate types", "Individual correlate types"))

# Moderators of GD  -----------------------------------------------------

#'# Moderation analyses
#'
#'## Substantive Moderators of GD

mods <- dat %>% select(c(percFemale, meanAge, sampleType, gdCriteria)) %>% names()
rmaModObjects <- rmaModTest <- vector(mode = "list", length(mods))
names(rmaModObjects) <- names(rmaModObjects) <- mods
for(i in 1:length(mods)){
  for(j in 1:length(corVect)){
    data <- dat %>% filter(correlate == corVect[j] & !is.na(yi))
    data <- data %>% mutate(percFemale = scale(as.numeric(percFemale)),
                            meanAge = scale(as.numeric(meanAge)),
                            sampleTypeBinary = factor(dplyr::recode(sampleType, `3` = 1, `4` = 1, `1` = 2, `2` = 2, `5` = NA_real_, `6` = NA_real_)),
                            gdCriteriaBinary = factor(dplyr::recode(gdCriteria, `1` = 1, `2` = 1, `3` = 2, `4` = 1, `5` = NA_real_)),
                            yi = yi_z, 
                            vi = vi_z)
    if(outlierSensitivity == TRUE){
      if(!is.null(infResults)){
        data <- data %>% filter(!result %in% as.numeric(infResults[[i]]$rowname))
      }
    }
    if(selInfSensitivity == TRUE){
      data <- data %>% filter(!selectiveInference == 1)
    }
    viMatrix <- data %$% impute_covariance_matrix(vi, cluster = study, r = rho)
    rmaModObject <- tryCatch(rma.mv(data$yi, V = viMatrix, mods = ~ data[,mods[i]], method = "REML", random = ~ 1|data[,"study"]/data[,"result"], sparse = TRUE), error = function(e) NULL)
    rmaModObjects[[i]][[j]] <- rmaModObject
    rmaModTest[[i]][[j]] <- c("Mod estimate" = round(rmaModObject$beta[2], 2), "Qm stat" = round(rmaModObject$QM, 2), "df" = as.integer(rmaModObject$QMdf[1]), "p" = round(rmaModObject$QMp, 3))
  }
  names(rmaModObjects[[i]]) <- names(rmaModTest[[i]]) <- corVect
}
names(rmaModTest) <- mods
(modResults <- lapply(as.data.frame(do.call(rbind, rmaModTest)), function(x){as.data.frame(x)}))

#'## Heatplot for substantive moderators

#' Heatplot for metric + dichotomized moderators based on effect size magnitudes
#' The first four correlates are the ones found act to protectively.

modHeatmapData <- cbind(
  do.call(rbind, lapply(modResults, function(x)x[1,])) %>% `colnames<-`(c("% female", "Mean age", "Sample type", "GD criteria")),
  do.call(rbind, lapply(modResults, function(x)x[4,])) %>% `colnames<-`(c("percFemaleP", "meanAgeP", "Sample type p", "GD criteria p")))
modHeatmapData <- modHeatmapData[c(1, 8, 17, 19, 24, 29, 2:7, 9:16, 18, 20:23, 25:28),]
modHeatmapData <- modHeatmapData %>% mutate(cellnotePercFemale = paste(modHeatmapData[,1], " (", ifelse(sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,5], 3))) == "", paste("<.001"), sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,5], 3)))), ")", sep = ""),
                                            cellnoteMeanAge = paste(modHeatmapData[,2], " (", ifelse(sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,6], 3))) == "", paste("<.001"), sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,6], 3)))), ")", sep = ""),
                                            cellnoteSampleType = paste(modHeatmapData[,3], " (", ifelse(sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,7], 3))) == "", paste("<.001"), sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,7], 3)))), ")", sep = ""),
                                            cellnoteGDcriteria = paste(modHeatmapData[,4], " (", ifelse(sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,8], 3))) == "", paste("<.001"), sub("^0+", "", sprintf("%.3f", round(modHeatmapData[,8], 3)))), ")", sep = ""))
modHeatmapData <- rbind(NA, modHeatmapData[1:6,], NA, modHeatmapData[7:nrow(modHeatmapData),])
rownames(modHeatmapData)[1] <- "––– Protective factors –––"
rownames(modHeatmapData)[8] <- "––– Risk factors ––––––––"
modHeatmapData[2:7,1:4] <- modHeatmapData[2:7,1:4]*-1 # Invert the protective factors so that the stronger positive relationship between the moderator and the outcome, the more red the color in the heatplot.
modHeatmapData[1,9:10] <- "β (p-value)"
modHeatmapData[1,11:12] <- "B (p-value)"
heatmap.2(as.matrix(modHeatmapData[,1:4]), 
          cellnote = modHeatmapData[,9:12], notecex = 1, cexCol = 1.2, srtCol = 45, cexRow = 1.2, key = FALSE, notecol = "black", dendrogram = "none", Rowv = FALSE, Colv = FALSE, trace = "none", col = bluered, tracecol = "#303030",
          lmat = rbind(c(0, 3), c(1, 2), c(2, 4)), lhei = c(.1, 10, 1), lwid = c(.5, .2))

# Methodological moderators -----------------------------------------------

#'## Methodological moderators

#'### F-test of equality of variances for designs employing restricted vs unrestricted samples
#'
#'#### Mean vi for restricted designs
(viRestricted <- dat %>% filter(sampleRestricted == 1) %$% mean(vi, na.rm = T))
dfRestricted <- dat %>% filter(sampleRestricted == 1, !is.na(vi)) %>% nrow() - 1
#'#### Mean vi for non-restricted designs
(viNonRestricted <- dat %>% filter(sampleRestricted == 0) %$% mean(vi, na.rm = T))
dfNonRestricted <- dat %>% filter(sampleRestricted == 0, !is.na(vi)) %>% nrow() - 1
#'#### F-test
(testEqVar <- c("F" = viRestricted/viNonRestricted, "p" = 1 - pf(viRestricted/viNonRestricted, df1 = dfRestricted, df2 = dfNonRestricted)))

#'### Year of Publication
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEpubYear <- summary(lmer(scale(sqrt(vi)) ~ scale(journalH5) + scale(pubYear) + (1|study), data = dat, REML = T))$coefficients)
#' Comment: all the variables were centered for easier interpretation of model coefficients. 
#'

#'#### Scatterplot year <-> precision
#'
#' Size of the points indicate the H5 index (the bigger the higher) of the journal that the ES is published in.
(yearPrecisionScatter <- dat %>%  ggplot(aes(pubYear, sqrt(vi))) + 
    geom_point(aes(size = journalH5), alpha = .70, colour = "#80afce") +
    geom_smooth(method = lm) +
    scale_x_continuous() +
    xlab("Year of publication") +
    ylab("Imprecision") +
    theme_bw() +
    theme(legend.position = "none"))

#'### Citations
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEcitations <- summary(lmer(scale(sqrt(vi)) ~ scale(pubYear) + scale(journalH5) + scale(citations) + (1|study), data = dat, REML = T))$coefficients)

#'#### Scatterplot precision <-> citations
#'
#' The relationship between precision (sqrt of variance) and number of citations (log).
(citationsPrecisionScatter <- dat %>% ggplot(aes(log(citations + 1), sqrt(vi))) + 
    geom_point(alpha = .70, colour = "#80afce") +
    geom_smooth(method = lm) +
    xlab("Citations (log)") +
    ylab("Precision") +
    theme_bw() +
    theme(legend.position = "none"))

#'### H5 index
#'
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
(LMEh5 <- summary(lmer(scale(sqrt(vi)) ~ scale(journalH5) + (1|study), data = dat, REML = T))$coefficients)

#'#### Scatterplot precision <-> journal H5
#'
#' The relationship between precision (sqrt of variance) and H5 index of the journal.
(h5PrecisionScatter <- dat %>% ggplot(aes(journalH5, sqrt(vi))) + 
    geom_point(alpha = .70, colour = "#80afce") +
    geom_smooth(method = lm) +
    xlab("H5 index") +
    ylab("Precision") +
    theme_bw() +
    theme(legend.position = "none"))

#'### Decline effect
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study.
meanAbsYi <- as.numeric(rma.mv(abs(yi), vi, data = dat, method = "REML", random = ~ 1|study/result)[1])
(declineEff <- summary(lmer(scale(meanAbsYi - abs(yi)) ~ scale(sqrt(vi)) + scale(pubYear) + (1|study), data = dat))$coefficients)

#'### Citation bias
#' Do more highly-cited studies report larger effect sizes?
(LMEcitationsYi <- summary(lmer(abs(yi) ~ scale(pubYear) + scale(journalH5) + scale(citations) + (1|study), data = dat, REML = T))$coefficients)
#' Additionally taking into account precision of the effect (vi)
(LMEcitationsYiVi <- summary(lmer(abs(yi) ~ scale(pubYear) + scale(journalH5) + scale(citations) + scale(vi) + (1|study), data = dat, REML = T))$coefficients)

# Sensitivity analyses ----------------------------------------------------

#'# Sensitivity analyses

#'## Corrections of statistical artifacts 
#'
#' Corrections for measurement error and selection effects
#' Indirect selection bias not adjusted for intelligence due to low k + missing reliability data
xyArtefactCorResult <- xArtefactCorResult <- vector(mode = "list", length(corVect))
names(xyArtefactCorResult) <- names(xArtefactCorResult) <- corVect
#+ include == FALSE
for(i in 1:length(corVect)){
  artefactCorObj <- suppressMessages(ma_r(ma_method = ifelse(indirectSel == TRUE, "ic", "ad"),rxyi = yi, n = ni, sample_id = label, 
                         rxx = rxxV1sample, ryy = rxxV2sample, ux = ux, uy = uy,
                         correct_rxx = TRUE, correct_ryy = TRUE, correct_rr_x = TRUE, correct_rr_y = TRUE, 
                         indirect_rr_x = indirectSel, indirect_rr_y = indirectSel,
                         data = dat %>% filter(correlate == corVect[i] & !is.na(yi) & !is.na(ni))))
  xyArtefactCorResult[[i]] <- data.frame(artefactCorObj$meta_tables[[1]][[ifelse(indirectSel == TRUE, 2, 3)]]$true_score)[c("k", "N", "mean_r", "mean_rho", "sd_rho", "CI_LL_95", "CI_UL_95")]
  xArtefactCorResult[[i]] <- data.frame(artefactCorObj$meta_tables[[1]][[ifelse(indirectSel == TRUE, 2, 3)]]$validity_generalization_x)[c("k", "N", "mean_r", "mean_rho", "sd_rho", "CI_LL_95", "CI_UL_95")]
}

#+ include == TRUE
#'#### Full results
(artCorResult <- list("Correcting for measurement error in both, GD and correlate" = do.call(rbind, xyArtefactCorResult),
                      "Correcting for measurement error only in GD" = do.call(rbind, xArtefactCorResult)))
#'#### Just the estimates + delta
cbind(rownames_to_column(artCorResult[[2]])[1], r = abs(artCorResult[[2]]$mean_r), rho = abs(artCorResult[[2]]$mean_rho), diff = abs(artCorResult[[2]]$mean_r) - abs(artCorResult[[2]]$mean_rho))

#'## Using Fisher's z instead of untransformed r
rmaObjectsZ <- rmaResultsZ <- vector(mode = "list", length(corVect))
names(rmaObjectsZ) <- names(rmaResultsZ) <- corVect
for(i in 1:length(corVect)){
  rmaObjectsZ[[i]] <- predict((dat %>% filter(correlate == corVect[i]) %>% mutate(yi = yi_z, vi = vi_z) %>% rmaCustom())[[1]], transf = transf.ztor)
}
do.call(rbind, rmaObjectsZ)[,c(1,3:8)]

#'## Numerical inconsistencies in reported p-values
#'
#'#### how many results were analyzed
nrow(statcheck) 
#'#### how many papers reported results in APA format
length(unique(statcheck$Source))/length(unique(dat[!is.na(dat$yi),]$paperID))
#'#### how many % statcheck errors
prop.table(table(statcheck$Error))*100
#'#### how many % statcheck errors affected the decision
table(statcheck$DecisionError)[2]/table(statcheck$Error)["TRUE"]*100
#'#### How many % papers contained statcheck errors
statcheck %>% filter(Error == TRUE) %>% select(Source) %>% unique() %>% nrow()/length(unique(statcheck$Source))*100

#'## p-curve for the full literature
(pCurveFull <- pcurvePerm(dat))
pcurveMod(metaResultPcurve, effect.estimation = FALSE, plot = TRUE)
title("p-curve for the full literature", cex.main = 1)

#'## Power for the full literature
(powerFull <- powerEst(dat, forBiasAdj = FALSE))
  
# Save workspace
save.image("workspace.RDS")

#' Session info
sessionInfo()

#' Script running time
endTime <- Sys.time()
endTime - startTime
