# Load libraries (and install if not installed already)
list.of.packages <- c("car", "reshape", "tidyverse", "tidyr", "psych", "metafor", "meta", "psychmeta", "dmetar", "esc", "lme4", "ggplot2", "knitr", "puniform", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp", "magrittr", "weightr", "clubSandwich", "ddpcr", "poibin", "robvis", "RoBMA", "gplots")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
select <- dplyr::select
funnel <- metafor::funnel
forest <- metafor::forest

# Determine the alpha level / use of one-tailed vs two-tailed test for p-uniform and PET-PEESE
if (test == "one-tailed") {
  alpha <- .10
} else if (test == "two-tailed") {
  alpha <- .05
}

# Meta-analysis -----------------------------------------------------------

# Meta analysis run on a filtered dataset

# Custom robust (RVE) multivariate RE meta-analytic model (using the CHE working model)
# Needs specific naming of ES, variances and data on clustering; yi = yi, vi = vi, study, result
# Using n/(n-p) small-sample correction for RVE SEs
rmaCustom <- function(data = NA){
  data <- data %>% filter(useMeta == 1)
  viMatrix <- data %$% impute_covariance_matrix(vi, cluster = study, r = rho)
  rmaObjectModBasedSE <- rma.mv(yi = yi, V = viMatrix, data = data, method = "REML", random = ~ 1|study/result, sparse = TRUE, slab = label)
  rmaObject <- robust.rma.mv(rmaObjectModBasedSE, cluster = data$study)
  return(list("RMA.MV object with RVE SEs with n/(n-p) small-sample correction" = rmaObject, 
              "RMA.MV object without cluster robust SEs" = rmaObjectModBasedSE))
}

# 95% prediction interval -------------------------------------------------
pi95 <- function(rmaObject = NA){
  pi95Out <- c("95% PI LB" = round(predict.rma(rmaObject[[1]])$cr.lb, 3), "95% PI UB" = round(predict.rma(rmaObject[[1]])$cr.ub, 3))
  pi95Out
}

# Heterogeneity -----------------------------------------------------------

heterogeneity <- function(rmaObject = NA){
  
  # Total heterogeneity - tau
  tau <- sqrt(sum(rmaObject[[1]]$sigma2))
  
  # I^2
  W <- diag(1/rmaObject[[1]]$vi)
  X <- model.matrix(rmaObject[[1]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2<- 100 * sum(rmaObject[[1]]$sigma2) / (sum(rmaObject[[1]]$sigma2) + (rmaObject[[1]]$k - rmaObject[[1]]$p)/sum(diag(P)))
  
  # Separate estimates of between- and within-cluster heterogeneity
  BW.hetero <- round(100 * rmaObject[[1]]$sigma2 / (sum(rmaObject[[1]]$sigma2) + (rmaObject[[1]]$k - rmaObject[[1]]$p)/sum(diag(P))), 2)
  
  studyID <- rmaObject[[1]]$mf.r[[1]]$study
  resultID <- rmaObject[[1]]$mf.r[[1]]$result
  resR <- rma.mv(yi = rmaObject[[1]]$yi, V = rmaObject[[1]]$vi, random = ~ 1|studyID/resultID)
  resF <- rma.mv(yi = rmaObject[[1]]$yi, V = rmaObject[[1]]$vi)
  
  # Jackson's approach to I^2
  JI2 <- round(c(100 * (vcov(resR)[1,1] - vcov(resF)[1,1])/vcov(resR)[1,1]), 2)
  
  # Intra-class correlation of underlying true effects
  icc <- round(rmaObject[[1]]$sigma2[1]/sum(rmaObject[[1]]$sigma2), 2)
  
  c("Tau" = tau,
    "I^2" = I2,
    "Jackson's I^2" = JI2,
    "Between-cluster heterogeneity" = BW.hetero[1],
    "Within-cluster heterogeneity" = BW.hetero[2],
    "ICC" = icc)
}

# Proportion of significant results ---------------------------------------

propSig <- function(p.values = NA){
  as.integer(table(p.values < .05)["TRUE"])/length(p.values < .05)
}

# Permutation p-curve -----------------------------------------------------

# Subseting only the effects that are focal for the published study
pcurvePerm <- function(data, esEstimate = FALSE, plot = FALSE, nIterations = nIterationsPcurve){
  resultIDpcurve <- list(NA)
  resultPcurve <- matrix(ncol = 11, nrow = nIterationsPcurve)
  set.seed(1)
  for(i in 1:nIterationsPcurve){
    datPcurve <- data[!duplicated.random(data$study) & data$focal == 1 & !is.na(data$yi) & !is.na(data$vi),]
    metaPcurve <- tryCatch(metagen(TE = yi, seTE = sqrt(vi), n.e = ni, data = datPcurve),
                           error = function(e) NULL)
    modelPcurve <- tryCatch(pcurveMod(metaPcurve, effect.estimation = esEstimate, N = datPcurve$ni, plot = plot), 
                            error = function(e) NULL)
    if(is.null(modelPcurve)){
      next
    } else {
      resultPcurve[i,] <- c("iterationNo" = i, "rightskew" = modelPcurve$pcurveResults[1,], "flatness" = modelPcurve$pcurveResults[2,])  
    }
    resultIDpcurve[[i]] <- datPcurve$result
  }
  colnames(resultPcurve) <- c("iterationNo", "rightskew.pBinomial", "rightskew.zFull", "rightskew.pFull", "rightskew.zHalf", "rightskew.pHalf", "flatness.pBinomial", "flatness.zFull", "flatness.pFull", "flatness.zHalf", "flatness.pHalf")
  medianResultPcurve <- resultPcurve %>% data.frame() %>% na.omit() %>% arrange(rightskew.zFull) %>% slice(ceiling(n()/2)) %>% unlist()
  metaResultPcurve <- tryCatch(metagen(TE = yi, seTE = sqrt(vi), data = data[data$result %in% unlist(resultIDpcurve[medianResultPcurve["iterationNo"]]),]), 
                               error = function(e) NULL)
  metaResultPcurve <<- metaResultPcurve
  pcurveMod(metaResultPcurve, effect.estimation = esEstimate, plot = plot)
}

# Multiple-parameter selection models -------------------------------------
# 4/3-parameter selection model (4PSM/3PSM)
selectionModel <- function(data, minNoPvals = minPvalues, nIteration = nIterations, fallback = FALSE, steps = c(.025, 1), deltas = NA){
  data <- data %>% filter(useMeta == 1)
  resultSM <- matrix(ncol = 8, nrow = nIteration)
  set.seed(1)
  for(i in 1:nIteration){
    dataSM <<- data[!duplicated.random(data$study) & data$focal == 1,]
    res <- tryCatch(rma(yi, vi,  data = dataSM), error = function(e) NULL)
    pTable <- table(cut(dataSM$p, breaks = c(0, .05, 0.5, 1)))
    if(fallback == TRUE | any(pTable < minNoPvals) | !anyNA(deltas)){
      threeFit <- tryCatch(selmodel(res, type = "stepfun", steps = steps, delta = deltas, alternative = "greater"),
                           error = function(e) NULL)
      threeOut <- if(is.null(threeFit)){
        next
      } else {
        round(c("est" = threeFit$beta, "se" = threeFit$se, "zvalue" = threeFit$zval, "pvalue" = threeFit$pval, "ciLB" = threeFit$ci.lb, "ciUB" = threeFit$ci.ub, "k" = threeFit$k, "steps" = length(threeFit$steps)), 3)
      }
      out <- threeOut
    } else { 
      fourFit <- tryCatch(selmodel(res, type = "stepfun", steps = c(.025, .5, 1), alternative = "greater"),
                          error = function(e) NULL)
      fourOut <- c("est" = fourFit$beta, "se" = fourFit$se, "zvalue" = fourFit$zval, "pvalue" = fourFit$pval, "ciLB" = fourFit$ci.lb, "ciUB" = fourFit$ci.ub, "k" = fourFit$k, "steps" = length(fourFit$steps))  
      if (is.null(fourFit)){
        out <- threeOut
      } else {
        out <- fourOut %>% round(., 3)
      }
    }
    resultSM[i,] <- out  
  }
  colnames(resultSM) <- c("est", "se", "zvalue", "pvalue", "ciLB", "ciUB", "k", "steps")
  resultSM <- resultSM %>% data.frame() %>% na.omit() %>% arrange(est) %>% slice(ceiling(n()/2)) %>% unlist()
  resultSM <<- resultSM
  resultSM
}

# Vevea & Woods step function model using a priori defined selection weights
veveaWoodsSM <- function(data, stepsDelta, nIteration = nIterationVWsensitivity){
  dataVW <- data[data$focal == 1 & data$useMeta == 1,]
  set.seed(1)
  tryCatch(do.call(rbind, lapply(stepsDelta[-1], function(delta) suppressWarnings(selectionModel(dataVW, steps = stepsDelta$steps, deltas = delta, nIteration = nIterationVWsensitivity)))),
           error = function(e) NULL)
}

# Robust Bayesian meta-analysis
bma <- function(data, seedNo = 1, chainsNo = robmaChains, nIterationBMA = robmaSamples){
  dataBMA <- data[data$focal == 1 & data$useMeta == 1,]
  set.seed(1)
  tryCatch(summary(dataBMA %>% filter(!is.na(yi) & !is.na(ni)) %$% RoBMA(r = yi, n = ni, study_names = label, seed = seedNo,
                         chains = chainsNo, sample = nIterationBMA, burnin = ifelse(2*nIterationBMA/5 < 50, 50, 2*nIterationBMA/5), parallel = TRUE)), 
           error = function(e) NULL)
}

# PET-PEESE ---------------------------------------------------------------

#PET-PEESE with 4/3PSM as the conditional estimator instead of PET. 
# Also implemented the modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/).
petPeese <- function(data, nBased = TRUE, selModAsCondEst = condEst){  # if nBased = TRUE, use the sample-size-based estimator, if FALSE, use the ordinary SE/var. If selModAsCondEst = TRUE, use the selection model as conditional estimator, otherwise use PET.
  data <- data %>% filter(useMeta == 1)
  viMatrix <- data %$% impute_covariance_matrix(vi, cluster = study, r = rho)  # compute the covariance matrix for the CHE working model
  
  if(nBased == TRUE){
    pet <<- robust.rma.mv(rma.mv(yi = yi ~ sqrt(nTerm), V = viMatrix, random = ~ 1|study/result, method = "REML", sparse = TRUE, data = data), cluster = data$study)
  } else {
    pet <<- robust.rma.mv(rma.mv(yi = yi ~ sqrt(vi), V = viMatrix, random = ~ 1|study/result, method = "REML", sparse = TRUE, data = data), cluster = data$study)
  }
  pet.out <- round(c(pet$b[1], pet$se[1], pet$zval[1], pet$pval[1], pet$ci.lb[1], pet$ci.ub[1]), 3)
  names(pet.out) <- c("PET estimate", "se", "zvalue", "pvalue", "ciLB", "ciUB")
  
  if(nBased == TRUE){
    peese <<- robust.rma.mv(rma.mv(yi = yi ~ nTerm, V = viMatrix, random = ~ 1|study/result, method = "REML", sparse = TRUE, data = data), cluster = data$study)
  } else {
    peese <<- robust.rma.mv(rma.mv(yi = yi ~ vi, V = viMatrix, random = ~ 1|study/result, method = "REML", sparse = TRUE, data = data), cluster = data$study)
  } 
  peese.out <- round(c(peese$b[1], peese$se[1], peese$zval[1], peese$pval[1], peese$ci.lb[1], peese$ci.ub[1]), 3)
  names(peese.out) <- c("PEESE estimate", "se", "zvalue", "pvalue", "ciLB", "ciUB")
  
  if(selModAsCondEst == TRUE){
    ifelse(resultSM["pvalue"] < alpha & ifelse(median(data$directionEffect, na.rm = T) == -1, -1, 1) * resultSM["est"] > 0,
           return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  } else {
    ifelse(pet$pval[1] < alpha & ifelse(median(data$directionEffect, na.rm = T) == -1, -1, 1) * pet$b[1] > 0,
           return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  
  # if(selModAsCondEst == TRUE){
  #   ifelse(resultSM["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultSM["est"] > 0,
  #          return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  # } else {
  #   ifelse(pet$pval[1] < alpha & ifelse(exists("side") & side == "left", -1, 1) * pet$b[1] > 0,
  #          return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  } 
}

# WAAP-WLS estimator ------------------------------------------------------

# Code by Felix Schonbrodt and Evan Carter; https://github.com/nicebread/meta-showdown/blob/master/MA-methods/6-WAAP.R
# based on Stanley, T. D., Doucouliagos, H., & Ioannidis, J. P. A. (2017). Finding the power to reduce publication bias. Statistics in Medicine, 54(3), 30–19. http://doi.org/10.1002/sim.7228

# WLS estimator

# Weighted least squares estimator: Stanley, T. D., & Doucouliagos, H. (2015). Neither fixed nor random: weighted least squares meta-analysis. Statistics in Medicine, 34(13), 2116–2127. http://doi.org/10.1002/sim.6481
WLS.est <- function(yi, vi, long = FALSE) {
  se <- sqrt(vi)
  yi.precision <- 1/se
  yi.stand <- yi/se
  l1 <- lm(yi.stand ~ 0 + yi.precision)
  s1 <- summary(l1)
  res <- data.frame(
    method = "WLS",
    term = "b0",
    estimate = 	s1$coefficients[1, 1],
    std.error = s1$coefficients[1, 2],
    statistic = s1$coefficients[1, 2],
    p.value = s1$coefficients[1, 4],
    conf.low = confint(l1)[1],
    conf.high = confint(l1)[2]
  )
  returnRes(res, long = FALSE)
}

# return object: type = 1: WAAP, type = 2: WLS (must be numeric, otherwise it distorts the structure of the results object)
waapWLS <- function(yi, vi, est = c("WAAP-WLS"), long = FALSE) {
  
  # 1. determine which studies are in the top-N set
  
  # FE (or, equivalently, WLS) as the proxy for ‘true’ effect.
  WLS.all  <- WLS.est(yi, vi, long=FALSE)
  true.effect <- WLS.all$estimate
  
  # only select studies that are adequatly powered (Stanley uses a two-sided test)
  powered <- true.effect/2.8 >= sqrt(vi)
  
  # 2. compute the unrestricted weighted average (WLS) rma of either all or only adequatly powered studies
  # Combined estimator: WAAP-WLS	
  kAdequate <- sum(powered,na.rm=T)
  
  if (kAdequate >= 2) {
    res <- WLS.est(yi[powered], vi[powered], long=FALSE)
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type=1, kAdequate=kAdequate))
  } else {
    res <- WLS.all
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type=2, kAdequate=kAdequate))
  }
  returnRes(res, long = FALSE)
}

# Median power for detecting SESOI and bias-corrected parameter estimates --------------

powerEst <- function(data = NA, forBiasAdj = TRUE){
  data <- data %>% filter(useMeta == 1)
  powerPEESE <- NA
  powerSM <- NA
  peeseEst <- petPeese(data)[1]
  power10r <- round(median(pwr::pwr.r.test(n = data[!is.na(data$ni),]$ni, r = .10)$power), 3)
  power30r <- round(median(pwr::pwr.r.test(n = data[!is.na(data$ni),]$ni, r = .30)$power), 3)
  power50r <- round(median(pwr::pwr.r.test(n = data[!is.na(data$ni),]$ni, r = .50)$power), 3)
  if(forBiasAdj == TRUE){
    powerPEESEresult <- round(median(pwr::pwr.r.test(n = data[!is.na(data$ni),]$ni, r = peeseEst)$power), 3)
    powerSMresult <- round(median(pwr::pwr.r.test(n = data[!is.na(data$ni),]$ni, r = resultSM["est"])$power), 3)
  }
  data.frame("Median power for detecting a SESOI of r = .10" = power10r,
    "Median power for detecting a SESOI of r = .30" = power30r,
    "Median power for detecting a SESOI of r = .50" = power50r,
    "Median power for detecting PET-PEESE estimate" = ifelse(forBiasAdj == TRUE, powerPEESEresult, "Not calculated"), 
    "Median power for detecting 4/3PSM estimate" = ifelse(forBiasAdj == TRUE, powerSMresult, "Not calculated"))
  
    # "Median power for detecting PET-PEESE estimate" = ifelse(forBiasAdj == TRUE, ifelse(peeseEst > 0, powerPEESEresult, paste("ES estimate in the opposite direction")), "Not calculated"), 
    # "Median power for detecting 4/3PSM estimate" = ifelse(forBiasAdj == TRUE, ifelse(resultSM["est"] > 0, powerSMresult, paste("ES estimate in the opposite direction")), "Not calculated"))
}

# Publication bias summary function-------------------------------

bias <- function(data = NA, rmaObject = NA, runRobMA = 1){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject[[1]]$yi, sqrt(rmaObject[[1]]$vi), method = "kendall")
  
  # Small-study effects correction
  # Vevea & Woods selection model
  resultsVeveaWoodsSM <- veveaWoodsSM(data, stepsDelta)
  
  # Robust Bayesian model-averaging approach
  bmaMod <- if(runRobMA == TRUE){bma(data)}
  
  # 3-parameter selection model
  resultSelModel <- selectionModel(data, minNoPvals = minPvalues, nIteration = nIterations, fallback = fallback)
  
  # PET-PEESE
  petPeeseOut <- petPeese(data)
  
  # WAAP-WLS
  waapWLSout <- data %>% filter(useMeta == 1) %$% waapWLS(yi, vi)
  waapWLSout[1, 9:10] <- waapWLSout[2, 9:10]
  waapWLSout <- waapWLSout[1,]
  
  # Permutation p-curve
  pcurvePermOut <- tryCatch(pcurvePerm(data, esEstimate = FALSE, plot = FALSE), error = function(e) NULL)
  
  # p-uniform* (van Aert & van Assen, 2021)
  resultPuniform <- matrix(ncol = 4, nrow = nIterations)
  set.seed(1)
  for(i in 1:nIterations){
    modelPuniform <- data[!duplicated.random(data$study) & data$focal == 1 & data$useMeta == 1 & !is.na(data$yi) & !is.na(data$vi),] %$% tryCatch(puni_star(yi = yi, vi = vi, alpha = alpha, side = ifelse(median(directionEffect, na.rm = T) == -1, "left", "right"), method = "ML"), error = function(e) NULL)
    resultPuniform[i,] <- ifelse(!is.null(modelPuniform), c("est" = modelPuniform[["est"]], "ciLB" = modelPuniform[["ci.lb"]], "ciUB" = modelPuniform[["ci.ub"]], "p-value" = modelPuniform[["pval.0"]]), NA)
  }
  colnames(resultPuniform) <- c("est", "ciLB", "ciUB", "pvalue")
  puniformOut <- resultPuniform %>% data.frame() %>% na.omit() %>% arrange(est) %>% slice(ceiling(n()/2)) %>% unlist()
  puniformOut <<- resultPuniform
  
  # Publication bias results
  return(list("ES-precision correlation" = esPrec,
              "4/3PSM" = resultSelModel,
              "Vevea & Woods SM" = resultsVeveaWoodsSM,
              "Robust BMA" = bmaMod$estimates,
              "PET-PEESE" = petPeeseOut,
              "WAAP-WLS" = waapWLSout,
              "p-uniform*" = puniformOut,
              "p-curve" = pcurvePermOut))
}

# Summary results ---------------------------------------------------------

maResults <- function(rmaObject = NA, data = NA, bias = T){
  list(
    "RMA results with model-based SEs" = rmaObject[[2]],
    "RVE SEs with Satterthwaite small-sample correction" = list("test" = coef_test(rmaObject[[2]], vcov = "CR2", cluster = data[data$useMeta == 1,]$study), "CIs" = conf_int(rmaObject[[2]], vcov = "CR2", cluster = data[data$useMeta == 1,]$study)),
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "Proportion of significant results" = propSig(data[data$useMeta == 1,]$p),
    "Publication bias" = if(bias ==T) {bias(data, rmaObject)} else {paste("Publication bias corrections not carried out")},
    "Power for detecting SESOI and bias-corrected parameter estimates" = if(bias ==T) {powerEst(data)} else {paste("Power for detecting bias-corrected parameter estimates not computed")})
}

biasResults <- function(rmaObject = NA, data = NA){
  list(
    "Publication bias" = bias(data, rmaObject),
    "Power for detecting SESOI and bias-corrected parameter estimates" = powerEst(data))
}

maResultsTable <- function(maResultsObject, metaAnalysis = TRUE, bias = TRUE){
  if(bias == TRUE & metaAnalysis == TRUE){
    noquote(c(
      "k" = as.numeric(maResultsObject[[1]]$k.all),
      "ES [95% CI]" = paste(round(as.numeric(maResultsObject[[2]]$test$beta), 2), " [", round(maResultsObject[[2]]$CIs$CI_L, 2), ", ", round(maResultsObject[[2]]$CIs$CI_U, 2), "]", sep = ""),
      "95% PI [LB, UB]" = paste("[", round(maResultsObject$`Prediction interval`[1], 2), ", ", round(maResultsObject$`Prediction interval`[2], 2), "]", sep = ""),
      "SE" = round(maResultsObject[[2]]$test$SE, 2),
      round(maResultsObject[[4]]["Tau"], 2),
      "I^2" = paste(round(maResultsObject[[4]]["I^2"], 0), "%", sep = ""),
      "% significant" = paste(round(maResultsObject[["Proportion of significant results"]]*100, 0), "%", sep = ""),
      "k for SMs" = maResultsObject[[6]][["4/3PSM"]]["k"],
      "3PSM est [95% CI]" = paste(round(maResultsObject[[6]][["4/3PSM"]]["est"], 2), " [", round(maResultsObject[[6]][["4/3PSM"]]["ciLB"], 2), ", ", round(maResultsObject[[6]][["4/3PSM"]]["ciUB"], 2), "]", sep = ""),
      "3PSM" = round(maResultsObject[[6]][["4/3PSM"]]["pvalue"], 3),
      "V&W [moderate/severe/extreme]" = if(is.numeric(maResultsObject[[6]][["Vevea & Woods SM"]])){paste(round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][1,1]), 2),"/",round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][2,1]), 2), "/", round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][3,1]), 2), sep = "")} else {NA},
      "RoBMA [95% CI]" = if(!is.null(maResultsObject[[6]][["Robust BMA"]])){paste(round(maResultsObject[[6]][["Robust BMA"]]$Mean[1], 2), " [", round(maResultsObject[[6]][["Robust BMA"]]$`0.025`[1], 2), ", ", round(maResultsObject[[6]][["Robust BMA"]]$`0.975`[1], 2), "]", sep = "")} else {NA},
      "PET-PEESE est [95% CI]" = paste(round(maResultsObject[[6]][["PET-PEESE"]][1], 2), " [", round(maResultsObject[[6]][["PET-PEESE"]][5], 2), ", ", round(maResultsObject[[6]][["PET-PEESE"]][6], 2), "]", sep = ""),
      "PET-PEESE" = round(maResultsObject[[6]][["PET-PEESE"]][4], 3),
      "Power to detect small SESOI" = maResultsObject[[7]][[1]],
      "Power to detect medium SESOI" = maResultsObject[[7]][[2]]
    ))
  } else if (metaAnalysis == FALSE & bias == TRUE) {
    noquote(c(
      "% significant" = paste(round(maResultsObject[["Proportion of significant results"]]*100, 0), "%", sep = ""),
      "k for SMs" = maResultsObject[[6]][["4/3PSM"]]["k"],
      "3PSM est [95% CI]" = paste(round(maResultsObject[[1]][["4/3PSM"]]["est"], 2), " [", round(maResultsObject[[1]][["4/3PSM"]]["ciLB"], 2), ", ", round(maResultsObject[[1]][["4/3PSM"]]["ciUB"], 2), "]", sep = ""),
      "3PSM" = round(maResultsObject[[1]][["4/3PSM"]]["pvalue"], 3),
      "V&W [moderate/severe/extreme]" = if(is.numeric(maResultsObject[[6]][["Vevea & Woods SM"]])){paste(round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][1,1]), 2),"/",round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][2,1]), 2), "/", round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][3,1]), 2), sep = "")} else {NA},
      "RoBMA [95% CI]" = if(!is.null(maResultsObject[[6]][["Robust BMA"]])){paste(round(maResultsObject[[6]][["Robust BMA"]]$Mean[1], 2), " [", round(maResultsObject[[6]][["Robust BMA"]]$`0.025`[1], 2), ", ", round(maResultsObject[[6]][["Robust BMA"]]$`0.975`[1], 2), "]", sep = "")} else {NA},
      "PET-PEESE est [95% CI]" = paste(round(maResultsObject[[1]][["PET-PEESE"]][1], 2), " [", round(maResultsObject[[1]][["PET-PEESE"]][5], 2), ", ", round(maResultsObject[[1]][["PET-PEESE"]][6], 2), "]", sep = ""),
      "PET-PEESE" = round(maResultsObject[[1]][["PET-PEESE"]][4], 3),
      "Power to detect small SESOI" = maResultsObject[[7]][[1]],
      "Power to detect medium SESOI" = maResultsObject[[7]][[2]]
    ))
  } else {
    noquote(c(
      "k" = as.numeric(maResultsObject[[1]]$k.all),
      "ES [95% CI]" = paste(round(as.numeric(maResultsObject[[2]]$test$beta), 2), " [", round(maResultsObject[[2]]$CIs$CI_L, 2), ", ", round(maResultsObject[[2]]$CIs$CI_U, 2), "]", sep = ""),
      "95% PI [LB, UB]" = paste("[", round(maResultsObject$`Prediction interval`[1], 2), ", ", round(maResultsObject$`Prediction interval`[2], 2), "]", sep = ""),
      "SE" = round(maResultsObject[[2]]$test$SE, 2),
      round(maResultsObject[[4]]["Tau"], 2),
      "I^2" = paste(round(maResultsObject[[4]]["I^2"], 0), "%", sep = "")
    ))
  }
}

# Return format function
# Code adapted from Carter, E. C., Schönbrodt, F. D., Hilgard, J., & Gervais, W. (2018). Correcting for bias in psychology: A comparison of meta-analytic methods. Retrieved from https://osf.io/rf3ys/.
# https://github.com/nicebread/meta-showdown/blob/master/MA-methods/7-Selection%20Models.R

# Return a result data frame either in wide or long format
returnRes <- function(res, long = TRUE, reduce = TRUE) {
  if (is.null(res)) return(NULL)
  # convert all factor columns to characters
  res %>% mutate_if(is.factor, as.character) -> res
  if (long == FALSE) {
    # return wide format
    return(res)
  } else {
    # transform to long format
    longRes <- melt(res, id.vars=c("method", "term"))
    if (reduce==TRUE & nrow(res) > 1) {longRes <- longRes %>% filter(!is.na(value)) %>% arrange(method, term, variable)}
    return(longRes)
  }
}

# t-test for summary statistics -------------------------------------------

tTestSummary <- function(mean1, mean2, sd1, sd2, n1, n2, withinSS = FALSE)
{
  if(withinSS == FALSE){
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 + n2 - 2)) 
    df <- n1 + n2 - 2
    t <- (mean1 - mean2)/se 
    out <- c(mean1 - mean2, se, t, 2*pt(-abs(t),df))    
    names(out) <- c("Difference in means", "SE", "t-statistic", "p-value")
    return(out)  
  } else if(withinSS == TRUE){
    se <- sqrt((sd1^2 + sd2^2 - 2*sd1*sd2*rmCor)/(n1 - 1)) 
    df <- n1 - 1
    t <- (mean1 - mean2)/se 
    out <- c(mean1 - mean2, se, t, 2*pt(-abs(t), df))    
    names(out) <- c("Difference in means", "SE", "t-statistic", "p-value")
    return(out)
  }
}

# Random selection of effects ---------------------------------------------

# Choose effects from a single study by random (for the purpose of permutation-based methods)
duplicated.random = function(x, incomparables = FALSE, ...)
{
  if (is.vector(x))
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else (is.matrix(x))
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
}