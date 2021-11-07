# Should plots be generated?
plots <- FALSE

# Missing data ------------------------------------------------------------
#'# Percentage of missing data overall
if(plots == TRUE){
  dat %>% missmap(rank.order = TRUE, margins = c(10, 0), legend = F)
}

# Outlier diagnostics -----------------------------------------------------
tab <- sort(table(dat$correlate))[sort(table(dat$correlate)) > kThreshold]
corVect <- names(tab)

rmaObjectsDiag <- diagResults <- infResults <- vector(mode = "list", length(corVect))
names(rmaObjectsDiag) <- names(diagResults) <- names(infResults) <- corVect
for(i in 1:length(corVect)){
  xLabel <- eval(substitute(corVect[i]))
  # Initial outlier diagnostics
  # Univariate MA
  rmaUni <- dat %>% filter(correlate == corVect[i] & !is.na(yi)) %$% rma(yi = yi, vi = vi, method = "REML", slab = result)
  # MA diagnostics
  if(plots == TRUE){
    baujat(rmaUni, symbol = "slab", xlab = xLabel)
    diagResults[[i]][["Baujat"]] <- recordPlot() 
  }
  # Influence diagnostics
  infl <- influence(rmaUni)
  diagResults[[i]][["Influence"]] <- infResults[[i]] <- as.data.frame(infl$inf) %>% rownames_to_column() %>% filter(rstudent > 2.58) %>% head(.,5)
  ### Plot the influence diagnostics
  if(plots == TRUE){
    plot(infl, slab.style = 2)
    title(xLabel, line = 3.2)
    diagResults[[i]][["Influence plot"]] <- recordPlot()
  }
  rmaObjectsDiag[[i]] <- rmaUni
}

#+eval = FALSE
#fit FE model to all possible subsets
# gosh.plot <- gosh(ma.uni, progbar = TRUE, subsets = 1000, parallel = "multicore")
# plot(gosh.plot, out = , breaks = 50) # Testing the influence of single outliers

#+eval = TRUE
# Outlier removal in case of a need
# Excluding improbably big effect sizes or ES with improbably small SE, i.e. excerting a big influence on the MA model due to combination of huge ES and small variance.
# Sensitivity analysis with the outlying ESs included will be reported as well.
# dat[c(),] <- NA
