# Please load the packages required by this script by sourcing functions.R in the mainScript.

# Read in the data
dat <- read_delim("igdMetaData.csv", ";", trim_ws = TRUE)

#str(dat)
#view(dat)
dat <- dat %>% modify_at(., .at = c("pubYear", "meanIGD", "meanControl", "sdIGD", "sdControl", "F", "t", "r", "r_s"), .f = ~as.numeric(as.character(.)))
#grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", dat$meanIGD)

# Some data wrangling to get the right type of data (formatting of the raw dataset in Excel introduces a lot of junk otherwise)
dat$pReported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$pReported)))

# Compute df2 from nTotal if not reported, gender ratio (% of female), create result ID, and initialize new variables
dat <- dat %>% mutate(df2 = ifelse(is.na(df2) & !is.na(nTotal) & is.na(chiSq), nTotal - 2, df2),
                      percFemale = nFemale/(nFemale + nMale),
                      study = paste(paperID, "/", studyID, sep = ""),
                      result = 1:nrow(.),
                      useMeta = ifelse(is.na(reasonNotUsed), 1, 0),
                      useCellN = NA,
                      yi = NA,
                      vi = NA)

# Between-ss design based on means & SDs ----------------------------------
# Specify the design, compute ni and p
dat <- escalc(measure = "SMD", m1i = meanIGD, m2i = meanControl, sd1i = sdIGD, sd2i = sdControl, n1i = nIGD, n2i = nControl, data = dat)
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(yi), "meansBtw", NA), 
                      nIGD = ifelse(finalDesign == "meansBtw" & is.na(nIGD) & !is.na(nTotal), nTotal/2, nIGD),
                      nControl = ifelse(finalDesign == "meansBtw" & is.na(nControl) & !is.na(nTotal), nTotal/2, nControl),
                      a = ((nIGD + nControl)^2)/(nIGD * nControl),
                      yiConv = yi/sqrt(yi^2 + a),
                      viConv = (a^2 * (nIGD + nControl)/(nIGD * nControl) + (yi^2)/(2 * (nIGD + nControl)))/(yi^2 + a)^3,
                      ni = ifelse(finalDesign == "meansBtw", nIGD + nControl, NA))
dat <- dat %>% rowwise %>% mutate(p = as.numeric(round(tTestSummary(meanIGD, meanControl, sdIGD, sdControl, nIGD, nControl, withinSS = FALSE)["p-value"], 5))) %>% data.frame()
dat$a <- NULL
dat$yi <- dat$vi <- NA

# Show the converted ESs
dat %>% filter(finalDesign == "meansBtw") %>% select(yiConv, viConv, p, df2, nIGD, nControl, ni)

# F-Test between with df1 == 1 ---------------------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = case_when(!is.na(F) & !is.na(df1) & !is.na(df2) & df1 == 1 ~ "F1"),
                      ni = ifelse(finalDesign == "F1", df2 + 2, ni),
                      p = ifelse(finalDesign == "F1", 1-pf(F, df1, df2), p))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$nIGD + dat$nControl) >= (dat$ni - 2) & (dat$nIGD + dat$nControl) <= (dat$ni + 2), 1, 0),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$yiConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "r")$es
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$viConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "r")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$yiConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = nIGD, grp2n = nControl, es.type = "r")$es
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$viConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = nIGD, grp2n = nControl, es.type = "r")$var

# Show the converted ESs
dat %>% filter(finalDesign == "F1") %>% select(yiConv, viConv, p, F, df2, ni, directionEffect, useCellN)

# t-tests between ---------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(t) & !is.na(df2), "tBtw", finalDesign),
                      ni = ifelse(finalDesign == "tBtw", df2 + 2, ni),
                      p = ifelse(finalDesign == "tBtw", 2*pt(abs(t), df2, lower.tail = FALSE), p))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$nIGD + dat$nControl) >= (dat$ni - 2) & (dat$nIGD + dat$nControl) <= (dat$ni + 2), 1, useCellN),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$yiConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = t, totaln = df2 + 2, es.type = "r")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$viConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = t, totaln = df2 + 2, es.type = "r")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$yiConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = t, grp1n = nIGD, grp2n = nControl, es.type = "r")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$viConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = t, grp1n = nIGD, grp2n = nControl, es.type = "r")$var

# Show the converted ESs
dat %>% filter(finalDesign == "tBtw") %>% select(yiConv, viConv, p, df2, ni, directionEffect, useCellN)

# Chi^2 between -------------------------------------------------------------

# Specify the design, compute ES, var, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse(is.na(yiConv) & !is.na(chiSq) & (!is.na(nIGD) & !is.na(nControl) | !is.na(nTotal)), "chiSqBtw", finalDesign),
                      ni = ifelse(finalDesign == "chiSqBtw", ifelse(!is.na(nIGD + nControl), nIGD + nControl, nTotal), ni),
                      p = ifelse(finalDesign == "chiSqBtw", 1 - pchisq(chiSq, 1), p))

dat[dat$finalDesign == "chiSqBtw" & !is.na(dat$finalDesign),]$yiConv <- dat %>% filter(dat$finalDesign == "chiSqBtw") %$% esc_chisq(chisq = abs(chiSq), totaln = ni, es.type = "r")$es
dat[dat$finalDesign == "chiSqBtw" & !is.na(dat$finalDesign),]$viConv <- dat %>% filter(dat$finalDesign == "chiSqBtw") %$% esc_chisq(chisq = abs(chiSq), totaln = ni, es.type = "r")$var

# Show the converted ESs
dat %>% filter(finalDesign == "chiSqBtw") %>% select(yiConv, viConv, chiSq, ni, design)


# Odds ratio --------------------------------------------------------------

# d_logOR <- logOR * (sqrt(3)/pi)
# d_logOR
# ## [1] 0.2235
# V_logOR <- .08
# V_d_logOR <- V_logOR * (3/(pi^2))
# V_d_logOR


# Correlation -------------------------------------------------------------
# Specify the design, convert rho to r, and compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse((!is.na(r) | !is.na(r_s)) & !is.na(df2), "cor", finalDesign),
                      r = ifelse(is.na(r) & !is.na(r_s), abs(2*sin(r_s*pi/6)), abs(r)),
                      ni = ifelse(finalDesign == "cor", df2 + 2, ni),
                      p = ifelse(finalDesign == "cor", 2*pt(abs(r*sqrt(df2 / (1 - r^2))), df2, lower.tail = FALSE), p))

# Compute the r effect size
dat <- escalc(measure = "COR", ri = r, ni = df2 + 2, data = dat)
dat <- dat %>% mutate(yi = yi*directionEffect)

# Show the converted ESs
dat %>% filter(finalDesign == "cor") %>% select(yi, vi, r, r_s, directionEffect, ni, p)

# # # PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
# dat <- dat %>% mutate(gConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * dReported, gConv),
#                       gVarConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * ((reportedOverallN)/(reportedOverallN/2 * reportedOverallN/2) + (dReported^2)/(2 * (reportedOverallN))), gVarConv)) 

# Variable computations
# # Multiply the ES by -1 if not in the opposite direction
dat <- dat %>% mutate(yi = ifelse(is.na(yi) & (!is.na(yiConv) & !is.na(directionEffect)), directionEffect * abs(yiConv), yi),
                      vi = ifelse(is.na(vi) & (!is.na(viConv) & !is.na(directionEffect)), viConv, vi),
                      p = ifelse(p < 1e-5, 1e-5, ifelse(p > 1 - 1e-5, 1 - 1e-5, p)),
                      label = paste(paperID, "/", studyID, "/", effectID, sep = ""),
                      ni = ifelse(is.na(ni), ifelse(is.na(nIGD + nControl), df2 + 2, nIGD + nControl), ni),
                      yi_z =  .5 * log((1 + yi)/(1 - yi)),
                      vi_z = 1/(ni - 3),
                      nTerm = 2/ni,
                      ux = sqrt((1-rxxV1reference)/(1-rxxV1sample)), 
                      uy = sqrt((1-rxxV2reference)/(1-rxxV2sample)))

dat <- as.data.frame(dat)

# dat %>% filter((!is.na(meanIGD) | !is.na(meanControl) | !is.na(sdIGD) | !is.na(sdControl) | !is.na(nIGD) | !is.na(nControl) | !is.na(F) |!is.na(t) | !is.na(r) | !is.na(chiSq) | !is.na(otherES)| !is.na(df1) | !is.na(df2) | !is.na(nTotal) | !is.na(r_s)) & (is.na(yi) | is.na(vi)) & !is.na(directionEffect) & is.na(otherES)) %>%
#    select(label, meanIGD,meanControl,sdIGD,sdControl,nIGD,nControl,F,t,r, r_s, chiSq,df1,df2,nTotal,r_s, directionEffect, yi, vi, gdCriteria, correlate) %>% view()

# dat %>% filter(!is.na(yi) & !is.na(vi) & (yi > .85 | yi < -.85)) %>%
#    select(label, meanIGD,meanControl,sdIGD,sdControl,nIGD,nControl,F,t,r, r_s, chiSq,df1,df2,nTotal,r_s, directionEffect, yi, vi, gdCriteria, correlate) %>% view()

# Remove outliers (based on the results from the maDiag script)
# dat <- dat %>% filter(!result %in% c())

# dat %>% select(label, yi, vi, yiConv, dReported, directionEffect, meanIGD, meanControl, sdExp, sdCtrl, seExp, seCtrl, nIGD, nControl, df2, F, t, r) %>% view()

# dat %>% filter(ni != nTotal) %>% select(ni, nTotal, df2, nIGD, nControl)
# dat %>% filter(abs(dReported - gConv) > .2) %>% select(result, gConv, dReported)
# 
# dat$studentSample <- ifelse(dat$sampleType == "student", 1, ifelse(dat$sampleType == "general", 0, NA))

# Remove outliers (based on the results from the maDiag script)
# dat <- dat %>% filter(!result %in% c(194))

