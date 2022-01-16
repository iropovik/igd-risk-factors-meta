# Risk and protective factors for (internet) gaming disorder: A meta-analysis

This GitHub repository contains Materials, R code, Data and Appendices for the paper "Risk and protective factors for (internet) gaming disorder: A meta-analysis".

{{TOC}}

# Appendix A

***
## Method
The protocol of this meta-analysis has been preregistered at PROSPERO (ID:[CRD42020187776](https://www.crd.york.ac.uk/prospero/display_record.php?RecordID=187776)), with a list of deviations from the protocol disclosed at ??. Data, R code, and analytical outputs can be found at ??.

### Search strategy
Two different strategies were used to identify relevant studies. First, after piloting the search string, we searched the following databases: Scopus, Web Of Science Core Collection, PsycInfo (via ??), PubMed (via ??), OSF Preprints, ThesisCommons, ProQuest Dissertations and Theses. The latter three databases indexed grey literature. The search was carried out on DATE?? and was limited to empirical studies written in English and published from 2013. Second, we carried out a forward citation search in Google Scholar by screening studies that cited 13 most widely used measures of GD (referenced in Appendix A) that complied with the DSM-5 or ICD-11 definition of IGD.

### Search string
Given the formal operationalization of the main domain of our interest, following the DSM-5 or ICD-11 definitions, the term “gaming disorder” was made a necessary constituent of the search string (we are unaware of any relevant synonyms of “gaming disorder”; neither “gaming disorder” nor “gaming” have any matching MeSH terms in the PubMED database), along with at least synonym related to concepts of risk or protective factors. Translations of the Boolean search string for all databases can be found in Appendix ??. As an example, the search string for Scopus was as follows: 
	 	 	 	

> **Scopus Boolean string**: TITLE-ABS-KEY("gaming disorder" AND ("risk factor" OR "factor" OR "risk" OR pron* OR predict* OR correlat* OR "role of" OR link* OR associat* OR relat* OR meta* OR cause* OR trigger* OR promot* OR vulnerab* OR protect* OR buffer* OR resilience* OR prevent* OR reduc*)) AND LANGUAGE(english) AND DOCTYPE(ar OR re) AND PUBYEAR AFT 2012

## Inclusion criteria
Conceptually, the present meta-analysis aimed to include all available studies carried out on all populations that measured and reported zero-order associations between GD and any specific risk or protective factors.

### Population
As the diagnostic criteria for gaming disorder apply equally to all gender and age groups (Petry et al., 2014 and empirically supported by e.g., Sigerson et al., 2017), we applied no exclusion criteria regarding the target population. The population of interest thus included general population, clinical sample, or groups of specific gamers (see e.g., King et al., 2020).

### Predictors/Exposures
The study had to assess at least one correlate of GD and report some type of a zero-order relationship with GD. The exposure in this meta-analysis was represented by various risk and protective factors of gaming disorder. King and Delfabbro (2019) categorized these factors into three main groups:
1)	individual factors: demographic factors such as age, gender, psychological or psychopatological factors (e.g., impulsivity, neuroticism, introversion, narcissism, aggressiveness, trait anxiety, deficient self-regulation, low self-esteem and self-efficacy), comorbidity, and family environment
2)	external factors: peer or familial influences, frequency or duration of gaming, gaming environment, and relational trauma
3)	gaming-related factors such as specific types of games and game features
 
### Outcome
The study had to assess gaming disorder, operationalized as a presence or an absence of the disorder (i.e., a binary format typical for the DSM-5 and ICD-11; the disorder is present if the subject endorse five or more criteria or have a summary score higher than the predetermined cut-off, mostly based on the latent class analysis) or as a degree of severity of gaming disorder symptoms (continuous scale, typical for research settings). As implied, we excluded studies where the definition and measurement of GD or IGD did not comply with the DSM-5 or ICD-11 (effectively most of the pre-2013 literature on gaming disorder). 

More specifically, the present meta-analysis was focused on Gaming Disorder as characterized in:
1. The recent 11th revision of the International Classification of Diseases (ICD-11) “by a pattern of persistent or recurrent gaming behavior (‘digital gaming’ or ‘video-gaming’), which may be online (i.e., over the internet) or offline, manifested by: (1) impaired control over gaming (e.g., onset, frequency, intensity, duration, termination, context); (2) increasing priority given to gaming to the extent that gaming takes precedence over other life interests and daily activities; and (3) continuation or escalation of gaming despite the occurrence of negative consequences. The behavior pattern is of sufficient severity to result in significant impairment in personal, family, social, educational, occupational or other important areas of functioning.” 
2. The fifth edition of the Diagnostic and Statistical Manual of Mental Disorders (DSM-5; under the name Internet Gaming Disorder) as “Persistent and recurrent use of the Internet to engage in games, often with other players, leading to clinically significant impairment or distress as indicated by five (or more) of the following in a 12-month period: (1) Preoccupation with Internet games. (2) Withdrawal symptoms when Internet gaming is taken away. (3) Tolerance— the need to spend increasing amounts of time engaged in Internet games. (4) Unsuccessful attempts to control the participation in Internet games. (5) Loss of interests in previous hobbies and entertainment as a result of, and with the exception of, Internet games. (6) Continued excessive use of Internet games despite knowledge of psychosocial problems. (7) Has deceived family members, therapists, or others regarding the amount of Internet gaming. (8) Use of Internet games to escape or relieve a negative mood (e.g., feelings of helplessness, guilt, anxiety). (9) Has jeopardized or lost a significant relationship, job, or educational or career opportunity because of participation in Internet games.” 

Before the inclusion of the Gaming Disorder into the DSM-5, many authors used different nomenclatures for problematic gaming such as video game addiction, computer game playing dependence, internet addiction disorder, video game dependency, problematic online gaming, or pathological video-game use (Pontes et al., 2014). The logical consequence of the inconsistent conceptualization was inconsistency in measurement tools as there were no standard diagnostic criteria (King & Delfabbro, 2014a; Petry & O’Brien, 2013; Pontes & Griffiths, 2014). This led to the development of different measures assessing different criteria using different operationalizations, often with no justification for some of the criteria at all (King & Delfabbro, 2014a, 2018a; King et al., 2013). Pontes and Griffiths (2014) called that a diagnostic confusion that probably led to underestimation in prevalence rates. As many prominent authors call for a unified approach in the study of problematic gaming (e.g., Griffiths et al., 2014), and the existence of standard diagnostic criteria is a prerequisite for a valid diagnosis, we have decided to limit our search only to publications that adopted DSM-5 or ICD-11 nomenclature (i.e., post-2013). Therefore, we used the established disorder name as defined in ICD-11 and DSM-5 as the search term. 
 
## Study selection
We aimed to include all available studies reporting a measure of association between gaming disorder and at least one risk factor. The set of records identified in the search phase was checked for duplicates and retractions. Screening of all the records by titles, abstracts, and keywords was carried out in Rayyan (Ouzzani et al., 2016) by two independent coders blinded to each other’s decisions. Full texts of records that were not marked for exclusion, were further screened for eligibility. None of the studies was excluded unless agreed upon by both coders (an additional coder was consulted in case of disagreement). Studies for which we were unable to obtain a full text either by database search or direct request from the corresponding author were excluded. The detailed flow of study selection is outlined in the Prisma Flowchart (Page et al., 2021; see Figure 1 in the main manuscript).

## Data extraction
Following the screening and selection of the studies, the data were extracted by two coders. After coding the random first 10% of all included studies (but not less than 30 studies), we checked for the inter-rater agreement using Cohen’s Kappa, subsetting only variables that entered some of the meta-analytic or meta-regression models (i.e., effect sizes, test statistics, data on moderators, aspects of the study design). If the inter-rater agreement for any of the coded variables was not substantial (κ < .6 for metric variables, % agreement < 80% for discrete variables), we discussed and revised the coding procedure for the given variable, did an additional check of Cohen’s κ and agreement, and possibly revisited the coding scheme at 20% of coded articles (but not less than 60 studies). All the disagreements were resolved by discussion and, if needed, by consulting a third coder. After checking for inter-rater agreement, the data for each study were coded by one of the two coders. Each coder kept track of any ambiguities or unresolved issues that required either to be discussed at coding meetings or contacting the original authors.

From each paper, we extracted all relevant effect sizes. We coded the following data for each coded effect: 
1. bibliographic information
2. study meta-data (H5 index of the journal, citation count for the paper, effect identifiers indicating hierarchical structure of the data, status of being published, study design, centrality of the given effect for the paper, author’s use of predictor selection based on statistical significance)
3. outcome and risk factor-related variables (gaming disorder measure used, way of operationalizing the gaming disorder, type of risk/protective factor)
4. gaming-related variables (genre of the game, platform)
5. sample characteristics (type of sample, mean age, gender proportion, restrictions in sampling)
6. effect data (various types of reported effect sizes and test statistics, sample sizes, degrees of freedom, reported p-value, direction of the effect, sample reliabilities for both IGD and correlate, and reliabilities or σ2 in the reference population for both variables). 

In case of missing data required to recompute the effect size, we contacted the corresponding author of the study. In case of no response in two weeks’ time, we followed up with a reminder. If, ultimately, we were unable to recover or soundly infer the missing data, the given effect was not included in the relevant meta-analytic models. Full coding scheme is available at: [??](osf.io).

### General coding principles
The following general coding principles were followed:
1. Authors’ claim that IGD was measure in accordance with DSM-5 or ICD-11 was not taken at face value. Instead, we always checked whether the actual measure that was used adhered to DSM-5 or ICD-11. If needed, this information was gathered from external sources, e.g., relevant validation studies. If IGD was assessed only based on a clinical judgment, we gave the authors the benefit of the doubt and included the effect if the authors stated that they followed DSM-5 or IGD-11.
2. Effects based on unordered multinominal predictors, omnibus effects in general, or estimates from path or SEM models (if the predicted IGD variable was endogenous) were excluded, only zero-order effects were included.
3. We always aimed to favor effects based on the maximum amount of information in the variables, i.e., we favored effects based on continuous rather than discrete or dichotomized variables. So for instance, when the authors reported conceptually the same effect using a correlation and a t-test based on a median split, we coded the correlation.
4. After coding all relevant effects, we made sure that none of the included variables was a constituent or a product of other variables reported in any given study.
5. If a variable was employed for participant matching or to restrict the variability in the sample in any way, its relationship with IGD was not coded.
6. For group designs, the comparator group had to be a control group, i.e., either a relevant general population or healthy controls, not a group sampled from a population having a different disorder (e.g., gambling disorder).
7. For longitudinal studies, we extracted the baseline data. If the predictor was an experimentally manipulated variable, we coded only the effect estimate prior to the treatment.

### Variable-specific coding rules:
1. **focal**: If the given effect was mentioned in the abstract, it was coded as focal. If the authors reported only the significant predictors of IGD in the abstract, all the tested predictors were equally deemed focal as they probably would have been reported in the abstract had they been significant.
2. **sampleRestricted**: Coding 1 if at least one of the groups (e.g., in group designs) is restricted by design of the study.
3. **meanAge**: If mean age only for the entire sample was reported, we took that as an approximation of mean age even if the effect was based only on a subset of participants. If the authors reported mean age only for subgroups and the effect was estimated on the full sample, we took a N-weighted average. If the authors reported only age range, we took the central value of that range as an approximation.
4. **sdIGD, sdControl**: If the authors reported only SEs, we computed the SDs.
5. **df2**: If the df for a chi-square test was not reported, we computed it based on the dimensions of the contingency table. In any other cases, df2 was not computed and only coded, if available.
6. **rxxV1sample, rxxV2sample**: In (rare) cases when the authors reported the sample reliabilities for the subscales only as a range, we took the middle of the range as the estimate of reliability.

### Effect size recomputation
If the included effect was not reported in the zero-order Pearson’s _r_ metric, we converted from the reported test statistics (_t_, _F_, χ^2), different effect size (Spearman’s _rho_) or the reported group means, _SDs_, and _Ns_. The computation and conversion of all effect sizes was carried out in code, using formulas laid out in Borenstein et al. (2009). For all group designs, we checked for an approximate concord between the group _Ns_ and the total sample size (_N_ +/-2). If these matched, we used the respective group _Ns_. Otherwise, we used the reported degrees of freedom to approximate the actual _N_, while assuming a balanced design. Also, if only the total sample size was reported, we also assumed a balanced design. Lastly, due to an arbitrary nature of the reported effect size directions (dependent upon scaling of the respective variables), we took absolute values and multiplied by the coded direction of the effect (-1 for a substantively negative relationship with GD, 1 for a positive relationship with GD).

## Analysis

### Strategy for data synthesis
The available quantitative evidence was synthesized using multilevel random-effects models utilizing the maximum-likelihood estimator. For any identified risk factor, we chose to refrain from doing any quantitative inference unless the given analytic model was based on more than 10 studies. We included all the relevant effects from every study. This introduced two types of dependencies among the effects – nesting of effects within papers and clustering due to estimation of effects based on the same sample. We modelled these dependencies using the robust sandwich-type variance estimation (RVE) with the “Correlated and hierarchical effects” working model (Pustejovsky & Tipton, 2020)[^This allowed us to simultaneously account for both types of dependencies among the effects, namely due to nesting of effects within studies and estimation of the effects based on the same participants. Because the data on sampling correlations among the effects tend to be unavailable, we assumed a constant sampling correlation of .5.]. To test for equality of effect sizes across the levels of the moderators studied, we used the robust HTZ-type Wald test (Pustejovsky & Tipton, 2020). The p-values for individual contrasts were adjusted using Holm’s method.

For each meta-analytic model, we estimated the mean weighted effect size, its 95% confidence intervals and prediction intervals, relative and absolute heterogeneity (_I^2_ and tau, with separate estimates of between- and within-cluster heterogeneity for the latter), and computed the proportion of significant results.

Prior to our analyses, we carried out an in-depth diagnosis of the random-effects meta-analytic model for each separate correlate. Specifically, we screened for influential outliers using the Baujat plot and influence diagnostics indices. Outliers exerting an excessive influence on the meta-analytic model (yielding a standardized residual > 2.58) were then excluded in a sensitivity analysis.

All the data, effect size conversion script, and the data analytic R script are freely available at the project’s Open Science Framework page [https://osf.io/9mzgr/](https://osf.io/9mzgr/). All models were fitted using restricted maximum-likelihood estimation using R packages _metafor_, version 2.5 (Viechtbauer, 2010) and _clubSandwich_, version 0.4.2. (Pustejovsky, 2020). The data analysis was carried out in R also using the following packages: _esc_ (Lüdecke, 2017), _tidyverse_ (Wickham et al., 2019), _lme4_ (Bates, Maechler, Bolker, Walker, 2015), _dmetar_ (Harrer et al., 2019), and _psych_ (Revelle, 2018).

### Moderation analyses
For each of the identified risk factors, we studied the following moderators of GD.   
 
1.	**Gender**. As shown by King and Delfabbro (2019) GD prevalence favors males: males play different and more riskier types of games (eSports, MMORPG, FPS), males are more at risk due to personality issues, attention problems, risk-taking, impulsivity and sensation-seeking, and lack help-seeking tendencies. 
2.	**Age**. According to King and Delfabbro (2019) adolescence is the most vulnerable developmental period. It is hypothesized that gaming can serve as a way of coping with developmental changes and problems adolescents have to face, which creates a potential for its excessive use. We will use two different time frames for adolescence (more traditional: <18 years/ ≥18 years (King & Delfabbro, 2019) and the recent one: <24 years/ ≥24 years (Sawyer et al., 2018)).
3. 	**Gaming style**. Whereas some authors consider competitive style of play (eSports) as riskier, others do not (see Chung et al., 2019; King & Delfabbro, 2019; Faust et al., 2013; Griffiths, 2009). We will differentiate between competitive and non-competitive styles of gaming. Either one of the two criteria has to be met for the subgroup to be considered as competitive gamers: a) the sample explicitly consists  of esports players, competitive video game players, electronic/virtual/digital sports/competition players, or professional players (see Reitman et al., 2019) or b) the sample implicitly consists of players which are routinely involved in training, leagues, ladders or tournaments both on the individual or team level, which correspond to the very definition of eSports (Hamari & Sjöblom, 2017).
4. 	**Game genre**. Most of the games are classified either as MMO, MMORPGs, RPGs, FPS or shooter games, RTS, sport, action, battle royale, simulator fighting, adventure, puzzle, racing, arcade, or card and board games. 
5. 	**Platform**. Different types of platforms, such as PC, console, mobile phones or tablets are associated with different accessibility of gaming.
6. 	**Sample type**. Coded either as online gamers, offline gamers, clinical sample, or non-gamers. For example, prevalence rates depend on the sample type, i.e. whether gamers or non-gamers were studied (Mihara & Higuchi, 2017).
7. 	**Operationalization of the gaming disorder**. Whether GD/IGD was measured with questionnaire and treated as interval variable (a degree of problematic gaming), a binary variable associated with some cut-off (DSM5 ≥ 5 criteria or a latent class analysis) or assessed by a clinician.

## Adjustment for bias
In the present meta-analysis, we attempted to adjust the included effects for (1) psychometric artifacts like measurement error, range restriction/enhancement, and collider bias, (2) publication bias, as well as (3) some other more specific forms of biases in the publication process.

### Psychometric artifacts
First, we carried out a correction for various psychometric artifacts that tend to bias individual effect sizes. In this respect, we employed the meta-analytic psychometric method by Hunter & Schmidt (2014), adjusting for the attenuation of the studied associations due to measurement error (unreliability, group misclassification). 

Apart from measurement error, we also accounted for selection effects of range restriction/enhancement and collider bias using methods described by Wiernik & Dahlke (2020). Here, we used the artifact-distribution method (Taylor series approximation), as implemented in the _psychmeta_ package (Dahlke & Wiernik, 2019). Specifically, we carried out adjustments for bivariate direct range restriction, coupled with a correction for measurement error in IGD (and not also in the covariate, as we were interested in the predictive power of the observed correlates). We also report an analysis, where we (1) adjusted also for the measurement error in the covariate (to examine the relationships at the construct level) and (2) adjusted for the indirect bivariate selection (the effect of conditioning on a collider). If the required sample and reference reliabilities or variances were unavailable for the given effect, no correction was carried out, effectively assuming reliability of 1 for both variables that the effect was based on.

Lastly, as a secondary approach to examine range restriction, we also estimated whether the population variances tended to be smaller for restricted samples (as judged by the raters), possibly leading to more easily detectable effect sizes, i.e., larger standardized effects for designs applied to restricted populations.

### Publication bias
More importantly, we also tried to adjust the meta-analytic estimates for publication bias. We used selection models and regression-based models, supplemented by several exploratory sensitivity analyses. Although bias adjustment methods assume a more realistic selection process, they may fail to recover the “true” magnitude of the studied effects under a number of realistic conditions. The estimates should thus rather be seen as approximations (see Ropovik et al., 2021). If the adjusted estimates from selection models markedly diverged from the crude meta-analytic estimates, then we primarily used bias-corrected estimates to guide our substantive inferences.

#### Selection models
Primarily, we used a permutation-based implementation of a 4-parameter selection model (see McShane et al., 2016). Selection models are a statistically principled family of models that can flexibly map the direct functional form of the biasing selection for statistical significance. They model models both, how the data were generated in the absence of any bias (population effect size and heterogeneity) as well as the selection mimicking publication process (based on the p-value). It employs the maximum likelihood estimation, trying to find the values of these multiple parameters under which the observed data are most likely, relative to any other set of possible parameter values (Iyengar & Greenhouse, 1988; McShane, Böckenholt, & Hansen, 2016).

In case there was too little data (less than 4 effects per a p-value interval), the estimation procedure automatically reverted to the 3-parameter selection model, leaving out the parameter representing the likelihood of a result being in the opposite direction. All selection models  analyzed only subset of effect deemed to be the study’s focal effects (reported in the abstract). The dependencies among the effects in selection models were handled using a permutation approach. For each set of effects, only a single focal outcome per study was randomly drawn, the model was repeatedly estimated in 500 iterations and averaged over this set of iterations by selecting the model with the median estimate. In Appendix ??, we also report robustness analyses, where we examined the sensitivity of adjusted effect size estimates to different assumptions about the biasing process. Specifically, we computed a series of Vevea and Woods (2005) 10-step function models with a priori defined selection weights, reflecting three different levels of severity of bias. Moreover, we employed the Bayesian model-averaging approach that is more robust to model misspecification, as it averages over a range of plausible selection models and the regression-based PET-PEESE model (Bartoš & Maier, 2020). The exact specification of the models can be seen in code at ??.

In Appendix ??, we also report the results of a one-parameter selection model, the p-uniform method (van Assen, van Aert, & Wicherts, 2015). The effect size estimation by p-uniform* is based on the fact that p-values follow a uniform distribution at the true effect size. To  ease the assumption of included effects being independent, we also implemented the p-uniform method using a permutation-based procedure (with 500 iterations). 

#### Regression-based models
As a secondary publication bias adjustment method, we used the PET-PEESE method (Stanley & Doucouliagos, 2014). PET-PEESE is a weighted conditional regression-based model regressing the effect size on an estimate of its precision. It is based on the notion that publication bias can be inferred from a pattern when less precise studies tend to show consistently stronger effect sizes on average. Namely, effects in small-_N_ studies need to be larger to outweigh the correspondingly large expected sampling variability and get detected. The regression slope in PET-PEESE then reflects the detectability of publication bias (possibly mixed with other small-study effects), while the intercept represents an estimate of the average effect size for a hypothetical, infinitely precise study (Stanley & Doucouliagos, 2014). In line with our previous work (see IJzerman et al., 2021; Kolek et al., 2022), we used a multi-level RVE-based implementation of PET-PEESE, with the same structure of random-effects as in the unadjusted models. Apart from that, we used √(2/_N_) as the predictor term in PET and a 2/_N_ term in PEESE, instead of standard error and variance as the proxies for precision. Unlike when using variance (which is calculated using the effect size), the use of _N_-based predictors does not induce an artificial correlation with the effect size, which is conceptually the estimand in these models. Models using _N_-based predictors therefore exhibit a markedly lower false-positive rate (Pustejovsky, 2017). As originally proposed, the employed models used PET as the conditional estimator for PET-PEESE. In case PET estimate was significant, PEESE was used, if PET was not significant, the model returned the PET estimate as the bias-corrected estimate. As a sensitivity analysis, we also conservatively employed the 3-parameter selection model as the conditional estimator for PET-PEESE instead of PET. The reason is that the former tends to perform better in terms of precision and error rates in similar analytic scenarios (Carter et al. 2019). 

In Appendix ??, we also report the results of another regression-based model, the Weighted Average of the Adequately Powered studies (WAAP-WLS; Stanley et al., 2016). The combined WAAP-WLS estimator tries to identify studies that are adequately powered to detect an average effect. If there are fewer than two such studies, the method falls back to the WLS estimator (Stanley, Doucouliagos, & Ioannidis, 2017). If there are at least two adequately powered studies, WAAP returns a WLS estimate based on effects from only those, adequately powered studies.

### Other forms of bias
Apart from examining the impact of psychometric artifacts and publication bias, we also explored other more specific biasing effects of the publication process. Specifically, we examined (1) the citation bias, testing whether more highly-cited studies report larger effect sizes; (2) decline bias, testing whether studies showing more extreme (possibly opposite) results appear early in the research line rather than late, as data accumulate (see Fanelli, Costas, & Ioannidis, 2017). Here, we employed the absolute deviation from the meta-analytic effect size as the outcome variable in a linear mixed-effects model, controlling for increasing precision over time (decreasing SE). Lastly, (3) we examined bias arising from the presence of conflict of interests, employing a subgroup analysis. Results of all these analyses and their interpretations can be found in Appendix ??.

## Quality of evidence assessment
Apart from answering substantive questions, we also attempted to assess the quality and integrity of the available evidence.

### Assessment of evidential value
First, we estimated the evidential value of the available studies using a permutation-based p-curve analysis (Simonsohn, Nelson, & Simmons, 2014), subsetting only one effect per study that was coded as focal (the result was reported in abstract in some way) for each iteration. This analysis aimed to provide some insight into the degree of selective reporting in the literature on different risk factors. In case there is evidential value in the given literature, then regardless of power, a right-skewed distribution of p-values can be expected. On the other hand, a left-skewed distribution of p-values may indicate a substantial prevalence of questionable research practices in the literature. 
We employed p-values recomputed based on reported test statistics. The within-study clustering of p-values was addressed employing a permutation-based procedure, iteratively sampling only a single focal p-value from each study (with 200 iterations), estimating the p-curve, and averaging over the set by selecting the model with the median z-score for the right-skew of the full p-distribution.

### Numerical inconsistencies in reported p-values
Second, all included papers were screened using a machine-based procedure for the presence of inconsistencies in reported p-values, employing the statcheck package (Epskamp & Nuijten, 2018). The method works as follows. First, pdf files are converted to plain text files. These are scanned in full for statistical results reported in APA style. Next, test statistics and degrees of freedom are extracted to recompute the p-values. Lastly, the recomputed p-values are compared to the reported p-values. Based on these extracted data, we computed how often were the reported p-values inconsistent with the reported test statistics and how many of those cases led to an error in inference.

### Median statistical power in the literature
Third, we computed the median statistical power to detect various smallest effect sizes of interest (_r_ = 0.10, 0.30, and 0.50), as well as median power to detect the bias-corrected 4PSM and PET-PEESE estimates.

## Additional sensitivity analyses
In Appendix ??, we also report the results of the following, more general sensitivity analyses. (1) We used Fisher’s _z_-transformed effect sizes instead of Pearson’s _r_. (2) We excluded all effect sizes coming from studies employing a selection inference approach, i.e., disregarding non-significant effects using, e.g., stepwise selection procedures or regularization techniques. Namely, these effect sizes are inflated, on average, leading to overestimation of the mean effect sizes.

# Appendix B

# Appendix C
