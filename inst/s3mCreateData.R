## *** Probabilistic Sensitivity Analysis
## Fitted Parameters
## 44% (95% CI 35–54) benign
## 17% (95% CI 7–26) G6
## Utilities
## MRR ERSPC
## *** One-way
## 44% (95% CI 35–54) benign
## 17% (95% CI 7–26) G6
## Utilities
## MRR ERSPC
## Discount rates
## | Utils | Cost |
## |-------+------|
## |     0 |    0 |
## |     3 |    3 |
## |     5 |    5 |
## |     0 |    3 |

## *** Special investigation for appendix
## Cost of panel
## 2-yearly

library(prostata)
library(dplyr)
library(parallel)

## path <- "../Data"
path <- "./"
makeModels <- TRUE
run_name <- "20181020-bx-staggered"

cost_parameters <- prostata:::FhcrcParameters$cost_parameters
## change the values!

makeModel <- function(screen = "introduced_screening_only",
                      n = 1e7, #1e7 is 3*20min on 2cores CHECK this before running!!!
                      start_screening=55,
                      stop_screening=70,
                      screening_interval=4,
                      introduction_year = 1997,
                      popDefault = introduction_year - start_screening, # single cohort 1942
                      discountRate.costs = 0.03,
                      discountRate.effectiveness=0.03,
                      formal_compliance=1.0, # ERSPC compliance (not observed in STHLM)
                      formal_costs=0.0, # Include GP cost for taking the test
                      utility_estimates=prostata:::FhcrcParameters$utility_estimates,
                      cost_parameters=cost_parameters,
                      production = prostata:::FhcrcParameters$production,
                      grade.clinical.rate.high = prostata:::FhcrcParameters$grade.clinical.rate.high,
                      utility_scale = prostata:::FhcrcParameters$utility_scale,
                      PSA_FP_threshold_nCa = prostata:::FhcrcParameters$PSA_FP_threshold_nCa,
                      PSA_FP_threshold_GG6 = prostata:::FhcrcParameters$PSA_FP_threshold_GG6,
                      PSA_FP_threshold_GG7plus = prostata:::FhcrcParameters$PSA_FP_threshold_GG7plus,
                      panelReflexThreshold = prostata:::FhcrcParameters$panelReflexThreshold,
                      mc.cores= if (Sys.info()["nodename"] == "vector.meb.ki.se") {
                                    floor(parallel::detectCores() / 4) # share on vector
                                } else { # greedy for all other systems
                                    parallel::detectCores()
                                },
                      panel=FALSE,
                      preserve = NULL,
                      includeBxrecords = FALSE,
                      ...) {

    function(screen, parms = NULL, ..., pop=popDefault, panel=panel, pres = preserve) {
        newparms <- list(start_screening = start_screening,
                         stop_screening = stop_screening,
                         screening_interval = screening_interval,
                         introduction_year = introduction_year,
                         discountRate.effectiveness = discountRate.effectiveness,
                         discountRate.costs = discountRate.costs,
                         formal_compliance = formal_compliance,
                         formal_costs = formal_costs,
                         utility_estimates = utility_estimates,
                         cost_parameters = cost_parameters,
                         production = production,
                         grade.clinical.rate.high = grade.clinical.rate.high,
                         utility_scale = utility_scale,
                         PSA_FP_threshold_nCa = PSA_FP_threshold_nCa,
                         PSA_FP_threshold_GG6 = PSA_FP_threshold_GG6,
                         PSA_FP_threshold_GG7plus = PSA_FP_threshold_GG7plus,
                         panelReflexThreshold = panelReflexThreshold,
                         includeBxrecords = includeBxrecords)

        ## if (not(all(names(newparms) %in% names(parms))))
        if (!is.null(parms)) newparms <- modifyList(newparms, parms)
        if (!is.null(pres)) newparms <- modifyList(newparms, pres)

        cat("Settings: ")
        cat("n:", n)
        cat(", screen:", screen)
        cat(", mc.cores:", mc.cores)
        cat(", pop:", pop)
        cat(", parms:", unlist(newparms)) #TODO also print variable names
        cat(", panel:", panel, "\n")
        set.seed(1234)
        callFhcrc(n = n, screen = screen, mc.cores = mc.cores,
                  ## This assumes recruitment in a single calendar year across ages
                  ## Study mean age @rand 60.2 we get 62.5 @screen.
                  pop = pop, flatPop = TRUE, parms = newparms, panel = panel, ...)
    }
}

modelSet <- function(model) {
    cat("NOTE: Processing a model set.\n")
    scenarioSet <- list(modelCtrl = model("noScreening", panel = FALSE),
                        modelPSA = model("introduced_screening_only", panel = FALSE),
                        modelS3M = model("introduced_screening_only", panel = TRUE),
                        modelS3M15 = model("introduced_screening_only", panel = TRUE,
                                           parms = list(PSA_FP_threshold_nCa = 4.38,
                                                        PSA_FP_threshold_GG6 = 3.95,
                                                        PSA_FP_threshold_GG7plus = 3.11,
                                                        panelReflexThreshold = 1.5)),
                        modelS3M20 = model("introduced_screening_only", panel = TRUE,
                                           parms = list(PSA_FP_threshold_nCa = 4.56,
                                                        PSA_FP_threshold_GG6 = 4.31,
                                                        PSA_FP_threshold_GG7plus = 3.28,
                                                        panelReflexThreshold = 2.0)))
} # this has precedence i.e. all is good expect sensitivity with these parameters see preserve

## ## These would need to be updated to be used again,
## ## see prostata:::FhcrcParameters$utility_estimates
## utility_favorable <- c("Invitation" = 1.0,
##                        "Formal PSA"= 1.0,
##                        "Formal panel"= 1.0,
##                        "Opportunistic PSA"= 1.0,
##                        "Opportunistic panel"= 1.0,
##                        "Biopsy"= 0.94,
##                        "Prostatectomy part 1"= 0.90,
##                        "Prostatectomy part 2"= 0.91,
##                        "Radiation therapy part 1"= 0.91,
##                        "Radiation therapy part 2"= 0.88,
##                        "Active surveillance"= 1.0,
##                        "Palliative therapy"= 0.24, # weird order heijnsdijk 2012
##                        "Terminal illness"= 0.24, # weird order heijnsdijk 2012
##                        "Metastatic cancer"= 0.855,
##                        "Death"= 0.0)
## utility_unfavorable <- c("Invitation" = 1.0,
##                        "Formal PSA"= 0.99,
##                        "Formal panel"= 0.99,
##                        "Opportunistic PSA"= 0.99,
##                        "Opportunistic panel"= 0.99,
##                        "Biopsy"= 0.87,
##                        "Prostatectomy part 1"= 0.56,
##                        "Prostatectomy part 2"= 0.70,
##                        "Radiation therapy part 1"= 0.71,
##                        "Radiation therapy part 2"= 0.61,
##                        "Active surveillance"= 0.85,
##                        "Palliative therapy"= 0.86, # weird order heijnsdijk 2012
##                        "Terminal illness"= 0.40, # weird order heijnsdijk 2012
##                        "Metastatic cancer"= 0.845,
##                        "Death"= 0.0)
## favUtilities = modelSet(makeModel(utility_estimates=utility_favorable)),
## unFavUtilities = modelSet(makeModel(utility_estimates=utility_unfavorable)),

nm <- list(c("Formal PSA", "Formal panel", "Opportunistic PSA", "Opportunistic panel"),
           "Biopsy",
           c("Prostatectomy part 1", "Prostatectomy part 2"),
           c("Radiation therapy part 1", "Radiation therapy part 2"),
           "Active surveillance",
           "Postrecovery period",
           "Palliative therapy",
           "Terminal illness")

{models <- list(baseScenarios = modelSet(makeModel()),
                intens2Yearly = modelSet(makeModel(screening_interval = 2)),
                intens8Yearly = modelSet(makeModel(screening_interval = 8)),
                discCost0Eff0 = modelSet(makeModel(discountRate.costs = 0.0,
                                                   discountRate.effectiveness = 0.0)),
                discCost5Eff5 = modelSet(makeModel(discountRate.costs = 0.05,
                                                   discountRate.effectiveness = 0.05)),
                discCost3Eff0 = modelSet(makeModel(discountRate.costs = 0.03,
                                                   discountRate.effectiveness = 0.0)),
                ## These should be the same but scenario/case dependent
                lowBiomarker =  modelSet(makeModel(preserve = list(PSA_FP_threshold_nCa = 3.89,
                                                                 PSA_FP_threshold_GG6 = 3.3))),
                highBiomarker = modelSet(makeModel(preserve = list(PSA_FP_threshold_nCa = 4.66,
                                                                   PSA_FP_threshold_GG6 = 4.2))),
                lowBiomarker_R15 = modelSet(makeModel(preserve = list(
                                                          PSA_FP_threshold_nCa = 4.06,
                                                          PSA_FP_threshold_GG6 = 3.57,
                                                          PSA_FP_threshold_GG7plus = 2.94,
                                                          panelReflexThreshold = 1.5))),
                highBiomarker_R15 = modelSet(makeModel(preserve = list(
                                                           PSA_FP_threshold_nCa = 4.76,
                                                           PSA_FP_threshold_GG6 = 4.3,
                                                           PSA_FP_threshold_GG7plus = 3.34,
                                                           panelReflexThreshold = 1.5))),
                lowBiomarker_R20 = modelSet(makeModel(preserve = list(
                                                          PSA_FP_threshold_nCa = 4.25,
                                                          PSA_FP_threshold_GG6 = 3.91,
                                                          PSA_FP_threshold_GG7plus = 3.06,
                                                          panelReflexThreshold = 2.0))),
                highBiomarker_R20 = modelSet(makeModel(preserve = list(
                                                           PSA_FP_threshold_nCa = 4.97,
                                                           PSA_FP_threshold_GG6 = 4.62,
                                                           PSA_FP_threshold_GG7plus = 3.51,
                                                           panelReflexThreshold = 2.0))),
                unFavUtilities_m20 = modelSet(makeModel(
                    utility_estimates = sapply(prostata:::FhcrcParameters$utility_estimates,
                                               function(x) max(1 - (1-x)*1.2, 0.0)))),
                favUtilities_p20 = modelSet(makeModel(
                    utility_estimates = c(sapply(prostata:::FhcrcParameters$utility_estimates,
                                                 function(x) 1 - (1-x) * 0.8)[-16], "Death"=0))),
                lowCosts = modelSet(makeModel(
                    cost_parameters = sapply(prostata:::FhcrcParameters$cost_parameters,
                                           function(x) x * 0.8),
                    production = transform(prostata:::FhcrcParameters$production, values = values * 0.8))),
                highCosts = modelSet(makeModel(
                    cost_parameters = sapply(prostata:::FhcrcParameters$cost_parameters,
                                           function(x) x * 1.2),
                    production = transform(prostata:::FhcrcParameters$production, values = values * 1.2))),
                lowBiopsyCost = modelSet(makeModel(cost_parameters = (function(biopsySEK) {
                    costs <- prostata:::FhcrcParameters$cost_parameters
                    costs["Biopsy"] <- biopsySEK
                    return(costs)}
                )(biopsySEK = 3500 - 1794))), # Total 5500 (biopsy + assessment)
                highBiopsyCost = modelSet(makeModel(cost_parameters = (function(biopsySEK) {
                    costs <- prostata:::FhcrcParameters$cost_parameters
                    costs["Biopsy"] <- biopsySEK
                    return(costs)}
                )(biopsySEK = 8500 - 1794))),  # Total 8000 (biopsy + assessment)
                baseScenariosAdd = modelSet(makeModel(utility_scale = 0.0)),
                discCost0Eff0Add = modelSet(makeModel(utility_scale = 0.0,
                                                      discountRate.costs = 0.0,
                                                      discountRate.effectiveness = 0.0)),
                ## baseScenariosMultiCohort = modelSet(makeModel(popDefault = (1997 - 70) :(1997 - 55), includeBxrecords = TRUE)),
                baseScenariosMultiCohort = modelSet(makeModel(popDefault = (1997 - 70) :(1997 - 55), includeBxrecords = TRUE)),
                discCost0Eff0MultiCohort = modelSet(makeModel(popDefault = (1997 - 70) :(1997 - 55),
                                                              discountRate.costs = 0.0,
                                                              discountRate.effectiveness = 0.0))## ,
                ## indUtilLow = lapply(setNames(nm, gsub("[^[:alnum:]]", "", nm) %>% gsub("^c", "", .)),
                ##                       function(x) {
                ##                           ut <- prostata:::FhcrcParameters$utility_estimates
                ##                           ut[x] <- sapply(x, function(x) max(1 - (1 - ut[x]) * 1.2, 0.0))
                ##                           modelSet(makeModel(utility_estimates = ut))}),
                ## indUtilHigh = lapply(setNames(nm, gsub("[^[:alnum:]]", "", nm) %>% gsub("^c", "", .)),
                ##                        function(x) {
                ##                            ut <- prostata:::FhcrcParameters$utility_estimates
                ##                            ut[x] <- sapply(x, function(x) 1 - (1-ut[x]) * 0.8)
                ##                            modelSet(makeModel(utility_estimates = ut))})
                )
                cat("NOTE: Processing full parameter set completed.\n")}

save(models, file = file.path(path, paste0("one-way-sensitivity-", run_name,".RData")))

## compare <- function(m1) {
##     ICER(models$baseScenarios$modelPSA, models[[m1]]$modelPSA)
## }

## lapply(c("lowBiomarker", "highBiomarker",
## "lowBiomarker_R15", "highBiomarker_R15", "lowBiomarker_R20",
## "highBiomarker_R20", "unFavUtilities_m20", "favUtilities_p20",
## "lowCosts", "highCosts", "lowBiopsyCost", "highBiopsyCost"), function(name) list(name = name, name=compare(name)))
