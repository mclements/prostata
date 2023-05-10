## Parameters that are specific to Trust's analyses

## Strategies

TrustParameters <- function(year=2020) {
    modifyList(ShuangParameters(year),
               list(start_screening = 45,
                    stop_screening = 75,
                    ## PSA<1 => 4-yearly; 1<=PSA<2 => 2-yearly; PSA>=2 => 1-yearly
                    risk_psa_threshold_lower = 1,
                    risk_psa_threshold_moderate = 2,
                    risk_lower_interval = 4,
                    risk_moderate_interval = 2,
                    risk_upper_interval = 1,
                    psaThreshold = 4,
                    psaThresholdBiopsyFollowUp = 4, # is this correct??
                    MRI_screen = TRUE,
                    cost_parameters =  c("Invitation" = 7.16                      # Invitation letter
                                         + 7.16,                                   # Results letter
                                         "Formal PSA" = 363.96                     # test sampling, primary care
                                         + 58.71                                   # PSA analysis
                                         + 0 * 1527.63,                            # No GP primary care
                                         "Formal panel" =  363.96                  # test sampling, primary care
                                         + 58.71                                   # PSA analysis not included in panel price
                                         + 2300                                    # From A23 Lab (Ola Steinberg (list price from A23 Lab)
                                         + 0 * 1527.63,                            # No GP for formal
                                         "Opportunistic PSA" = 58.71               # PSA analysis
                                         + 0.2 * 1527.63,                          # GP primary care
                                         "Opportunistic panel" = 58.71             # PSA analysis not included in panel price
                                         + 2300                                    # From BergusMedical (official lab for Sthlm3)
                                         + 0.2 * 1527.63,                          # GP primary care
                                         "Biopsy" = 3010                           # Systematic biopsy cost (SBx)
                                         + 4335.30,                                # Pathology of biopsy
                                         "MRI" = 3500,                             # MRI cost
                                         "Combined biopsy" = 3010*1.5              # Biopsy cost (SBx|TBx) ?Double the price
                                         + 4335.30,                                # Pathology of biopsy
                                         "Assessment" = 1460,                      # Urologist and nurse consultation
                                         "Prostatectomy" = 118009.91               # Robot assisted surgery
                                         + 6446.51*20*0.25                         # Radiation therapy
                                         + 1460*1,                                 # Urology and nurse visit
                                         "Radiation therapy" = 6446.51*20          # Radiation therapy
                                         + 3992.65*1                               # Oncologist new visit
                                         + 1721.90*1                               # Oncologist further visit
                                         + 400*20                                  # Nurse visit
                                         + 69035.66*0.2,                           # Hormone therapy
                                         "Active surveillance - yearly - w/o MRI" = 1460    # Urology visit and nurse visit
                                         + 363.96*3                                # PSA sampling
                                         + 58.71*3                                 # PSA analysis
                                         + 3010*0.33                               # Systematic biopsy (SBx)
                                         + 4335.30*0.33,                           # Pathology of biopsy
                                         "Active surveillance - yearly - with MRI" = 1493.43   # Urology visit and nurse visit
                                         + 363.96*3                                # PSA sampling
                                         + 58.71*3                                 # PSA analysis
                                         + 3500*0.33                               # MRI cost
                                         + 3010*1.5*0.33                           # Biopsy cost (SBx|TBx)
                                         + 4335.30*0.33,                           # Pathology of biopsy
                                         "ADT+chemo" = 145216,                     # Drug treatment for metastasis
                                         "Post-Tx follow-up - yearly first" = 1460 # Urologist and nurse consultation
                                         + 363.96                                  # PSA test sampling
                                         + 58.71,                                  # PSA analysis
                                         "Post-Tx follow-up - yearly after" = 363.96  # PSA test sampling
                                         + 58.71                                   # PSA analysis,
                                         + 146,                                    # Telefollow-up by urologist
                                         "Palliative therapy - yearly" = 165293.35,# Palliative care cost
                                         "Terminal illness" = 165293.35*0.5,       # Terminal illness cost
                                         "Opportunistic DRE" = 349                 # test sampling, primary care
                                         + 0.2 * 1539                              # GP primary care
                                         )))
}

if (FALSE) {
    library(prostata)
    set.seed(12345)
    germany_2018 = callFhcrc(n=1e5, screen="germany_2018",
                             parms=prostata:::TrustParameters(),
                             mc.core=6)
    set.seed(12345)
    regular_screen_4 = callFhcrc(n=1e5, screen="regular_screen",
                                 parms=modifyList(prostata:::TrustParameters(),
                                                  list(screening_interval=4)),
                                 mc.cores=6)
    set.seed(12345)
    noScreening = callFhcrc(n=1e5, screen="noScreening",
                            parms=prostata:::TrustParameters(),
                            mc.cores=6)

    quickTab = function(list) {
        do.call(rbind,
                lapply(list, function(item) summary(item)[c("QALE","healthsector.costs")]))
    }

    tab = quickTab(list(noScreening=noScreening,
                       germany_2018=germany_2018,
                       regular_screen_4=regular_screen_4))
    plot(tab)
    text(tab, labels=rownames(tab), pos=c(4,2,2))

    ## Expected number of PSA tests = number of DRE
    ## The number of "toScreen" is _only_ the PSA tests for those with a negative DRE

    set.seed(12345)
    temp = callFhcrc(n=1e5, screen="germany_2018",
                     parms=prostata:::TrustParameters(),
                     mc.core=6, nLifeHistories=1e5)

    subset(temp$lifeHistories, event == "toScreenInitiatedBiopsy" & begin == 45) |> head()
    subset(temp$lifeHistories, event == "toScreenDiagnosis" & end<=50) |> nrow()
    table(subset(temp$lifeHistories, begin == 45)$event)

    subset(temp$lifeHistories, event == "toDRE" & begin == 45)
    table(temp$lifeHistories$event)
    subset(temp$lifeHistories, event == "toScreenInitiatedBiopsy" & end == 45)
    subset(temp$lifeHistories, id == 1882)
    table(temp$lifeHistories$event)

    ## reconstruction for ERSPC first wave at age 55 years
    set.seed(12345)
    temp_4 = callFhcrc(n=1e5, screen="regular_screen",
                       parms=modifyList(prostata:::TrustParameters(),
                                        list(MRI_screen=FALSE,screening_interval=4, start_screening=55, stop_screening=75)),
                       mc.cores=6, nLifeHistories=1e5)

    table(temp_4$lifeHistories$event)
    tmp = subset(temp_4$lifeHistories, end < 55.07693 & event == "toScreen")
    xtabs(~ext_state,tmp)

    tmp = subset(temp_4$lifeHistories, end < 55.07693 & end >= 55 & event == "toScreenInitiatedBiopsy") # biopsies at first screen
    sensitivity = function(psa) prostata:::dre$sensitivity[findInterval(psa, prostata:::dre$psa_low)]
    specificity = function(psa) prostata:::dre$specificity[findInterval(psa, prostata:::dre$psa_low)]
    sensitivity = function(psa) psa*0+0.51
    specificity = function(psa) psa*0+0.59

    sensitivity = function(psa) psa*0+0.993
    specificity = function(psa) psa*0+0.0257
    ## sensitivity(c(0.5,1.5,1500))
    tmp = transform(tmp,
                    pdrepos = ifelse(ext_state=="Healthy", 1-specificity(psa), sensitivity(psa)),
                    pdreneg = ifelse(ext_state=="Healthy", specificity(psa), 1-sensitivity(psa)))
    ppv = function(x) sum(x[2:4])/sum(x)
    xtabs(pdrepos~ext_state, tmp) |> ppv()
    xtabs(pdreneg~ext_state, tmp) |> ppv()

    objective = function(beta) {
        sensitivity = beta[1]
        specificity = beta[2]
        tmp = transform(tmp,
                        pdrepos = ifelse(ext_state=="Healthy", 1-specificity, sensitivity),
                        pdreneg = ifelse(ext_state=="Healthy", specificity, 1-sensitivity))
        ppv = function(x) sum(x[2:4])/sum(x)
        (ppv(xtabs(pdrepos~ext_state, tmp))-0.486)^2+
            (ppv(xtabs(pdreneg~ext_state, tmp))-0.224)^2
    }
    optim(c(0.6,0.7),objective, method="L-BFGS-B", lower=0, upper=1)
        
    
    xtabs(~ext_state, tmp)
}
