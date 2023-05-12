## Parameters that are specific to Trust's analyses

## Strategies

TrustParameters <- function(year=2020, MRI_screen=TRUE) {
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
                    MRI_screen = MRI_screen,
                    active_surveillance_cost_scale_first_two_years = if (MRI_screen) 1214.41/455.83 else 693.40/283.49,
                    cost_parameters =  c("Invitation" = 0,                         # Invitation letter
                                         "Formal PSA" = 29.15,                        # 
                                         "Formal panel" =  29.15,                     # 
                                         "Opportunistic DRE" = 29.15,              # Cost for DRE alone
                                         "Opportunistic PSA" = 43.95-29.15,        # NB: PSA analysis *less* the cost for DRE
                                         "Opportunistic panel" = 29.15,               # 
                                         "Biopsy" = 358.11,                        # Systematic biopsy
                                         "MRI" = 121.01,                           # MRI cost
                                         "Combined biopsy" = 472.18 - 121.01,      # Biopsy cost (NB: assumes a preceding MRI)
                                         "Assessment" = 19.54,                     # Urologist consultation
                                         "Prostatectomy" = 10927.37,               # Robot assisted surgery
                                         "Radiation therapy" = 8872.83,            # Radiation therapy
                                         "Active surveillance - yearly - w/o MRI" = 283.49, # Urology visit and nurse visit
                                         "Active surveillance - yearly - with MRI" = 455.83, # Urology visit and nurse visit
                                         "ADT+chemo" = 145216/10,                  # Drug treatment for metastasis
                                         "Post-Tx follow-up - yearly first" = 2*410.14 - 144.57, # Urologist and nurse consultation (double counts for costs in the first two years and adjusts)
                                         "Post-Tx follow-up - yearly after" = 144.57,  # PSA test sampling
                                         "Palliative therapy - yearly" = 165293.35/10, # Palliative care cost
                                         "Terminal illness" = 165293.35*0.5/10       # Terminal illness cost
                                         ),
                    ## Jansen et al 2014 (TTO)
                    background_utilities =
                        data.frame(lower = c(0, 18, 25, 35, 45, 55, 65, 75),
                                   upper = c(18, 25, 35, 45, 55, 65, 75, 1.0e55),
                                   utility = c(1, 0.972, 0.973, 0.966, 0.945, 0.922, 0.891, 0.839)),
                    ## default treatment assignment: 2021 guidelines (assuming CM from age 75)
                    prtx =
                        data.frame(DxY=2008,
                                   Age=c(0,0,0,75,75,75),
                                   G=as.integer(c(6,7,8,6,7,8)-6),
                                   CM=c(1,0,0,1,1,1),
                                   RP=c(0,0.5,0.5,0,0,0),
                                   RT=c(0,0.5,0.5,0,0,0)),
                    weibull_onset = TRUE,
                    ## fitted parameters
                    weibull_onset_shape=exp(-0.0754039582652198),
                    weibull_onset_scale=exp(5.2319848747848),
                    beta7=exp(-2.35969408469582),
                    beta8=exp(-1.30383055358671),
                    mu0=c(0.003233, 0.000193, 0.0001, 0.000096, 0.000113, 0.000075, 0.00008, 0.000102, 0.000082, 0.000063, 0.000058, 0.000063, 0.000077, 0.000092, 0.000113, 0.000119, 0.000231, 0.000266, 0.000328, 0.000434, 0.000407, 0.000397, 0.000402, 0.000448, 0.000457, 0.000439, 0.000437, 0.00043, 0.000401, 0.000512, 0.00054, 0.00055, 0.000595, 0.000683, 0.000713, 0.000774, 0.000928, 0.001085, 0.001051, 0.001192, 0.001304, 0.001444, 0.001411, 0.001653, 0.001817, 0.001847, 0.002028, 0.002231, 0.002671, 0.002995, 0.003273, 0.003604, 0.003926, 0.004485, 0.004965, 0.005712, 0.00638, 0.006836, 0.007844, 0.008599, 0.009509, 0.010741, 0.011794, 0.013032, 0.014189, 0.015759, 0.016911, 0.018351, 0.020393, 0.021884, 0.0229, 0.025254, 0.026996, 0.029723, 0.031317, 0.036081, 0.039332, 0.041678, 0.046748, 0.049817, 0.059204, 0.064676, 0.071599, 0.082113, 0.092783, 0.110798, 0.120641, 0.14056, 0.159357, 0.183627, 0.204152, 0.226824, 0.253886, 0.281857, 0.316669, 0.346189, 0.375991, 0.402264, 0.441481, 0.494507, 0.582567, 0.492796, 0.550753, 0.717067, 0.516318, 0.80061, 0.673725, 0.699847, 0.693939, 0.980952)))
}




## read_string = function(string, header=TRUE, sep="\t", ...)
##     read.table(text=string, header=header, sep=sep, ...)

germany_observed_tables = 
    list(prtx = data.frame(DxY = 2008,
                           Age = c(50L, 
                                   50L, 50L, 60L, 60L, 60L, 70L, 70L, 70L, 80L, 80L, 80L),
                           G = c(0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L),
                           CM = c(0.18333333, 0.11052632, 0.325, 0.26923077, 0.12541806, 0.375, 0.3625498,
                                  0.31076389, 0.53040541, 0.53932584, 0.56737589, 0.73188406), 
                           RP = c(0.65, 0.81578947, 0.65, 0.50480769, 0.77926421, 0.5625, 
                                  0.39043825, 0.50694444, 0.39527027, 0.3258427, 0.25531915, 
                                  0.20289855),
                           RT = c(0.16666667, 0.07368421, 0.025, 0.22596154, 
                                  0.09531773, 0.0625, 0.24701195, 0.18229167, 0.07432432, 0.13483146, 
                                  0.17730496, 0.06521739)))


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
