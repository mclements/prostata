##==============================================================================
## FILENAME: daniela.R
## PROJECT:  Cost-effectiveness of prostate cancer screening using MRI: based on evidence from the STHLM3-MRI and Göteborg-2 trial (Study 2)
## PURPOSE:  This script is to store the parameters needed to feed the model, allowing for selection between strategies
##    Parameters to include:
##      Test characteristics (varying from diff strategies)
##      Cost per procedure (weighted unit costs 2024 out of three regions. Stockholm, Skåne, Västra Götaland)
##      Production
##        Employment proportion*salary*working hours*social fee
##      Utility value 
##      Utility duration
##      Mortality rate
##      Currency rate
##      Participation/Rescreening participation/Biopsy compliance rate

## AUTHOR:   Daniela Skalt

## CREATED:	2024-01-15
## UPDATED: 2024-01-22
## UPDATED: 2024-04-04
## UPDATED: 2025-08-04
## UPDATED: 2025-08-19
## R VERSION: R-4.4.3

##==============================================================================

DanielaParameters <- function(year=2024,...) {
    
    stopifnot(year %in% 2024)
  
    # Parameters including: cost, lost production, utility value and utility duration
    params <-
      list(
        Andreas = FALSE,
        ## Weighted costs, 2024, krona
           #Weighted unit costs of the three regions Stockholm, Skåne and Västra Götaland (calculated based on the population size and then used for the model) 
        cost_parameters =   c("Invitation" = 7                                # Invitation letter
                              + 7,                                            # Result letter (referral to MRI/urologist)
                              "Formal PSA" =                                  # PSA test (OPT setting)
                                65.38                                         # PSA test sampling
                              + 216.32                                        # PSA analysis
                              + 2045.09*0,                                    # No GP primary care
                              "Opportunistic PSA" =                           # PSA test (Opportunistic setting) 
                                65.38                                         # PSA test sampling
                              + 216.32                                        # PSA analysis
                              + 2045.09*0.2,                                  # GP primary care
                              "Assessment" = 2305.78,                         # Urologist and nurse consultation 
                              "MRI" = 3548.67,                                # MRI cost, OPT-specific
                              "MRIpos" = 2711,                                # Writing on the image  
                              #"MRI" =                                        # MRI cost, opportunistic
                              "Biopsy" =                                     # SBx cost (used if no MRI, not relevant in our case)
                               2777.79                                     # Systematic biopsy (not used, because all strategies use MRI)
                               + 4367.54                                     # Pathology incl digital format
                               + 2261.19,                                    # Urologist visit regarding result of biopsy
                              "Combined biopsy" =                             # TBx cost (Bad name, modify to actual combined biopsy when needed in sensitivity analysis)
                              ##  2711*0.38                                   # Writing on the image  
                              2912.30                                         # TBx  
                              + 4367.54                                       # Pathology incl digital format - 10 cores (The only difference with TBx/SBx)   
                              + 2261.19                                       # Urology visit regarding result of biopsy
                              + 3796*0.006,                                   # MRI in high risk PCa patients for lymph node evaluation
                              #"Combined biopsy" =                            # SBx+TBx cost (Combined = 1.5*SBx; default)(only used in sensitivity analysis)
                              ##   2711*0.38                                  # Writing on the image, according to OPT report only 38% (not used anymore) 
                              # 2777.79*1.5                                   # 1.5 times the price of SBx
                              # + 4791.39                                     # Pathology incl digital format - 12 cores
                              # + 2261.19                                     # Urology visit regarding result of biopsy
                              # + 3796*0.006,                                 # MRI in high risk PCa patients for lymph node evaluation
                              #"Active surveillance - yearly - w/o MRI"       # not used in our study because all strategies use MRI
                              "Active surveillance - yearly - with MRI" =     # Active surveillance with TBx
                                3992.04*1                                     # Urologist and nurse consultation
                              + 355.40*3                                      # PSA test sampling, three times a year  
                              + 64.64*3                                       # PSA analysis, three times a year
                              + 3805.10*0.33                                  # MRI scan, not OPT-specific, every 3 year
                              + 2679*0.33                                     # Writing on the image
                              + 3311.66*0.33                                  # TBx
                              + 4687.98*0.33,                                 # Pathology incl digital format
                              #SBx/TBx, not OPT-specific
                              # + 3177.15*0.33*1.5                            # SBx/TBx: Double the price of SBx (relevant for sensitivity analyses)
                              # + 5111.82*0.33,                               # Pathology incl digital format
                              "Prostatectomy" = 
                                128627.65*1                                   # Robot assisted surgery
                              + 7499.60*20*0.25                               # Radiation therapy (20 times during the whole course, 25% recurrence) 
                              + 3992.04*1,                                    # Specialist and nurse consultation
                              "Radiation therapy" = 
                                7499.60*20                                    # Radiation therapy, 20 times during the whole course
                              + 4635.23*1                                     # Oncologist new visit // Included data - new
                              + 3965.32*1                                     # Oncologist further visit // Included data - new
                              + 557.20*20                                     # Nurse visit, 20 times
                              + 85318.09*0.2,                                 # Hormone therapy, 20% patients
                              "ADT+chemo" = 179466.61,                        # Metastasis (hormone + chemo therapy)
                              "Post-Tx follow-up - yearly first" = 
                                3992.04                                       # Specialist and nurse consultation
                              + 355.40                                        # PSA test sampling
                              + 64.64,                                        # PSA test analysis
                              "Post-Tx follow-up - yearly after" = 
                                1352.27                                       # Tele follow-up
                              + 355.40                                        # PSA test sampling
                              + 64.64,                                        # PSA test analysis
                              "Palliative therapy - yearly" = 204279.40,      # Palliative care cost-yearly
                              "Terminal illness" = 102139.70),                # Terminal illness cost-yearly
             
            
            production = data.frame(ages = c(0, 54, 64, 74),                  # Age group: 45-54,55-64,65-74,74+
                                    values=apply(
                                        rbind(0.888,0.781,0.207, 0), #modify  # employment proportion, full time 2024 (link is not working anymore, use another one)
                                        1,
                                        function(empl) empl*
                                                       29.9/40*    #check     # average working hours 2024 -> is it total employment not stratified by age and sex?
                                                       41600*12*              # average monthly salary 2024
                                                       1.373)),               # including non-optional social fees 2024
            
            ## average salary 2024
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0110__AM0110B/LonYrkeUtbildningAN/table/tableViewLayout1/ 
            ## employment proportion 2023 --> can't find the page, will leave 2023 for now
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0401__AM0401A/NAKUBefolkning2Ar/
            ## average working hours 2024
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0401__AM0401S/NAKUFaktMedArbtidAr/table/tableViewLayout1/
            ## Social fee 2024
            ## https://www.ekonomifakta.se/Fakta/skatt/Skatt-pa-arbete/Sociala-avgifter-over-tid/
            ## https://www.scb.se/hitta-statistik/statistik-efter-amne/priser-och-ekonomiska-tendenser/priser/konsumentprisindex-kpi/pong/tabell-och-diagram/konsumentprisindex-kpi/inflation-i-sverige/
            ## yearly CPI up to 2024
            ## https://www.scb.se/hitta-statistik/statistik-efter-amne/priser-och-ekonomiska-tendenser/priser/konsumentprisindex-kpi/pong/tabell-och-diagram/konsumentprisindex-kpi/kpi-faststallda-tal-1980100/
            ## exchange rate, Riksbank 2024
            ## https://www.riksbank.se/en-gb/statistics/interest-rates-and-exchange-rates/search-interest-rates-and-exchange-rates/?s=g130-SEKEURPMI&s=g130-SEKUSDPMI&a=Y&from=2024-01-02&to=2024-12-30&fs=3#result-section
      
            
            
            ## Based on Hao et al 2022, JAMA Oncology
            lost_production_years= c(
                "Formal PSA"                                  =  2/24/365.25,
                "Opportunistic PSA"                           =  2/24/365.25,
                #"Formal panel"                               =  2/24/365.25,
                #"Opportunistic panel"                        =  2/24/365.25,
                "MRI"                                         =  2/24/365.25,
                "Biopsy"                                     =  2/24/365.25,         # SBx, not needed in our analysis
                "Combined biopsy"                             =  2/24/365.25,
                "Assessment"                                  =  2/24/365.25,                        
                "Prostatectomy"                               =  6/52,
                "Radiation therapy"                           =  8/52,
                #"Active surveillance - yearly - w/o MRI"     =  3*2/24/365.25        # PSA tests
                #                                              + 1*2/24/365.25        # Urologist visit
                #                                              + 0.33*2/24/365.25,    # Biopsy (SBx)
                "Active surveillance - yearly - with MRI"     =  3*2/24/365.25        # PSA tests
                                                               + 1*2/24/365.25        # Urologist visit
                                                               + 0.33*2/24/365.25     # MRI
                                                               + 0.33*2/24/365.25,    # TBx OR Combined biopsy (TBx|SBx)
                "Post-Tx follow-up - yearly"                  =  2/24/365.25,         # PSA tests 
                "Premature mortality"                         =  0,                   # Keep it following Shuang's code
                "Long-term sick leave"                        =  0.0768*67.52/365.25,
                ## "Early retirement"                         =  0.00203*235.5/365.25,# 0.203% employed PCa patients (50-64) have early retirement (based on 2016 data)
                "Terminal illness"                            =  6/12),                 # 7.68% employed PCa patients (50-64) have long-term sick leave (based on 2016 data)
                       
                                                                                              
            

            ## Based on Hao et al 2022 and its references
            utility_estimates = c(
                "Invitation"                                  = 1,                             # Heijnsdijk 2012
                "Formal PSA"                                  = 0.99,                          # Heijnsdijk 2012
                "Opportunistic PSA"                           = 0.99,                          # Heijnsdijk 2012
                "Formal panel"                               = 0.99,                          # Heijnsdijk 2012
                "Opportunistic panel"                      = 0.99,                          # Heijnsdijk 2012
                "Biopsy"                                     = 0.90,                           # Heijnsdijk 2012
                "Cancer diagnosis"                            = 0.80,                          # Heijnsdijk 2012
                "Prostatectomy part 1"                        = 0.860,                         # Magnus 2019 (Krahn 2009, Ku 2009)
                "Prostatectomy part 2"                        = 0.900,                         # Magnus 2019 (Krahn 2009, Ku 2009)
                "Radiation therapy part 1"                    = 0.890,                         # Krahn 2009
                "Radiation therapy part 2"                    = 0.920,                         # Krahn 2009
                "Active surveillance"                         = 0.980,                         # Loeb 2018
                "Postrecovery period"                         = 0.930,                         # Magnus 2019 (Avila 2014, Bremner 2014, Krahn 2013, Ku 2009)
                "ADT+chemo"                                   = 0.803,                         # ADT+chemo Krahn 2003
                "Palliative therapy"                          = 0.680,                         # *15D value; Magnus 2019 (Farrkila 2014, Torvinen 2013)
                "Terminal illness"                            = 0.40,                          # Heijnsdijk 2012
                "Death"                                       = 0.00),
              
            
            
            ## Utility duration is given in years, based on Hao et al 2022
            utility_duration = c(
                "Invitation"                                  = 0.0,
                "Formal PSA"                                  = 1/52,
                "Opportunistic PSA"                           = 1/52,                           
                "Formal panel"                               = 1/52,                           
                "Biopsy"                                     = 3/52,
                "Opportunistic panel"                      = 1/52,
                "Combined biopsy"                             = 3/52,
                "Cancer diagnosis"                            = 1/12,
                "Prostatectomy part 1"                        = 2/12,
                "Prostatectomy part 2"                        = 10/12,
                "Radiation therapy part 1"                    = 2/12,
                "Radiation therapy part 2"                    = 10/12,
                "Active surveillance"                         = 7,
                "Postrecovery period"                         = 9,
                "ADT+chemo"                                   = 1.5,                                              
                "Palliative therapy"                          = 12/12,        
                "Terminal illness"                            = 6/12),
            
        
                # Strategy 2: PSA3 + MRI3 + TBx (STHLM3-MRI)
                pMRIposG0=0.184292645141448,                    # Pr(MRI+ | ISUP 0 || undetectable)
                pMRIposG1=0.316652231007673,                    # Pr(MRI+ | ISUP 1 && detectable)
                pMRIposG2=0.836863636375978,                    # Pr(MRI+ | ISUP 2+ && detectable)
                pMRIposG0_r2=0.0710412948807811,                # Pr(MRI+ | ISUP 0 || undetectable)
                pMRIposG1_r2=0.101787194055774,                 # Pr(MRI+ | ISUP 1 && detectable)
                pMRIposG2_r2=0.656413754041609,                 # Pr(MRI+ | ISUP 2+ && detectable)
                pTBxG0ifG1_MRIpos=0.542857142857143,            # Pr(TBx gives ISUP 0 | ISUP 1, MRI+)
                pTBxG0ifG2_MRIpos=0.0621468926553672,           # Pr(TBx gives ISUP 0 | ISUP 2-3, MRI+) -- NB: actually G2 and G3
                pTBxG0ifG4plus_MRIpos=0.0621468926553672,       # Pr(TBx gives ISUP 0 | ISUP 4+, MRI+)
                pTBxG0ifG1_MRIpos_r2=0.333333333333333,         # Pr(TBx gives ISUP 0 | ISUP 1, MRI+)
                pTBxG0ifG2_MRIpos_r2=0.0666666666666667,        # Pr(TBx gives ISUP 0 | ISUP 2-3, MRI+) -- NB: actually G2 and G3
                pTBxG0ifG4plus_MRIpos_r2=0.0666666666666667,    # Pr(TBx gives ISUP 0 | ISUP 4+, MRI+)
                #pSBxG0ifG1_MRIpos=    0.06322724               # Pr(SBx gives ISUP 0 | ISUP 1, MRI+)
                #psBxG0ifG2_MRIpos=    0.09903489               # Pr(SBx gives ISUP 0 | ISUP 2)
                #pSBxG0ifG1_MRIpos_r2= 0.03391762               # Pr(SBx gives ISUP 0 | ISUP 1, MRI+) 
                #psBxG0ifG2_MRIpos_r2= 0.11669378               # Pr(SBx gives ISUP 0 | ISUP 2)
                # 3. Strategy 3: PSA3 + MRI3 + TBx (G2)
                # pMRIposG0= 0.231610631,                     # Pr(MRI+ | ISUP 0 || undetectable)
                # pMRIposG1= 0.454778423,                     # Pr(MRI+ | ISUP 1 && detectable)
                # pMRIposG2= 0.869292977,                     # Pr(MRI+ | ISUP 2+ && detectable)
                # pMRIposG0_r2= 0.041978136,                  # Pr(MRI+ | ISUP 0 || undetectable)
                # pMRIposG1_r2= 0.351186586,                  # Pr(MRI+ | ISUP 1 && detectable)
                # pMRIposG2_r2= 0.422062828,                  # Pr(MRI+ | ISUP 2+ && detectable)
                # pTBxG0ifG1_MRIpos= 0.138888889,             # Pr(TBx gives ISUP 0 | ISUP 1, MRI+)
                # pTBxG0ifG2_MRIpos= 0.018867925,             # Pr(TBx gives ISUP 0 | ISUP 2-3, MRI+) -- NB: actually G2 and G3
                # pTBxG0ifG4plus_MRIpos= 0.018867925,         # Pr(TBx gives ISUP 0 | ISUP 4+, MRI+)
                # pTBxG0ifG1_MRIpos_r2= 0.350000000,          # Pr(TBx gives ISUP 0 | ISUP 1, MRI+)
                # pTBxG0ifG2_MRIpos_r2= 0.000000000,          # Pr(TBx gives ISUP 0 | ISUP 2-3, MRI+) -- NB: actually G2 and G3
                # pTBxG0ifG4plus_MRIpos_r2= 0.000000000,      # Pr(TBx gives ISUP 0 | ISUP 4+, MRI+)  
                #pSBxG0ifG1_MRIpos=    xxxx                 # Pr(SBx gives ISUP 0 | ISUP 1, MRI+)
                #psBxG0ifG2_MRIpos=    xxxx                 # Pr(SBx gives ISUP 0 | ISUP 2)
                #pSBxG0ifG1_MRIpos_r2= xxxx                 # Pr(SBx gives ISUP 0 | ISUP 1, MRI+) 
                #psBxG0ifG2_MRIpos_r2= xxxx                 # Pr(SBx gives ISUP 0 | ISUP 2)
                currency_rate = 1/11.4322,               # Riksbanken 2024, Krona-EUR
                #currency rate = 1/10.5614               # Riksbanken 2024, Krona-USD
                screeningParticipation = 0.80,   # probability of actually having the first PSA test, 
                rescreeningParticipation = 0.80, # probability of actually having the re-screening PSA tests
                ##source:
                # 85% based on participation in Breast cancer screening
                # Probability of re-screening was 84% in Finnish trial, so we decide to take the same value as for the first round (85%),
                # Finnish trial first screening round was 69%
                biopsyCompliance = 0.95,
                    
        
            ## https://www.mortality.org/ 2020-2024, 1*5 Death rates for Swedish male, #accessed on July 30th, 2025
            mu0=c(0.002263, 0.000183, 0.000107, 0.000101, 0.000112, 0.000057, 0.000076, 0.000069, 0.000072, 
                  0.000075, 0.000056, 0.000083, 0.000087, 0.000112, 0.000144, 0.000169, 0.000228, 0.000274, 
                  0.000353, 0.000521, 0.000485, 0.000553, 0.000547, 0.000572, 0.000544, 0.000631, 0.000577,
                  0.000586, 0.000730, 0.000653, 0.000611, 0.000615, 0.000642, 0.000592, 0.000685, 0.000728,
                  0.000703, 0.000725, 0.000787, 0.000819, 0.000899, 0.000927, 0.000968, 0.001083, 0.001088,
                  0.001204, 0.001212, 0.001554, 0.001492, 0.001688, 0.001907, 0.002130, 0.002319, 0.002428,
                  0.003032, 0.003314, 0.003773, 0.003891, 0.004365, 0.005027, 0.005646, 0.006001, 0.006846,
                  0.007719, 0.008783, 0.009607, 0.010895, 0.011612, 0.013254, 0.014639, 0.015929, 0.018023,
                  0.019435, 0.021822, 0.024277, 0.026833, 0.029521, 0.033638, 0.037920, 0.042076, 0.048215,
                  0.054000, 0.062416, 0.071760, 0.081748, 0.094843, 0.106288, 0.122754, 0.142449, 0.161181, 
                  0.185782, 0.212148, 0.237919, 0.264058, 0.301515, 0.331116, 0.370883, 0.416982, 0.441188, 
                  0.466331, 0.545220, 0.591121, 0.560239, 0.668402, 0.718562, 0.726336, 0.798525, 0.834050, 0.540379))
        
        
    modifyList(params, list(...))

}


## Background health state values using Teni 2021
DanielaTables <- list(background_utilities =
                                 data.frame(lower=c(0, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
                                            upper=c(30, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 1.0e55),
                                            utility=c(1, 0.925, 0.938, 0.930, 0.924, 0.914, 0.910, 0.910, 0.915, 0.909, 0.892, 0.865, 0.831, 0.803, 0.751)))

