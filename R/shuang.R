## =======================================================================================================================
## FILENAME: shuang.R
## PROJECT: cost-effectiveness of prostate cancer screening in Sweden
## PURPOSE: To provide input data as parameters for the simulation
## AUTHOR: Shuang Hao

## UPDATED: 2022-01-27

ShuangParameters <- function(year=2018) {

    stopifnot(year %in% 2018:2020)
    
    base2018=list(
        Andreas = FALSE, # general flag to use Shuang's parameters

        ## Swedish governmental report on organised PSA testing (p.22):
        ## https://www.socialstyrelsen.se/SiteCollectionDocuments/2018-2-13-halsoekonomisk-analys.pdf
        ## Based on the Swedish south region 2017:
        ## http://sodrasjukvardsregionen.se/avtal-priser/regionala-priser-och-ersattningar-foregaende-ar/
        ## S3M cost from Karolinska University Laboratory (KUL):
        ## https://www.karolinska.se/KUL/Alla-anvisningar/Anvisning/10245

        ## Latest cost for 2018
        cost_parameters =  c("Invitation" = 7                         # Invitation letter
                             + 7,                                      # Results letter
                             "Formal PSA" = 355.82                     # test sampling, primary care
                             + 57.4                                    # PSA analysis
                             + 0 * 1493.43,                            # No GP primary care
                             "Formal panel" =  355.82                  # test sampling, primary care
                             + 57.4                                    # PSA analysis not included in panel price
                             + 3300                                    # From BergusMedical (official lab for Sthlm3)
                             + 0 * 1493.43,                            # No GP for formal
                             "Opportunistic PSA" = 57.4                # PSA analysis
                             + 0.2 * 1493.43,                          # GP primary care
                             "Opportunistic panel" = 57.4              # PSA analysis not included in panel price
                             + 3300                                    # From BergusMedical (official lab for Sthlm3)
                             + 0.2 * 1493.43,                          # GP primary care
                             "Biopsy" = 3010                           # Systematic biopsy cost (SBx)
                             + 4238.25,                                # Pathology of biopsy
                             "MRI" = 3500,                             # MRI cost
                             "Combined biopsy" = 3010*1.5              # Biopsy cost (SBx|TBx) ?Double the price
                             + 4238.25,                                # Pathology of biopsy
                             "Assessment" = 1460,                      # Urologist and nurse consultation
                             "Prostatectomy" = 121170.69               # Robot assisted surgery
                             + 6302.19*20*0.25                         # Radiation therapy
                             + 1460*1,                                 # Urology and nurse visit
                             "Radiation therapy" = 6302.19*20          # Radiation therapy
                             + 3903.27*1                               # Oncologist new visit
                             + 1683.36*1                               # Oncologist further visit
                             + 400*20                                  # Nurse visit
                             + 67490.20*0.2,                           # Hormone therapy
                             "Active surveillance - yearly - w/o MRI" = 1460    # Urology visit and nurse visit
                             + 355.82*3                                # PSA sampling
                             + 57.4*3                                  # PSA analysis
                             + 3010*0.33                               # Systematic biopsy
                             + 4238.25*0.33,                           # Pathology of biopsy
                             "Active surveillance - yearly - with MRI" = 1460   # Urology visit and nurse visit
                             + 355.82*3                                # PSA sampling
                             + 57.4*3                                  # PSA analysis
                             + 3500*0.33                               # MRI cost
                             + 3010*1.5*0.33                           # Biopsy cost (SBx|TBx)
                             + 4238.25*0.33,                           # Pathology of biopsy
                             "ADT+chemo" = 71579.64*1.5,               # NEW: Chemo and hormone therapy
                             "Post-Tx follow-up - yearly first" = 1460 # Urologist and nurse consultation
                             + 355.82                                  # PSA test sampling
                             + 57.4,                                   # PSA analysis
                             "Post-Tx follow-up - yearly after" = 355.82  # PSA test sampling
                             + 57.4                                    # PSA analysis,
                             + 146,                                    # Telefollow-up by urologist
                             "Palliative therapy - yearly" = 161593.05, # Palliative care cost
                             "Terminal illness" = 161593.05*0.5),      # Terminal illness cost

        ## Swedish governmental report on organised PSA testing (p.23):
        ## https://www.socialstyrelsen.se/SiteCollectionDocuments/2018-2-13-halsoekonomisk-analys.pdf
        ## Swedish official statitics on mean salary for general population at working age
        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/arbetsmarknad/loner-och-arbetskostnader/lonestrukturstatistik-hela-ekonomin/pong/tabell-och-diagram/genomsnittlig-manadslon-efter-sektor/
        ## Percentage of social and employee contribution on salary 37.13% for abetare
        ## https://www.ekonomifakta.se/Fakta/Skatter/Skatt-pa-arbete/Sociala-avgifter-over-tid/
        ## Exchange rate from SEK to EUR from the national bank
        ## https://www.riksbank.se/sv/statistik/sok-rantor--valutakurser/arsgenomsnitt-valutakurser/?y=2018&m=12&s=Comma&f=y
        ## Consumer price index
        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/priser-och-konsumtion/konsumentprisindex/konsumentprisindex-kpi/pong/tabell-och-diagram/konsumentprisindex-kpi/kpi-faststallda-tal-1980100/
        ## Employment proportion
        ## http://www.statistikdatabasen.scb.se/pxweb/sv/ssd/START__AM__AM0401__AM0401A/NAKUBefolkning2Ar/?loadedQueryId=63385&timeType=from&timeValue=2001

        ## Updated to year 2018
        production = data.frame(ages = c(0, 54, 64, 74),
                                values=apply(
                                    rbind(0.891, 0.781, 0.169, 0), # employment portion, full time
                                    1,
                                    function(empl) empl*
                                                   30.7/40*                       # average working hours
                                                   34600*12*                      # average salary
                                                   1.3713)),                      # including non-optional social fees
        lost_production_years= c("Formal PSA"=2/24/365.25,
                                 "Formal panel"=2/24/365.25,
                                 "Opportunistic PSA"=2/24/365.25,
                                 "Opportunistic panel"=2/24/365.25,
                                 "MRI"=2/24/365.25,
                                 "Biopsy"=2/24/365.25,
                                 "Combined biopsy"=2/24/365.25,
                                 "Assessment"=2/24/365.25,         # do we need this?
                                 "Prostatectomy"=6/52,
                                 "Radiation therapy"=8/52,
                                 "Active surveillance - yearly - w/o MRI"=
                                     3 * 2/24/365.25                 # PSA tests
                                 + 1 * 2/24/365.25                 # Urologist visit
                                 + 0.33 * 2/24/365.25,             # Biopsy (SBx)
                                 "Active surveillance - yearly - with MRI"=
                                     3 * 2/24/365.25                 # PSA tests
                                 + 1 * 2/24/365.25                 # Urologist visit
                                 + 0.33 * 2/24/365.25              # MRI
                                 + 0.33 * 2/24/365.25,             # Combined biopsy (TBx|SBx)
                                 "Post-Tx follow-up - yearly"  =
                                     2/24/365.25,                    # PSA tests (Same for first year and following years)
                                 "Premature mortality" = 0,       # To discuss, depending on the age of death
                                 "Long-term sick leave" = 0.0768*67.52/365.25, # 7.68% employed PCa patients (50-64) have long-term sick leave (based on 2016 data)
                                 ## "Early retirement" = 0.00203*235.5/365.25, # 0.203% employed PCa patients (50-64) have early retirement (based on 2016 data)
                                 "Terminal illness" = 6/12),       # Should delete, should be reflected from sick leave or disability pension

        ## Latest review 2019 based on Heijnsdijk 2012, Magnus 2019 and extended review by Shuang - PORPUS-U
        utility_estimates = c("Invitation" = 1,                    # Heijnsdijk 2012
                              "Formal PSA" = 0.99,                 # Heijnsdijk 2012
                              "Formal panel" = 0.99,               # Heijnsdijk 2012
                              "Opportunistic PSA" = 0.99,          # Heijnsdijk 2012
                              "Opportunistic panel" = 0.99,        # Heijnsdijk 2012
                              "Biopsy" = 0.90,                     # Heijnsdijk 2012
                              "Cancer diagnosis" = 0.80,           # Heijnsdijk 2012
                              "Prostatectomy part 1" = 0.860,      # Magnus 2019 (Krahn 2009, Ku 2009)
                              "Prostatectomy part 2" = 0.900,      # Magnus 2019 (Krahn 2009, Ku 2009)
                              "Radiation therapy part 1" = 0.890,  # Krahn 2009
                              "Radiation therapy part 2" = 0.920,  # Krahn 2009
                              "Active surveillance" = 0.980,       # Loeb 2018
                              "Postrecovery period" = 0.930,       # Magnus 2019 (Avila 2014, Bremner 2014, Krahn 2013, Ku 2009)
                              "ADT+chemo" = 0.803,                 # Krahn 2003
                              "Palliative therapy" = 0.680,        # *15D value; Magnus 2019 (Farrkila 2014, Torvinen 2013)
                              "Terminal illness" = 0.40,           # Heijnsdijk 2012
                              "Death" = 0.00),
        ## Utility duration is given in years.
        utility_duration = c("Invitation" = 0.0,
                             "Formal PSA" = 1/52,
                             "Formal panel" = 1/52,
                             "Opportunistic PSA" = 1/52,
                             "Opportunistic panel" = 1/52,
                             "Biopsy" = 3/52,
                             "Combined biopsy" = 3/52,
                             "Cancer diagnosis" = 1/12,
                             "Prostatectomy part 1" = 2/12,
                             "Prostatectomy part 2" = 10/12,
                             "Radiation therapy part 1" = 2/12,
                             "Radiation therapy part 2" = 10/12,
                             "Active surveillance" = 7,
                             "Postrecovery period" = 9,
                             "ADT+chemo" = 1.5,                   # Assumption!!!
                             "Palliative therapy" = 12/12,        # Palliative therapy
                             "Terminal illness" = 6/12),
        pMRIposG0=0.4515496,          # Pr(MRI+ | ISUP 0 || undetectable) 2020-03-23
        pMRIposG1=0.7145305,          # Pr(MRI+ | ISUP 1 && detectable) 2020-03-23
        pMRIposG2=0.9305352,          # Pr(MRI+ | ISUP 2+ && detectable) 2020-03-23
        pSBxG0ifG1=0.1402583,         # Pr(SBx gives ISUP 0 | ISUP 1) 2020-03-23
        pSBxG0ifG2=0.1032593,         # Pr(SBx gives ISUP 0 | ISUP 2) 2020-03-23
        pSBxG1ifG2=0.119,             # Pr(SBx gives ISUP 1 | ISUP 2) (not used)
        pTBxG0ifG1_MRIpos=0,          # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) #Updated 2020-03-24 from Bx -> TBx
        pTBxG0ifG2_MRIpos=0,          # Pr(TBx gives ISUP 0 | ISUP 2, MRI+) #Updated 2020-03-24 from Bx -> TBx
        pTBxG1ifG2_MRIpos=0,          # Pr(TBx gives ISUP 1 | ISUP 2, MRI+) (not used) #Updated 2020-03-24 from Bx -> TBx
        currency_rate = 1/10.2567,    # Riksbanken 2018

        mu0=c(0.002644, 0.000207, 9.5e-05, 0.00014, 0.000114, 5.6e-05, 6.4e-05, 9e-05, 5.5e-05, 6.4e-05,
              6.1e-05, 7.1e-05, 9.1e-05, 0.000146, 0.000113, 0.000157, 0.000231, 0.000311, 0.000374, 0.000507,
              0.000519, 0.000649, 0.000614, 0.000678, 0.000713, 0.000677, 0.000698, 0.000698, 0.000714, 0.000766,
              0.000669, 0.000764, 0.00074, 0.000691, 0.000722, 0.000685, 0.000738, 0.000745, 0.000773, 0.000831,
              0.000971, 0.001013, 0.001022, 0.00119, 0.001322, 0.001345, 0.001672, 0.0018, 0.002111, 0.002327,
              0.002626, 0.002744, 0.002851, 0.003482, 0.003782, 0.004206, 0.004645, 0.004845, 0.005491, 0.006256,
              0.00687, 0.007793, 0.008251, 0.009231, 0.010119, 0.011176, 0.012616, 0.013612, 0.01474, 0.016834,
              0.018531, 0.020063, 0.021987, 0.024831, 0.028604, 0.031887, 0.03597, 0.040838, 0.045803, 0.050884,
              0.058464, 0.06515, 0.074592, 0.08552, 0.096554, 0.10965, 0.123894, 0.140601, 0.155434, 0.181008,
              0.201892, 0.227595, 0.25149, 0.28064, 0.30848, 0.344113, 0.366119, 0.419757, 0.436817, 0.485551,
              0.569208, 0.57603, 0.622803, 0.600511, 0.786461, 0.823125)
        ## 2010-2014 death rates from human mortality database https://www.mortality.org/,
        ## find "Sweden" and click on 1x5 death rates https://www.mortality.org/hmd/SWE/STATS/Mx_1x5.txt,
        ## Access date: 2020-03-26
        ## Note on 2020-12-09: due to the update from the Human Mortality database, data before age 80 kept the same from last access (0.05884)
        ## Note on 2020-12-09: Data from age80 changed slightly
    )
    
    base2019 <-
        modifyList(base2018,

                   ## Swedish governmental report on organised PSA testing (p.23 6.2 Organiserad prostatacancertestning):
                   ## https://kunskapsbanken.cancercentrum.se/diagnoser/prostatacancer/vardprogram/f
                   ## Based on the Swedish south region 2017:
                   ## http://sodrasjukvardsregionen.se/avtal-priser/regionala-priser-och-ersattningar-foregaende-ar/
                   ## S3M unit cost ("list price") from A23 Lab (Ola Steinberg):
                   ## https://a23lab.se/stockholm3/

                   ## Latest cost for 2019
                   list(cost_parameters =  c("Invitation" = 7.12                      # Invitation letter
                                             + 7.12,                                   # Results letter
                                             "Formal PSA" = 362.16                     # test sampling, primary care
                                             + 58.42                                   # PSA analysis
                                             + 0 * 1520.08,                            # No GP primary care
                                             "Formal panel" =  362.16                  # test sampling, primary care
                                             + 58.42                                   # PSA analysis not included in panel price
                                             + 2300                                    # From A23 Lab (Ola Steinberg (list price from A23 Lab)
                                             + 0 * 1520.08,                            # No GP for formal
                                             "Opportunistic PSA" = 58.42               # PSA analysis
                                             + 0.2 * 1520.08,                          # GP primary care
                                             "Opportunistic panel" = 58.42             # PSA analysis not included in panel price
                                             + 2300                                    # From BergusMedical (official lab for Sthlm3)
                                             + 0.2 * 1520.08,                          # GP primary care
                                             "Biopsy" = 3063.71                        # Systematic biopsy cost (SBx)
                                             + 4313.88,                                # Pathology of biopsy
                                             "MRI" = 3562.45,                          # MRI cost
                                             "Combined biopsy" = 3063.71*1.5           # Biopsy cost (SBx|TBx) ?Double the price
                                             + 4313.88,                                # Pathology of biopsy
                                             "Assessment" = 1486.05,                   # Urologist and nurse consultation
                                             "Prostatectomy" = 117426.74               # Robot assisted surgery
                                             + 6414.65*20*0.25                         # Radiation therapy
                                             + 1486.05*1,                              # Urology and nurse visit
                                             "Radiation therapy" = 6414.65*20          # Radiation therapy
                                             + 3972.92*1                               # Oncologist new visit
                                             + 1713.40*1                               # Oncologist further visit
                                             + 407.14*20                               # Nurse visit
                                             + 68694.51*0.2,                           # Hormone therapy
                                             "Active surveillance - yearly - w/o MRI" = 1486.05    # Urology visit and nurse visit
                                             + 362.16*3                                # PSA sampling
                                             + 58.42*3                                 # PSA analysis
                                             + 3063.71*0.33                            # Systematic biopsy (SBx)
                                             + 4313.88*0.33,                           # Pathology of biopsy
                                             "Active surveillance - yearly - with MRI" = 1486.05   # Urology visit and nurse visit
                                             + 362.16*3                                # PSA sampling
                                             + 58.42*3                                 # PSA analysis
                                             + 3562.45*0.33                            # MRI cost
                                             + 3063.71*1.5*0.33                        # Biopsy cost (SBx|TBx)
                                             + 4313.88*0.33,                           # Pathology of biopsy
                                             "ADT+chemo" = 72856.91*1.5,               # NEW: Chemo and hormone therapy
                                             "Post-Tx follow-up - yearly first" = 1486.05 # Urologist and nurse consultation
                                             + 362.16                                  # PSA test sampling
                                             + 58.42,                                  # PSA analysis
                                             "Post-Tx follow-up - yearly after" = 362.16  # PSA test sampling
                                             + 58.42                                   # PSA analysis,
                                             + 148.61,                                 # Telefollow-up by urologist
                                             "Palliative therapy - yearly" = 164476.53,# Palliative care cost
                                             "Terminal illness" = 164476.53*0.5),      # Terminal illness cost

                        ## Swedish governmental report on organised PSA testing (p.23):
                        ## https://www.socialstyrelsen.se/SiteCollectionDocuments/2018-2-13-halsoekonomisk-analys.pdf
                        ## Swedish official statitics on mean salary for general population at working age
                        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/arbetsmarknad/loner-och-arbetskostnader/lonestrukturstatistik-hela-ekonomin/pong/tabell-och-diagram/genomsnittlig-manadslon-efter-sektor/
                        ## Average working weekly working hours in Sweden 2019
                        ## https://https://www.statista.com/statistics/528482/sweden-average-weekly-working-hours/
                        ## Percentage of social and employee contribution on salary 37.00% for abetare
                        ## https://www.ekonomifakta.se/Fakta/Skatter/Skatt-pa-arbete/Sociala-avgifter-over-tid/
                        ## Exchange rate from SEK to EUR from the national bank
                        ## https://www.riksbank.se/sv/statistik/sok-rantor--valutakurser/arsgenomsnitt-valutakurser/?y=2018&m=12&s=Comma&f=y
                        ## Consumer price index
                        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/priser-och-konsumtion/konsumentprisindex/konsumentprisindex-kpi/pong/tabell-och-diagram/konsumentprisindex-kpi/kpi-faststallda-tal-1980100/
                        ## Employment proportion
                        ## http://www.statistikdatabasen.scb.se/pxweb/sv/ssd/START__AM__AM0401__AM0401A/NAKUBefolkning2Ar/?loadedQueryId=63385&timeType=from&timeValue=2001

                        ## Updated to year 2019
                        production = data.frame(ages = c(0, 54, 64, 74),
                                                values=apply(
                                                    rbind(0.89, 0.779, 0.175, 0), # employment portion, full time
                                                    1,
                                                    function(empl) empl*
                                                                   30.2/40*                    # average working hours
                                                                   34000*12*                   # average salary
                                                                   1.3700)),                   # including non-optional social fees

                        ## Updated currency exchange rate
                        currency_rate = 1/10.5892,    # Riksbanken 2019

                        mu0=c(0.002438,0.000177,0.000126,0.000068,0.000119,0.000099,0.000070,0.000042,0.000068,0.000023,
                              0.000050,0.000067,0.000089,0.000101,0.000131,0.000147,0.000166,0.000291,0.000415,0.000482,
                              0.000569,0.000602,0.000602,0.000582,0.000720,0.000823,0.000673,0.000735,0.000747,0.000699,
                              0.000698,0.000664,0.000758,0.000766,0.000718,0.000813,0.000847,0.000797,0.000785,0.000913,
                              0.001033,0.000985,0.001118,0.001084,0.001194,0.001278,0.001440,0.001436,0.001682,0.002005,
                              0.002083,0.002207,0.002719,0.002982,0.003331,0.003674,0.004121,0.004430,0.004757,0.005319,
                              0.006225,0.006676,0.007804,0.008453,0.009510,0.010282,0.011912,0.012557,0.013540,0.015319,
                              0.016865,0.018438,0.020668,0.022392,0.025004,0.028228,0.031763,0.035649,0.040985,0.045583,
                              0.052091,0.059322,0.067212,0.076669,0.087961,0.101087,0.115534,0.131556,0.147308,0.169331,
                              0.194948,0.215649,0.242328,0.270844,0.300826,0.335924,0.371850,0.411620,0.432771,0.482348,
                              0.526314,0.536479,0.665380,0.477720,0.845950,0.821827,1.135842,2.321872,6.000000)))
    ## 2015-2019 death rates from human mortality database https://www.mortality.org/,
    ## find "Sweden" and click on 1x5 death rates https://www.mortality.org/hmd/SWE/STATS/Mx_1x5.txt,
    ## Access date: 2020-12-08
    ## Note 2020-12-09: mu0 from age 0 to age 109

    base2020 <-
        modifyList(base2019,
                   ## Swedish governmental report on organised PSA testing (p.23 6.2 Organiserad prostatacancertestning):
                   ## https://kunskapsbanken.cancercentrum.se/diagnoser/prostatacancer/vardprogram/f
                   ## Based on the Swedish south region 2017:
                   ## http://sodrasjukvardsregionen.se/avtal-priser/regionala-priser-och-ersattningar-foregaende-ar/
                   ## S3M unit cost ("list price") from A23 Lab (Ola Steinberg):
                   ## https://a23lab.se/stockholm3/

                   ## Latest cost for 2020
                   list(cost_parameters =  c("Invitation" = 7.16                      # Invitation letter
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
                                             "Terminal illness" = 165293.35*0.5),      # Terminal illness cost

                        ## Swedish governmental report on organised PSA testing (p.23):
                        ## https://www.socialstyrelsen.se/SiteCollectionDocuments/2018-2-13-halsoekonomisk-analys.pdf
                        ## Swedish official statistics on mean salary for general population at working age
                        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/arbetsmarknad/loner-och-arbetskostnader/lonestrukturstatistik-hela-ekonomin/pong/tabell-och-diagram/genomsnittlig-manadslon-efter-sektor/
                        ## Average working weekly working hours in Sweden 2019
                        ## https://https://www.statista.com/statistics/528482/sweden-average-weekly-working-hours/
                        ## Percentage of social and employee contribution on salary 37.00% for abetare
                        ## https://www.ekonomifakta.se/Fakta/Skatter/Skatt-pa-arbete/Sociala-avgifter-over-tid/
                        ## Exchange rate from SEK to EUR from the national bank
                        ## https://www.riksbank.se/sv/statistik/sok-rantor--valutakurser/arsgenomsnitt-valutakurser/?y=2018&m=12&s=Comma&f=y
                        ## Consumer price index
                        ## https://www.scb.se/hitta-statistik/statistik-efter-amne/priser-och-konsumtion/konsumentprisindex/konsumentprisindex-kpi/pong/tabell-och-diagram/konsumentprisindex-kpi/kpi-faststallda-tal-1980100/
                        ## Employment proportion
                        ## http://www.statistikdatabasen.scb.se/pxweb/sv/ssd/START__AM__AM0401__AM0401A/NAKUBefolkning2Ar/?loadedQueryId=63385&timeType=from&timeValue=2001

                        ## Updated to year 2020
                        production = data.frame(ages = c(0, 54, 64, 74),
                                                values=apply(
                                                    rbind(0.884, 0.778, 0.186, 0), # employment portion, full time
                                                    1,
                                                    function(empl) empl*
                                                                   29.2/40*                    # average working hours 2020
                                                                   36100*12*                   # average salary 2020
                                                                   1.3720)),                   # including non-optional social fees 2020

                        ## Updated currency exchange rate
                        currency_rate = 1/10.4867,    # Riksbanken 2020

                        ## Updated test characteristics 2022-01-11 (based on logit transform)
                        pMRIposG0=0.184292645141448,          # Pr(MRI+ | ISUP 0 || undetectable)
                        pMRIposG1=0.316652231007673,          # Pr(MRI+ | ISUP 1 && detectable)
                        pMRIposG2=0.836863636375978,          # Pr(MRI+ | ISUP 2+ && detectable)
                        pSBxG0ifG1=0.0632272407816892,         # Pr(SBx gives ISUP 0 | ISUP 1)
                        pSBxG0ifG2=0.0990348918673406         # Pr(SBx gives ISUP 0 | ISUP 2)
                        ))

    if(year==2018) base2018 else (if(year==2019) base2019 else base2020)
}

## Update the background health state values using Burström (2001, for men in Sweden) - Previously Burström (2006, for general population in Stockholm)
ShuangTables <- list(background_utilities =
                         data.frame(lower=c(0, 20, 30, 40, 50, 60, 70, 80),
                                    upper=c(20, 30, 40, 50, 60, 70, 80, 1.0e55),
                                    utility=c(1, 0.91, 0.90, 0.86, 0.84, 0.83, 0.81, 0.74)))

## Code to compare results
compareParameters <- function(year1,year2) {
    p1 = ShuangParameters(year1)
    p2 = ShuangParameters(year2)
    n1 = names(p1)
    n2 = names(p2)
    if (any(diff <- setdiff(n1,n2))) print(diff)
    if (any(diff <- setdiff(n2,n1))) print(diff)
    for (name in n1[n1 %in% n2]) {
        if (length(p1[[name]]) != length(p2[[name]])) {
            print(name)
            cat("Different lengths:\n")
            print(p1[[name]])
            print(p2[[name]])
        } else if(any(p1[[name]] != p2[[name]])) {
            print(name)
            cat("Differences:\n")
            print(p1[[name]]-p2[[name]])
        }
    }
}

## compareParameters(2018,2019)
## compareParameters(2019,2020)
## ShuangParameters(2018)$cost_parameters
## ShuangParameters(2019)$cost_parameters
## ShuangParameters(2020)$cost_parameters
