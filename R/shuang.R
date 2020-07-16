ShuangParameters <- list(

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
    pTBxG0ifG1_MRIpos=0,           # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) #Updated 2020-03-24 from Bx -> TBx
    pTBxG0ifG2_MRIpos=0,           # Pr(TBx gives ISUP 0 | ISUP 2, MRI+) #Updated 2020-03-24 from Bx -> TBx
    pTBxG1ifG2_MRIpos=0,           # Pr(TBx gives ISUP 1 | ISUP 2, MRI+) (not used) #Updated 2020-03-24 from Bx -> TBx
    currency_rate = 1/10.2567,     # Riksbanken 2018
    
    mu0=c(0.002644, 0.000207, 9.5e-05, 0.00014, 0.000114, 5.6e-05, 6.4e-05, 
          9e-05, 5.5e-05, 6.4e-05, 6.1e-05, 7.1e-05, 9.1e-05, 0.000146, 
          0.000113, 0.000157, 0.000231, 0.000311, 0.000374, 0.000507, 0.000519, 
          0.000649, 0.000614, 0.000678, 0.000713, 0.000677, 0.000698, 0.000698, 
          0.000714, 0.000766, 0.000669, 0.000764, 0.00074, 0.000691, 0.000722, 
          0.000685, 0.000738, 0.000745, 0.000773, 0.000831, 0.000971, 0.001013, 
          0.001022, 0.00119, 0.001322, 0.001345, 0.001672, 0.0018, 0.002111, 
          0.002327, 0.002626, 0.002744, 0.002851, 0.003482, 0.003782, 0.004206, 
          0.004645, 0.004845, 0.005491, 0.006256, 0.00687, 0.007793, 0.008251, 
          0.009231, 0.010119, 0.011176, 0.012616, 0.013612, 0.01474, 0.016834, 
          0.018531, 0.020063, 0.021987, 0.024831, 0.028604, 0.031887, 0.03597, 
          0.040838, 0.045803, 0.050884, 0.058464, 0.06515, 0.074592, 0.08552, 
          0.096554, 0.10965, 0.123894, 0.140601, 0.155434, 0.181008, 0.201892, 
          0.227595, 0.25149, 0.28064, 0.30848, 0.344113, 0.366119, 0.419757, 
          0.436817, 0.485551, 0.569208, 0.57603, 0.622803, 0.600511, 0.786461, 
          0.823125)                # 2010-2014 death rates from human mortality database https://www.mortality.org/, 
                                   # find "Sweden" and click on 1x5 death rates https://www.mortality.org/hmd/SWE/STATS/Mx_1x5.txt, 
                                   # Access date: 2020-03-26
)
