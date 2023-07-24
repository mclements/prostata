EdnaParameters <- function(base = ShuangParameters())
  modifyList(base,
             list(includeEventHistories=FALSE,
                  beta7=0.039052,
                  beta8=0.245344,
                  fullUptakePortion = 0.6,
                  startFullUptake = 1932.0, # arbitrary
                  yearlyUptakeIncrease = 0.0, # no change
                  uptakeStartAge=45.0,
                  startUptakeMixture = 1950.0,
                  endUptakeMixture = 1950.0, # = startUptakeMixture
                  screeningIntroduced = 1995.0,
                  shapeA = 4.8,
                  scaleA = 15.0,
                  shapeT = 4.8,  # = shapeA
                  scaleT = 15.0, # = scaleA
                  cap_pScreened = c(0.3235438, # 50-54
                      0.3531057, # 55-59
                      0.3582405, # 60-64
                      0.3236445),# 65-59
                  screeningParticipation = 0.85,
                  biopsyCompliance = 0.85,
                  eol=1,
                  MRI_screen = TRUE, MRI_clinical = FALSE,
                  MRInegSBx= FALSE,              # No SBx for MRI- (by default)
                  RP_mortHR = 0.63,  #(95% CI, 0.21 to 1.93)
                  weibull_onset = TRUE,
                  ## weibull_onset_shape = 0.8164337, # updated 2021-07-27
                  ## weibull_onset_scale= 92.3742713, # updated 2021-07-27
                  weibull_onset_shape = 0.8130842, # updated 2021-07-28
                  weibull_onset_scale= 92.1281835, # updated 2021-07-28
                  discountRate.effectiveness = 0.035,
                  discountRate.costs = 0.035,
                  frailty = TRUE,
                  ##other_variance = 0,
                  grs_risk_threshold = 0.075,
                  cost_parameters =  c("Invitation" = 0                         # Invitation letter
                      + 0,                                      # Results letter
                      "Formal PSA" = 21 #from NICE guideline inflated to 2021 prices
                      + 0                                    # PSA analysis
                      + 0 * 0,                            # No GP primary care
                      "Formal panel" =  0                # test sampling, primary care
                      + 0                                    # PSA analysis not included in panel price
                      + 0                                    # From BergusMedical (official lab for Sthlm3)
                      + 0 * 0,                            # No GP for formal
                      "Opportunistic PSA" = 21               # PSA analysis
                      + 0 * 0,                          # GP primary care
                      "Opportunistic panel" = 0              # PSA analysis not included in panel price
                      + 0                                    # From BergusMedical (official lab for Sthlm3)
                      + 0 * 0,                          # GP primary care 
                      "Biopsy" = 581                           # Callendar 2021
                      + 0,                                # Pathology of biopsy #Callender 2021
                      "MRI" = 339,                             # MRI cost  #Callendar 2021
                      "Combined biopsy" = 581              # Biopsy cost (TBx), still called Combined biopsy
                      + 0,                                # Pathology of biopsy
                      "Assessment" = 545.03,                      # Callender 2021
                      "Prostatectomy" = 9808               # Robot assisted surgery #Callendar 2021
                      +  0*0*0                         # Radiation therapy
                      + 0*1,                                 # Urology and nurse visit
                      "Radiation therapy" = 6462*1          # Radiation therapy  #Callendar 2021
                      + 0*1                               # Oncologist new visit
                      + 0*1                               # Oncologist further visit
                      + 0*20                                  # Nurse visit
                      + 0*0.2,                           # Hormone therapy
                      "Active surveillance - yearly - w/o MRI" =
                      + 21*3                                # PSA sampling
                      + 105.19*2                                  # Urological appointment
                      + 581*0.33                               # Systematic biopsy
                      + 0*0.33,                           # Pathology of biopsy
                      "Active surveillance - yearly - with MRI" = #5052   # Urology visit and nurse visit
                      + 21*3                                # PSA sampling
                      + 105.19*2                                  # Urological appointment
                      + 339*0.33                               # MRI cost
                      + 581*1*0.33                           # Biopsy cost (SBx|TBx)
                      + 0*0.33,                           # Pathology of biopsy
                      "ADT+chemo" = 0,               # NEW: Chemo and hormone therapy #NICE model  LHRH treatment:Decapeptyl 11.25   injection (3-month dose)
                      "Post-Tx follow-up - yearly first" = 0 # Urologist and nurse consultation
                      + 0                                 # PSA test sampling
                      + 0,                                   # PSA analysis
                      "Post-Tx follow-up - yearly after" = 0  # PSA test sampling
                      + 0                                    # PSA analysis,
                      + 0,                                    # Telefollow-up by urologist
                      "Palliative therapy - yearly" = 7383, # Palliative care cost #Callendar 2021
                      "Terminal illness" = 7383/2,      # Terminal illness cost #
                      "Polygenic risk stratification" = 25),# Callender et al (2021) with exchange rate of approximately 12
                  background_utilities <-
                  data.frame(lower=c(0, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
                             upper=c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 1.0e55),
                             utility=c(1, 0.934, 0.922, 0.922, 0.905, 0.905, 0.849, 0.849, 0.804, 0.804, 0.785, 0.785,
                                 0.734, 0.734)),
                  ## Latest review 2019 based on Heijnsdijk 2012, Magnus 2019 and extended review by Shuang - PORPUS-U
                  utility_estimates = 
                                        #c("Invitation" = 1,                    # Heijnsdijk 2012
                                        #                   "Formal PSA" = 1,                 # Heijnsdijk 2012
                                        #                  "Formal panel" = 1,               # Heijnsdijk 2012
                                        #                 "Opportunistic PSA" = 1,          # Heijnsdijk 2012
                                        #                "Opportunistic panel" = 1,        # Heijnsdijk 2012
                                        #               "Biopsy" = 0.90,                     # Heijnsdijk 2012
                                        #              "Cancer diagnosis" = 0.80,           # Heijnsdijk 2012
                                        #             "Prostatectomy part 1" = 0.860,      # Magnus 2019 (Krahn 2009, Ku 2009)
                                        #            "Prostatectomy part 2" = 0.900,      # Magnus 2019 (Krahn 2009, Ku 2009)
                                        #           "Radiation therapy part 1" = 0.890,  # Krahn 2009
                                        #          "Radiation therapy part 2" = 0.920,  # Krahn 2009
                                        #         "Active surveillance" = 0.980,       # Loeb 2018
                                        #        "Postrecovery period" = 0.930,       # Magnus 2019 (Avila 2014, Bremner 2014, Krahn 2013, Ku 2009)
                                        #       "ADT+chemo" = 0.803,                 # Krahn 2003
                                        #      "Palliative therapy" = 0.680,        # *15D value; Magnus 2019 (Farrkila 2014, Torvinen 2013)
                                        #     "Terminal illness" = 0.40,           # Heijnsdijk 2012
                                        #    "Death" = 0.00),
                  ## Latest review 2019 based on Heijnsdijk 2012, Magnus 2019 and extended review by Shuang - EQ-5D
                  c("Invitation" = 1,                   # Heijnsdijk 2012
                    "Formal PSA" = 0.99,                 # Heijnsdijk 2012
                    "Formal panel" = 0.99,               # Heijnsdijk 2012
                    "Opportunistic PSA" = 0.99,          # Heijnsdijk 2012
                    "Opportunistic panel" = 0.99,        # Heijnsdijk 2012
                    "Biopsy" = 0.90,                     # Heijnsdijk 2012
                    "Combined biopsy" = 0.90,            # Heijnsdijk 2012
                    "Cancer diagnosis" = 0.80,           # Heijnsdijk 2012
                    "Prostatectomy part 1" = 0.829,      # Hall 2015
                    "Prostatectomy part 2" = 0.893,      # Extended review by SH (Glazener 2011, Korfarge 2005, Hall 2015)
                    "Radiation therapy part 1" = 0.818,  # Hall 2015
                    "Radiation therapy part 2" = 0.828,  # Extended review by SH (Korfarge 2005, Hall 2015)
                    "Active surveillance" = 0.9,         # Loeb 2018
                    "Postrecovery period" = 0.861,       # Extended review by SH (Torvinen 2013, Waston 2016)
                    "ADT+chemo" = 0.727,                 # Extended review by SH (Hall 2019, Diels 2015, Loriot 2015, Skaltsa 2014, Wu 2007, Chi 2018)
                    "Palliative therapy" = 0.62,         # Magnus 2019
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
                  pTBxG0ifG1_MRIpos=0.2474775,          # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) #Updated 2020-03-24 from Bx -> TBx
                  pTBxG0ifG2_MRIpos=0.06570613,          # Pr(TBx gives ISUP 0 | ISUP 2, MRI+) #Updated 2020-03-24 from Bx -> TBx
                  pTBxG1ifG2_MRIpos=0,          # Pr(TBx gives ISUP 1 | ISUP 2, MRI+) (not used) #Updated 2020-03-24 from Bx -> TBx
                  currency_rate = 1/1,    # Riksbanken 2018 #What is this?
                  risk_psa_threshold=1.5, # PSA threshold for risk-stratified screening
                  risk_lower_interval=6, # re-screening interval for lower risk
                  risk_upper_interval=4, # re-screening interval for higher risk
                                        #start_screening = 55.0, # start of organised screening
                                        #stop_screening = 69.0,  # end of organised screening
                  screening_interval = 2.0,
                  mu0=c(0.0039015,
                      0.000228,
                      0.0001295,
                      0.0000995,
                      0.0000825,
                      0.0000855,
                      0.000085,
                      0.0000655,
                      0.0000655,
                      0.000056,
                      0.000069,
                      0.0000755,
                      0.0000825,
                      0.0001035,
                      0.000111,
                      0.000143,
                      0.000187,
                      0.0002375,
                      0.0003135,
                      0.000324,
                      0.000349,
                      0.000362,
                      0.0003665,
                      0.0003635,
                      0.000387,
                      0.000426,
                      0.0004215,
                      0.0004565,
                      0.0005045,
                      0.000526,
                      0.0005705,
                      0.0006145,
                      0.000644,
                      0.0007075,
                      0.0007565,
                      0.0008275,
                      0.0008955,
                      0.0010465,
                      0.0009965,
                      0.0011255,
                      0.0012155,
                      0.001328,
                      0.0014455,
                      0.0015865,
                      0.0017045,
                      0.001886,
                      0.002026,
                      0.0021955,
                      0.002346,
                      0.002566,
                      0.002774,
                      0.002982,
                      0.003232,
                      0.003411,
                      0.003696,
                      0.003977,
                      0.0044655,
                      0.004836,
                      0.005312,
                      0.005773,
                      0.0063245,
                      0.0069025,
                      0.0077435,
                      0.008446,
                      0.0091055,
                      0.0100065,
                      0.0109515,
                      0.0119085,
                      0.013035,
                      0.0142925,
                      0.0153615,
                      0.0168065,
                      0.0187835,
                      0.0214235,
                      0.023645,
                      0.0264195,
                      0.029666,
                      0.033073,
                      0.037241,
                      0.04129,
                      0.0462225,
                      0.0518525,
                      0.057735,
                      0.0658325,
                      0.074345,
                      0.083563,
                      0.0949735,
                      0.1060235,
                      0.1198965,
                      0.1343965,
                      0.1469585,
                      0.164905,
                      0.1820145,
                      0.199694,
                      0.2212765,
                      0.244611,
                      0.2687395,
                      0.2855855,
                      0.308576,
                      0.339533,
                      0.3638745,
                      0.3638745,
                      0.3638745,
                      0.3638745,
                      0.3638745,
                      0.3638745
                      ),
                                        # 2017-2019 death, rates from UK life tables,
                  background_utilities =
                  data.frame(lower=c(0, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
                             upper=c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 1.0e55),
                             utility=c(1, 0.934, 0.922, 0.922, 0.905, 0.905, 0.849, 0.849, 0.804, 0.804, 0.785, 0.785,
                                 0.734, 0.734)), #Updated to UK
                  rescreening=uk_rescreening,
                  prtx =
                  data.frame(DxY=2016,
                             Age=c(50,50,50,55,55,55,60,60,60,65,65,65,70,70,70,75,75,75,80,80,80,85,85,85),
                             G=as.integer(c(6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8)-6),
                             CM=c(0.854149,0.854149,0.8641953,0.8641953,0.857022556,0.8641953,0.941532765,0.941532765,0.258915744,0.258915744,0.299382813,0.299382813,0.334049772,0.40071331,0.720422758,0.720422758,0.131734693,0.131734693,0.160002774,0.160002774,0.157891579,0.299838682,0.691923648,0.691923648),
                             RP=c(0.121796351,0.121796351,0.076248984,0.076248984,0.062706767,0.076248984,0.036925616,0.036925616,0.6028688,0.6028688,0.404516093,0.404516093,0.198613215,0.047024925,0.018454266,0.018454266,0.538923635,0.538923635,0.339997226,0.339997226,0.151101511,0.014589653,0.012621807,0.012621807),
                             RT=c(0.024054328,0.024054328,0.059555716,0.059555716,0.080270677,0.059555716,0.021541619,0.021541619,0.130861376,0.130861376,0.296101094,0.296101094,0.467337013,0.552261765,0.261122976,0.261122976,0.329341672,0.329341672,0.5,0.5,0.69100691,0.685571665,0.295454545,0.295454545))
                  ))
