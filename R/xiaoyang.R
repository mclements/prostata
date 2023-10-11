##==============================================================================
## FILENAME: Xiaoyang_Du_Study3.R
## PROJECT: Cost-effectiveness of Artificial Intelligence Assisted Prostate Pathology
## PURPOSE: Update input data as parameters for the simulation
## AUTHOR: Xiaoyang Du

## CREATED:	2023-03-27
## UPDATED: 2023-09-04

## INPUT DATA: 				
## OUTPUT: 	

## R VERSION: RStudio 2022.12.0+353 
##==============================================================================


XiaoyangParameters <- function(year=2022, Biopsy_cost=3010 + 4543.81, Pathology_cost=4543.81, pPathPath=1/3, pAICorenegG0=0.809,
                               pAICorenegG1=0.51, pAICorenegG2=0.383, pAICorenegG4plus=0.216, ...) {
    
    stopifnot(year %in% 2022)

    params <-
        list(
            Andreas = FALSE,
            ## Cost for 2022
            cost_parameters =  c("Opportunistic PSA" = 1975*0.2            # GP visit, Primary care
                                 + 61.54*1,                                # PSA test analysis
                                 ## "Opportunistic panel" = 61.54*1        # PSA analysis (not used)
                                 ## + 2300                                 # From BergusMedical (official lab for Sthlm3) (need to update the price to 2022)
                                 ## + 1975*0.2,                            # GP primary care
                                 "Assessment" = 1460,                      # Specialist and nurse consultation
                                 ## "MRI" = 4666,                          # MRI cost (no used)
                                 "Biopsy" = Biopsy_cost,                   # SBx (only used if no MRI)
                                 ## reminder: "Combined biopsy" should have been "Biopsy after MRI"
                                 ## "Combined biopsy" = 3010*1.5           # SBx+TBx cost (Combined = 1.5*SBx; default)(no used)
                                 ## + 4543.81,                             # MRI cost (no used)
                                 ## "Combined biopsy" = 2990 + 4543.81,    # TBx cost (yes, it's a stupid name:())(no used)
                                 "Pathology" = Pathology_cost,             # Pathology of biopsy 
                                 "AI pathology" = 10*10.6317,             # Assumed, 10 euro
                                 "Active surveillance - yearly - w/o MRI" = 1460*1  # Specialist and nurse consultation
                                 + 403*3                                   # PSA test sampling (3 times per year)
                                 + 61.54*3                                 # PSA test analysis (3 times per year)
                                 + 3010/3                                  # SBx cost (once every 2-3 year)
                                 + 4543.81/3,                              # Pathology of biopsy
                                 ## + 1000,                                   # AI Pathology
                                 ## "Active surveillance - yearly - with MRI" = 1460*1  # Specialist and nurse consultation (not used)
                                 ## + 403*3                                # PSA test sampling (3 times per year) (not used)
                                 ## + 61.54*3                              # PSA test analysis (3 times per year) (not used)
                                 ## + 4666/3                               # MRI cost (once every 2-3 years) (not used)
                                 ## + 3010*1.5/3                           # SBx+TBx cost (not used)
                                 ## + 2990/3,                              # TBx cost (not used)
                                 ## + 4543.81/3,                           # Pathology of biopsy (not used)
                                 ## + 1000,                                # AI Pathology (not used)
                                 "Prostatectomy" = 130653.32*1             # Robot assisted surgery
                                 + 7137.18*20*0.25                         # Radiation therapy (20 times during the whole course, 25% recurrence)
                                 + 1460*1,                                 # Specialist and nurse consultation
                                 "Radiation therapy" = 7137.18*20          # Radiation therapy (20 times during the whole course)
                                 + 3210*1                                  # Oncologist new visit
                                 + 1753*1                                  # Oncologist further visit
                                 + 400*20                                  # Nurse visit (20 times)
                                 + 76432.04*0.2,                           # Hormone therapy (20% patients)
                                 "ADT+chemo" = 160774.24,                  # Metastasis (hormone + chemo therapy)
                                 "Post-Tx follow-up - yearly first" = 1460 # Urologist and nurse consultation
                                 + 403                                     # PSA test sampling
                                 + 61.54,                                  # PSA test analysis
                                 "Post-Tx follow-up - yearly after" = 146  # Tele follow-up
                                 + 403                                     # PSA test sampling
                                 + 61.54,                                  # PSA test analysis
                                 "Palliative therapy - yearly" = 183002.65,            # Palliative care cost-yearly
                                 "Terminal illness" = 91501.32),           # Terminal illness cost-yearly
            
            production = data.frame(ages = c(0, 54, 64, 74),
                                    values=apply(
                                        rbind(0.889, 0.773, 0.192, 0), # employment proportion, full time 2022
                                        1,
                                        function(empl) empl*
                                                       30.0/40*                       # average working hours 2022
                                                       38300*12*                      # average salary 2022 
                                                       1.372)),                      # including non-optional social fees 2022
            ## average salary 2022
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0110__AM0110B/LonYrkeUtbildningA/table/tableViewLayout1/
            ## employment proportion 2022 
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0401__AM0401A/NAKUBefolkning2Ar/
            ## average working hours 2022
            ## https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__AM__AM0401__AM0401S/NAKUFaktMedArbtidAr/table/tableViewLayout1/
            
            lost_production_years= c(## "Formal PSA"=2/24/365.25,
                ## "Formal panel"=2/24/365.25,
                "Opportunistic PSA"=2/24/365.25,
                ## "Opportunistic panel"=2/24/365.25,
                ## "MRI"=2/24/365.25,
                "Biopsy"=2/24/365.25,
                ## "Combined biopsy"=2/24/365.25,
                "Assessment"=2/24/365.25,         
                "Prostatectomy"=6/52,
                "Radiation therapy"=8/52,
                "Active surveillance - yearly - w/o MRI"=
                    3 * 2/24/365.25                  # PSA tests
                + 1 * 2/24/365.25                  # Urologist visit
                + 0.33 * 2/24/365.25,              # Biopsy (SBx)
                ## "Active surveillance - yearly - with MRI"=
                ## 3 * 2/24/365.25                 # PSA tests
                ## + 1 * 2/24/365.25               # Urologist visit
                ## + 0.33 * 2/24/365.25            # MRI
                ## + 0.33 * 2/24/365.25,           # Combined biopsy (TBx|SBx)
                "Post-Tx follow-up - yearly"  =
                    2/24/365.25,                     # PSA tests (Same for first year and following years)
                "Premature mortality" = 0,         # Keep it following Shuang's code
                "Long-term sick leave" = 0.0768*67.52/365.25, # Keep it following Shuang's code: 7.68% employed PCa patients (50-64) have long-term sick leave (based on 2016 data)
                ## "Early retirement" = 0.00203*235.5/365.25, # 0.203% employed PCa patients (50-64) have early retirement (based on 2016 data)
                "Terminal illness" = 6/12),
            

            ## Based on Hao et al 2022 and its references
            utility_estimates = c(## "Invitation" = 1,                 # Heijnsdijk 2012
                ## "Formal PSA" = 0.99,              # Heijnsdijk 2012
                ## "Formal panel" = 0.99,            # Heijnsdijk 2012
                "Opportunistic PSA" = 0.99,          # Heijnsdijk 2012
                ## "Opportunistic panel" = 0.99,     # Heijnsdijk 2012
                "Biopsy" = 0.90,                     # Heijnsdijk 2012
                "Cancer diagnosis" = 0.80,           # Heijnsdijk 2012
                "Prostatectomy part 1" = 0.860,      # Magnus 2019 (Krahn 2009, Ku 2009)
                "Prostatectomy part 2" = 0.900,      # Magnus 2019 (Krahn 2009, Ku 2009)
                "Radiation therapy part 1" = 0.890,  # Krahn 2009
                "Radiation therapy part 2" = 0.920,  # Krahn 2009
                "Active surveillance" = 0.980,       # Loeb 2018
                "Postrecovery period" = 0.930,       # Magnus 2019 (Avila 2014, Bremner 2014, Krahn 2013, Ku 2009)
                "ADT+chemo" = 0.803,                 # ADT+chemo Krahn 2003
                "Palliative therapy" = 0.68,         # *15D value; Magnus 2019 (Farrkila 2014, Torvinen 2013)
                "Terminal illness" = 0.40,           # Heijnsdijk 2012
                "Death" = 0.00),
            ## Utility duration is given in years, based on Hao et al 2022
            utility_duration = c(## "Invitation" = 0.0,
                ## "Formal PSA" = 1/52,
                ## "Formal panel" = 1/52,
                "Opportunistic PSA" = 1/52,
                ## "Opportunistic panel" = 1/52,
                "Biopsy" = 3/52,
                ## "Combined biopsy" = 3/52,
                "Cancer diagnosis" = 1/12,
                "Prostatectomy part 1" = 2/12,
                "Prostatectomy part 2" = 10/12,
                "Radiation therapy part 1" = 2/12,
                "Radiation therapy part 2" = 10/12,
                "Active surveillance" = 7,
                "Postrecovery period" = 9,
                "ADT+chemo" = 1.5,                   
                "Palliative therapy" = 12/12,        
                "Terminal illness" = 6/12),
            
            pSBxG0ifG1=0.056600317,                  # Pr(SBx gives ISUP 0 | ISUP 1) STHML3MRI; Van Der Leest et al 2019
            pSBxG0ifG2=0.109053889,                  # Pr(SBx gives ISUP 0 | ISUP 2-3) STHML3MRI; Van Der Leest et al 2019
            pSBxG0ifG4plus=0.109053889,              # Pr(SBx gives ISUP 0 | ISUP 4+) STHML3MRI; Van Der Leest et al 2019
            ## pMRIposG0=0.184011278,                # Pr(MRI+ | ISUP 0) (not used) STHML3MRI; Van Der Leest et al 2019
            ## pMRIposG1=0.270423739,                # Pr(MRI+ | ISUP 1) (not used) STHML3MRI; Van Der Leest et al 2019
            ## pMRIposG2=0.922763676,                # Pr(MRI+ | ISUP 2+) (not used) STHML3MRI; Van Der Leest et al 2019
            ## PBxG0=0,                              # Pr(Bx+ | ISUP 0) (Bx+ refers to any biopsy, assumed SBx and TBx are accurate) (not used)
            ## pTBxG0ifG1_MRIpos=0,                  # Pr(SBx+TBx gives ISUP 0 | ISUP 1, MRI+) (TBx refers to either TBx or SBx/TBx) (not used)
            ## pTBxG0ifG2_MRIpos=0,                  # Pr(SBx+TBx gives ISUP 0 | ISUP 2+, MRI+) (not used)
            ## pTBxG0ifG1_MRIpos= 0.534883721,       # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) STHML3MRI; Van Der Leest et al 2019 (not used)
            ## pTBxG0ifG2_MRIpos= 0.063636364,       # Pr(TBx gives ISUP 0 | ISUP 2+, MRI+) STHML3MRI; Van Der Leest et al 2019 (not used)
            ## pHumanpos=1,                          # Pr(Human path+ | Cancer+) (Assumption for base case: human pathologist is 100% accurate)
            ## pHumanNeg=1,                          # Pr(Human path- | Cancer-)
            pAIposG1=0.992,                          # Pr(AI+ | ISUP 1, sensitivity = 0.99) Henrik's data output
            pAIposG2=1,                              # Pr(AI+ | ISUP 2-3, sensitivity = 0.99) Henrik's data output
            pAIposG4plus=1,                          # Pr(AI+ | ISUP 4+, sensitivity = 0.99) Henrik's data output
            pAICorenegG0=pAICorenegG0,                   # Pr(Core- | AI assisted, ISUP 0) (used for cost reduction calculation)
            pAICorenegG1=pAICorenegG1,                    # Pr(Core- | AI assisted, ISUP 1) (used for cost reduction calculation)
            pAICorenegG2=pAICorenegG2,                   # Pr(Core- | AI assisted, ISUP 2-3) (used for cost reduction calculation)
            pAICorenegG4plus=pAICorenegG4plus,               # Pr(Core- | AI assisted, ISUP 4+)  (used for cost reduction calculation)
            pPathPath=pPathPath,                        # Proportion of pathology cost due to the pathologist
            pReducedBxCostG0=pAICorenegG0*pPathPath*Pathology_cost/Biopsy_cost,         # Proportional reduction in cost, Henrik's data output
            pReducedBxCostG1=pAICorenegG1*pPathPath*Pathology_cost/Biopsy_cost,         # Proportional reduction in cost, Henrik's data output
            pReducedBxCostG2=pAICorenegG2*pPathPath*Pathology_cost/Biopsy_cost,         # Proportional reduction in cost, Henrik's data output
            pReducedBxCostG4plus=pAICorenegG4plus*pPathPath*Pathology_cost/Biopsy_cost, # Proportional reduction in cost, Henrik's data output
            ## pReducedBxCostG0=pAICorenegG0/3*4543.81/(4543.81+3010),       # Proportional reduction in cost, Henrik's data output
            ## pReducedBxCostG1=pAICorenegG1/3*4543.81/(4543.81+3010),        # Proportional reduction in cost, Henrik's data output
            ## pReducedBxCostG2=pAICorenegG2/3*4543.81/(4543.81+3010),       # Proportional reduction in cost, Henrik's data output
            ## pReducedBxCostG4plus=pAICorenegG4plus/3*4543.81/(4543.81+3010),   # Proportional reduction in cost, Henrik's data output
            currency_rate = 1/10.6317,               # Riksbanken 2022, Krona-EUR
            
            ## https://www.mortality.org/ 2020-2022, 1*5 Death rates for Swedish male
            mu0=c(0.002279, 0.000212, 0.000104, 0.000097, 0.000106, 0.000042, 0.000062, 0.000094, 0.000062, 
                  0.000073, 0.000057, 0.000087, 0.000077, 0.0001, 0.000117, 0.000156, 0.000245, 0.000321, 
                  0.000392, 0.000597, 0.00056, 0.000577, 0.000558, 0.000562, 0.000613, 0.000671, 0.00058,
                  0.000608, 0.000765, 0.000621, 0.000645, 0.000664, 0.000608, 0.000672, 0.000638, 0.000749,
                  0.000723, 0.000765, 0.00075, 0.000751, 0.000967, 0.000926, 0.001022, 0.001005, 0.001085,
                  0.001238, 0.001213, 0.001575, 0.001472, 0.001813, 0.001988, 0.002322, 0.002397, 0.002414,
                  0.003256, 0.003456, 0.003915, 0.004065, 0.004574, 0.005295, 0.006063, 0.006382, 0.007304,
                  0.007884, 0.009241, 0.010023, 0.01143, 0.012196, 0.013824, 0.015285, 0.016601, 0.018574,
                  0.019871, 0.022534, 0.024974, 0.027412, 0.030283, 0.035432, 0.038687, 0.04384, 0.051071,
                  0.055954, 0.067053, 0.074813, 0.085092, 0.100737, 0.110326, 0.1276, 0.14978, 0.163505, 
                  0.193686, 0.222016, 0.245564, 0.274162, 0.306668, 0.334743, 0.367593, 0.420501, 0.453341, 
                  0.479057, 0.606353, 0.5903, 0.614028, 0.681244, 0.751617, 0.788766, 0.613565, 0.370761, 0.149514))

    modifyList(params, list(...))

}


## background health state values using Teni 2021
XiaoyangTables <- list(background_utilities =
                           data.frame(lower=c(0, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
                                      upper=c(30, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 1.0e55),
                                      utility=c(1, 0.925, 0.938, 0.93, 0.924, 0.914, 0.91, 0.91, 0.915, 0.909, 0.892, 0.865, 0.831, 0.803, 0.751)))


