EdnaParameters <- function(base = ShuangParameters())
    modifyList(base,
               list(

                   ## screening parameters
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

                   eol=1 # placemarker
               ))
