EdnaParameters <- function(base = ShuangParameters)
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

                   eol=NULL # placemarker
                   
               ))
    
