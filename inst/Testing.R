library(prostata)
## re-run with the old parameters to check if you get the same results


cost_parameters <- prostata:::FhcrcParameters$cost_parameters # copy
cost_parameters["Opportunistic PSA"] <- 350 # update
cost_parameters["Opportunistic panel"] <- 350 # update

sim1 <- callFhcrc(n=1e6, screen="regular_screening", mc.cores=2, 
                    parms=list(start_screening=55, stop_screening=69, screening_interval=4,
                               panel=1,
                               ## which costs are used for GP-based testing?
                               PSA_FP_threshold_nCa=4.21, # reduce FP in no cancers with PSA threshold
                               PSA_FP_threshold_GG6=3.76, # reduce FP in GG 6 with PSA threshold
                               PSA_FP_threshold_GG7plus=3, # reduce FP in GG >= 7 with PSA threshold
                               panelReflexThreshold = 1.0,
                               cost_parameters=cost_parameters))
## sim1.5
## sim2
sim0 <- callFhcrc(n=1e6, screen="noScreen", mc.cores=2)
simPSA <- callFhcrc(n=1e6, screen="regular_screening", mc.cores=2, 
                  parms=list(start_screening=55, stop_screening=69, screening_interval=4,
                             panel=0,
                             ## which costs are used for GP-based testing?
                             cost_parameters=cost_parameters))
## n=1e7? 1e8?

