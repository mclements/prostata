#' prostata
#'
#' @section Introduction:
#'
#' The packages provides implementations of discrete event simulation in both R
#' and C++.
#'
#' @section Background:
#'
#' Prostata is a natural history model of prostate cancer model which extends
#' the model developed by Ruth Etzioni and colleagues at the
#' [[http://www.fredhutch.org][Fred Hutchinson Cancer Research Center
#' (FHCRC)]]. This is a well validated prostate cancer natural history model
#' developed within the
#' [[http://cisnet.cancer.gov/prostate/profiles.html][Cancer Intervention and
#' Surveillance Modeling Network (CISNET)]]. We here provide our extensions of
#' the model with extended states and finer calibration. The model is provided
#' as an R package, depending on our
#' [[https://github.com/mclements/microsimulation][microsimulation]] frame work.
#' We specifically developed this package for modelling the cost-effectiveness
#' of prostate cancer screening, where many (e.g. 10^7) men are followed from
#' birth to death.
#'
#' @section Running the simulation:
#'
#' There are a number of available testing scenarios. They determine testing
#' frequencies and re-testing intervals over calendar period and ages.
#'
#' @docType package
#' @name prostata-package
#' @aliases prostata
#' @author Mark Clements \email{mark.clements@ki.se
#' Andreas}
#' @references \url{https://github.com/mclements/prostata}
#' @references \url{https://github.com/mclements/microsimulation}
#' @seealso \code{\link{Rcpp}}
#' @useDynLib prostata, .registration=TRUE
#' @import microsimulation
#' @importFrom Rcpp evalCpp compileAttributes
NULL

if(.Platform$OS.type == "unix") {
    #' @importFrom utils packageName
    pkg <- packageName()
    ## http://r.789695.n4.nabble.com/How-to-construct-a-valid-seed-for-l-Ecuyer-s-method-with-given-Random-seed-td4656340.html

    #' @rdname set.user.Random.seed
    #' @param seed PARAM_DESCRIPTION
    #' @keywords internal
    #' @importFrom microsimulation set.user.Random.seed
    set.user.Random.seed <- function (seed)
        microsimulation::set.user.Random.seed(seed,pkg)

    #' @rdname next.user.Random.substream
    #' @keywords internal
    #' @importFrom microsimulation next.user.Random.substream
    next.user.Random.substream <- function ()
        microsimulation::next.user.Random.substream(pkg)

    #' @rdname user.Random.seed
    #' @keywords internal
    #' @importFrom microsimulation user.Random.seed
    user.Random.seed <- function()
        microsimulation::user.Random.seed(pkg)
}

## initial values for the FHCRC model
FhcrcParameters <- list(
    revised_natural_history=TRUE,
    ## panel=FALSE,
    grade.clinical.rate.high=0.3042454700,
    startFullUptake = 1932.0,
    fullUptakePortion = 0.9,
    yearlyUptakeIncrease = 0.03,
    startUptakeMixture = 1945.0,
    endUptakeMixture = 1960.0,
    screeningIntroduced = 1995.0,
    shapeA = 3.8,
    scaleA = 15.0,
    shapeT = 2.16,
    scaleT = 11.7,
    tau2 = 0.0829, # log PSA measurement variance  - normal
    susceptible = susceptible <- 1.0, # portion susceptible
    g0=0.0005 / susceptible, # onset parameter
    g3p = exp(-6.9353063), # T3+ parameter
    gm = exp(-6.4824485), # metastatic parameter
    gc=0.0015, # clinical diagnosis parameter
    thetac=19.1334, # clinical diagnosis parameter after metastatic
    mubeta0=-1.609, # mean of beta0, where beta0 is the log PSA intercept at age 35 years
    sebeta0=0.2384, # SE of beta0
    mubeta1=0.04463, # mean of beta1, where beta1 is the log PSA slope prior to cancer onset
    sebeta1=0.0430, # SE of beta1
    mubeta2=c(0.0397,0.1678, 0.0), # base::grade: the gleason specific change in slope of log PSA after cancer onset
    sebeta2=c(0.0913,0.3968, 0.0), # base::grade: variance of beta2
    rev_mubeta2=c(0.051, 0.129, 0.1678), # ext::grade: same as above for extended gleason grade (6-, 7, 8+)
    rev_sebeta2=c(0.064, 0.087, 0.3968), # ext::grade
    alpha7 = log(0.2), # log of the proportion gleason 7 at age 35
    beta7 = 0.0642775, # slope of log proportion of gleason 7
    alpha8 = log(0.002), # log of the proportion of gleason 8+ at age 35
    beta8 = 0.1864655, # slope of log proportion of gleason 8+
    RR_T3plus=2.0, # prostate cancer mortality rate ratio comparing T3+ with T1-T2, add lit reference
    ## mubeta2.scale=1.0, # cf. 2.1
    ## beta.rho=0.62,
    c_txlt_interaction = 1.0, # treatment lead time interaction, default to no effect
    c_baseline_specific = 1.0, # scale factor on survival, default to no effect
    c_benefit_value0 = 10, # (value -> reduction): 0.04 -> 10%; 0.18 -> 20%; 10 -> 28%, generalised stage-shift parameter
    sxbenefit = 1.0, # hazard rate ratio for screening benefit, defaults to no effect
    c_benefit_type = 0, # 0=stage-shift (=> c_benefit_value0=10), 1=lead-time based (=> c_benefit_value1=0.1)
    c_benefit_value1 = 0.1, # approximate increase in cure for each year of lead-time, used for the lead-time based screening effect
    RP_mortHR = 0.56, # mortality hazard ratio for surgery over watchful-waiting from SPCG-4 Bill-Axelson 2014
    screeningParticipation = 0.75, # probability of actually having the first PSA test
    rescreeningParticipation = 0.95, # probability of actually having the re-screening PSA tests
    biopsyCompliance = 0.856, # Formal biopsy compliance, Schröder ERSPC 2014
    biopsySensitivityTimeProportionT1T2 = 0.5272495, # time portion when T1-T2 cancers are sensitivity to biopsies (expit from calibration). The remaining part, starting at onset, is not detectable.
    studyParticipation = 50.0/260.0, # observed fraction of population who participated in STHLM3 study
    nLifeHistories = 10L, screen = 0L, ## integers
    psaThreshold = 3.0,
    psaThresholdBiopsyFollowUp = 4.0, # revised PSA threshold for negative follow-up screen
    biomarker_model = 0, # biomarker_model = 0 random, biomarker_model = 1 psa/risk based correction of FP
    PSA_FP_threshold_nCa=4.15, # reduce FP in no cancers with PSA threshold
    PSA_FP_threshold_GG6=3.41, # reduce FP in GG 6 with PSA threshold
    ## Natural history calibration
    rTPF=1.0,
    rFPF=0.6,
    c_low_grade_slope=-0.006,
    stockholmTreatment = TRUE,
    discountRate.effectiveness = 0.03,
    discountRate.costs = 0.03,
    full_report = 1.0,
    formal_costs = 0.0,
    formal_compliance = 0.0,
    start_screening = 50.0, # start of organised screening
    stop_screening = 70.0,  # end of organised screening
    screening_interval = 2.0, # screening interval for regular_screening and introduced_screening
    introduction_year = 2015.0, # year to start organised screening for introduced_screening  & cancel the opportunistic testing under stopped_screening
    mu0=c(0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05,
        0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012,
        0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448,
        0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651,
        0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077,
        0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04,
        0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217,
        0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491,
        0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914,
        0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693,
        0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182,
        0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776,
        0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436,
        0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768,
        0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538,
        0.782609),
    hr_locoregional = transform(expand.grid(age=c(50,60,70), ext_grade=0:2,
                                            psa10=0:1),
                                hr = c(0.2009515, 0.4671682,
                                       1.0793264, 2.0419524,
                                       1.4795704, 1.8301608,
                                       1.7294489, 1.0701277,
                                       1.4378105, 0.7065310,
                                       1.0955767, 1.2093063,
                                       5.0236260, 3.2827077,
                                       2.1218456, 1.1981716,
                                       0.9827433, 0.8948585)),
    hr_metastatic = data.frame(age = c(50, 60, 70),
                               hr = c(0.9911762, 0.8275815, 0.7313553)),
    biopsy_sensitivity = data.frame(Year = c(1987, 1988, 1989, 1990, 1991, 1992,
                                             1993, 1994, 1995, 1996, 1997, 1998,
                                             1999, 2000),
                                    Sensitivity = c(0.502, 0.574, 0.642, 0.702,
                                                    0.749, 0.782, 0.803, 0.814,
                                                    0.821, 0.830, 0.844, 0.866,
                                                    0.894, 0.925)),
    cost_parameters = c("Invitation" = 50,
                        "Formal PSA" = 130,
                        "Formal panel" = 730,
                        "Opportunistic PSA" = 1910,
                        "Opportunistic panel" = 2510, #N.B. This one is new and should be used
                        "Biopsy" = 12348,
                        "Prostatectomy" = 117171,
                        "Radiation therapy" = 117171,
                        "Active surveillance" = 141358,
                        "Cancer death" = 585054),
    ## IHE doesn't use the postrecovery period (as reported in the Heijnsdijk 2012 reference), should we?
    production = data.frame(ages = c(0, 55, 65, 75),
                            values=c(467433.137375286, 369392.309986899, 45759.6141748681, 0.0)),
    lost_production_proportions= c("Formal PSA"=0.0011,
                                   "Formal panel"=0.0011,
                                   "Opportunistic PSA"=0.0025,
                                   "Opportunistic panel"=0.0025,
                                   "Biopsy"=0.0044,
                                   "Prostatectomy"=0.1083,
                                   "Radiation therapy"=0.1250,
                                   "Active surveillance"=0.0833,
                                   "Metastatic cancer"=0.7602),
    utility_estimates = c("Invitation" = 1,
                          "Formal PSA" = 0.99,
                          "Formal panel" = 0.99,
                          "Opportunistic PSA" = 0.99,
                          "Opportunistic panel" = 0.99,
                          "Biopsy" = 0.90,
                          "Prostatectomy part 1" = 0.67,
                          "Prostatectomy part 2" = 0.77,
                          "Radiation therapy part 1" = 0.73,
                          "Radiation therapy part 2" = 0.78,
                          "Active surveillance" = 0.97,
                          "Palliative therapy" = 0.60,
                          "Terminal illness" = 0.40,
                          "Metastatic cancer" = 0.85,
                          "Death" = 0.00),
    ## Utility duration is given in years.
    utility_duration = c("Invitation" = 0.0,
                         "Formal PSA" = 1/52,
                         "Formal panel" = 1/52,
                         "Opportunistic PSA" = 1/52,
                         "Opportunistic panel" = 1/52,
                         "Biopsy" = 3/52,
                         "Prostatectomy part 1" = 2/12,
                         "Prostatectomy part 2" = 10/12,
                         "Radiation therapy part 1" = 2/12,
                         "Radiation therapy part 2" = 10/12,
                         "Active surveillance" = 7,
                         "Palliative therapy" = 30/12,
                         "Terminal illness" = 6/12)
)
IHE <- list(prtx=data.frame(Age=50.0,DxY=1973.0,G=1:2,CM=0.6,RP=0.26,RT=0.14)) ## assumed constant across ages and periods
ParameterNV <- FhcrcParameters[sapply(FhcrcParameters,class)=="numeric" & sapply(FhcrcParameters,length)==1]
## ParameterIV <- FhcrcParameters[sapply(FhcrcParameters,class)=="integer" & sapply(FhcrcParameters,length)==1]

#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 15 rows and 3 variables:
#' \describe{
#'   \item{\code{psa}}{double COLUMN_DESCRIPTION}
#'   \item{\code{age}}{double COLUMN_DESCRIPTION}
#'   \item{\code{compliance}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"swedenOpportunisticBiopsyCompliance"
swedenOpportunisticBiopsyCompliance <- data.frame(
    psa = c(3, 5, 10, 3, 5, 10, 3, 5, 10, 3, 5, 10, 3, 5, 10),
    age = c(40, 40, 40, 50, 50, 50, 60, 60, 60, 70, 70, 70, 80, 80, 80),
    compliance = c(0.3764045, 0.5680751, 0.7727273, 0.3110770, 0.5726548, 0.7537372, 0.2385155, 0.4814588, 0.6929770, 0.1754264, 0.3685056, 0.5602030, 0.1629213, 0.2697368, 0.5010052))

#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 15 rows and 3 variables:
#' \describe{
#'   \item{\code{psa}}{double COLUMN_DESCRIPTION}
#'   \item{\code{age}}{double COLUMN_DESCRIPTION}
#'   \item{\code{compliance}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"swedenFormalBiopsyCompliance"
swedenFormalBiopsyCompliance <- cbind(expand.grid(psa=c(3,5,10),age=seq(40,80,10)),
                                compliance=FhcrcParameters$biopsyCompliance)
#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 24 rows and 6 variables:
#' \describe{
#'   \item{\code{DxY}}{double COLUMN_DESCRIPTION}
#'   \item{\code{Age}}{double COLUMN_DESCRIPTION}
#'   \item{\code{G}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{CM}}{double COLUMN_DESCRIPTION}
#'   \item{\code{RP}}{double COLUMN_DESCRIPTION}
#'   \item{\code{RT}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"stockholmTreatment"
stockholmTreatment <-
    data.frame(DxY=2008,
               Age=c(50,50,50,55,55,55,60,60,60,65,65,65,70,70,70,75,75,75,80,80,80,85,85,85),
               G=as.integer(c(6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8)-6),
               CM=c(0.231023,0.044872,0.25,0.333333,0.09542,0.328358,0.409439,0.12825,0.348101,0.479167,0.178182,0.401639,0.689013,0.359143,0.56087,0.876543,0.744444,0.809524,1,0.970711,0.952096,0.9375,1,1),
               RP=c(0.700623,0.815839,0.5,0.553592,0.740111,0.326226,0.49981,0.6866,0.291437,0.409879,0.552483,0.279235,0.210949,0.318374,0.179363,0.041152,0.058652,0.049689,0,0.009763,0.023952,0.0625,0,0),
               RT=c(0.068354,0.13929,0.25,0.113074,0.164469,0.345416,0.090751,0.185151,0.360462,0.110954,0.269335,0.319126,0.100038,0.322482,0.259767,0.082305,0.196903,0.140787,0,0.019526,0.023952,0,0,0))

secularTrendTreatment2008OR <-
    data.frame(year=as.numeric(1988:2009),
               OR=c(0.194409,0.155874,0.130704,0.154104,0.148039,0.237752,0.299491,0.330934,0.407047,0.377968,0.443079,0.551047,0.686101,0.700657,0.846134,0.93956,1.03026,1.191726,1.067736,1.012136,1,0.911471))


#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 52 rows and 5 variables:
#' \describe{
#'   \item{\code{age5}}{double COLUMN_DESCRIPTION}
#'   \item{\code{total_cat}}{double COLUMN_DESCRIPTION}
#'   \item{\code{shape}}{double COLUMN_DESCRIPTION}
#'   \item{\code{cure}}{double COLUMN_DESCRIPTION}
#'   \item{\code{scale}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"rescreening"
rescreening <- data.frame(age5 = c(30, 30, 30, 30, 35, 35, 35, 35, 40, 40,
40, 40, 45, 45, 45, 45, 50, 50, 50, 50, 55, 55, 55, 55, 60, 60,
60, 60, 65, 65, 65, 65, 70, 70, 70, 70, 75, 75, 75, 75, 80, 80,
80, 80, 85, 85, 85, 85, 90, 90, 90, 90), total_cat = c(0, 1,
3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0,
1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10,
0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10), shape = c(1.04293463521838,
0.825156738737092, 0.899003822110993, 0.589358680965725, 1.05948319345651,
0.895913335495644, 0.84619233996788, 0.636706369588122, 1.16028493785731,
0.97919993225436, 0.68848363112076, 0.667566641635962, 1.25727757690268,
1.11021543939268, 0.717059585208267, 0.724574127291869, 1.3247305288869,
1.18510576695807, 0.858171614262325, 0.749998577524088, 1.3340371479645,
1.20792880820659, 0.920717960779785, 0.852052206426261, 1.32268361303857,
1.24450960892428, 0.988946105418831, 0.919376727690091, 1.31156083562523,
1.26972967056357, 1.03454490825271, 0.954911145157684, 1.27126976718353,
1.26045834296745, 1.06684366805728, 0.968088650561122, 1.18494141305856,
1.20593215449012, 1.05998296831204, 0.989315997334925, 1.12395322065035,
1.13733607943079, 1.04483735235382, 0.988385269358649, 1.1220716456468,
1.09128664878198, 1.01326552242001, 0.948730052972523, 1.08315738475919,
1.09407460180658, 1.05197207574662, 0.978397483932468), cure = c(0.585284664126963,
0.545117838495959, 0.402804166722055, 0.14283955595853, 0.398430641565923,
0.327274451509188, 0.148468826128857, 0.00333583482942569, 0.257384245965907,
0.208165119313289, 0.157385208407254, 0.0696401481749125, 0.166635358789789,
0.133819199845029, 0.0610524867384188, 0.0535645284203467, 0.134003946557746,
0.114600739520109, 0.0537234127141967, 0.0308512461181499, 0.104553332825899,
0.0862550079726627, 0.0399447271614142, 0.0259518317972345, 0.10432761596883,
0.0823122881511086, 0.0355479853523502, 0.0190476453614591, 0.103738700348739,
0.0816458941897609, 0.0412696810847696, 0.0319927106467962, 0.108815314898247,
0.0867189417281315, 0.0472560235806836, 0.0300347697835532, 0.124252895869383,
0.111313062025094, 0.0675699701502124, 0.0373036736422192, 0.15590486866906,
0.153015363026187, 0.0974421665330981, 0.0563274705516883, 0.228204960789236,
0.2051853623161, 0.169506663277659, 0.0887344403643672, 0.334252927026757,
0.295865002768107, 0.24644310712191, 0.151125704749819), scale = c(3.33440594470802,
3.56052439703143, 0.247924244305486, 0.150095610340916, 4.05027295845979,
3.78800058764034, 0.35102275649495, 0.486681740716857, 3.62361340129802,
3.48003996256557, 0.661105185961319, 0.247063080430436, 2.80272048406229,
2.6202749242671, 0.764421647692264, 0.268131101357167, 2.29273349177322,
2.11664063949949, 0.749125964202465, 0.417294230771058, 2.06795900869592,
1.80352742984593, 0.78399749774511, 0.451625737966194, 1.82650979495336,
1.62873839529095, 0.806784448546511, 0.478644431229979, 1.61927028533258,
1.46189762327418, 0.810126410765397, 0.500531738601558, 1.53243134315106,
1.42674067036488, 0.869134221688119, 0.569591312385404, 1.44388696689722,
1.41174931932889, 0.938154963786188, 0.631891327069348, 1.47052388447475,
1.4543783566239, 1.01084363142901, 0.644940081254986, 1.29365253055963,
1.4079170674224, 1.03720449704243, 0.643101192871478, 1.08541958643012,
1.29524623033074, 1.02176143186057, 0.673175882772333))
#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 136 rows and 2 variables:
#' \describe{
#'   \item{\code{cohort}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{pop}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"pop1"
pop1 <- data.frame(cohort=2035:1900,
                    pop=c(rep(17239,32), 16854, 16085, 15504, 15604, 16381, 16705,
                       16762, 16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234,
                       13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 12435,
                       12955, 13013, 13125, 13065, 12249, 11103, 9637, 9009, 8828,
                       8350, 7677, 7444, 7175, 6582, 6573, 6691, 6651, 6641, 6268,
                       6691, 6511, 6857, 7304, 7308, 7859, 7277, 8323, 8561, 7173,
                       6942, 7128, 6819, 5037, 6798, rep(6567,46)))

#' @title List tables for the simulation
#' @description DATASET_DESCRIPTION
#' \strong{all_cause_mortality}
#' @format A data frame with 12120 rows and 4 variables:
#' \describe{
#'   \item{\code{BirthCohort}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Age}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{AnnualMortality}}{double COLUMN_DESCRIPTION}
#'   \item{\code{Survival}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{biopsy_frequency}
#' @format A data frame with 3 rows and 7 variables:
#' \describe{
#'   \item{\code{PSA.beg}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{PSA.end}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{X55.59}}{double COLUMN_DESCRIPTION}
#'   \item{\code{X60.64}}{double COLUMN_DESCRIPTION}
#'   \item{\code{X65.69}}{double COLUMN_DESCRIPTION}
#'   \item{\code{X70.74}}{double COLUMN_DESCRIPTION}
#'   \item{\code{X75.79}}{double COLUMN_DESCRIPTION}
#'}
#'
#' \strong{biopsy_sensitivity}
#' @format A data frame with 14 rows and 2 variables:
#' \describe{
#'   \item{\code{Year}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Sensitivity}}{double COLUMN_DESCRIPTION}
#'}
#'
#' \strong{dre}
#' @format A data frame with 4 rows and 4 variables:
#' \describe{
#'   \item{\code{psa.low}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{psa.high}}{double COLUMN_DESCRIPTION}
#'   \item{\code{sensitivity}}{double COLUMN_DESCRIPTION}
#'   \item{\code{specificity}}{double COLUMN_DESCRIPTION}
#'}
#'
#' \strong{prob_grade7}
#' @format A data frame with 51 rows and 2 variables:
#' \describe{
#'   \item{\code{slope}}{double COLUMN_DESCRIPTION}
#'   \item{\code{Pr.grade.7.}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{pradt}
#' @format A data frame with 8208 rows and 6 variables:
#' \describe{
#'   \item{\code{Tx}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Age}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{DxY}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Grade}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{NoADT}}{double COLUMN_DESCRIPTION}
#'   \item{\code{ADT}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{prtx}
#' @format A data frame with 2664 rows and 6 variables:
#' \describe{
#'   \item{\code{Age}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{DxY}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{G}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{CM}}{double COLUMN_DESCRIPTION}
#'   \item{\code{RP}}{double COLUMN_DESCRIPTION}
#'   \item{\code{RT}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{survival_dist}
#' @format A data frame with 42 rows and 3 variables:
#' \describe{
#'   \item{\code{Grade}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Time}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Survival}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{survival_local}
#' @format A data frame with 294 rows and 5 variables:
#' \describe{
#'   \item{\code{AgeLow}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{AgeHigh}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Grade}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Time}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Survival}}{double COLUMN_DESCRIPTION}
#'}
#' \strong{prob_earlystage}
#' @format A data frame with 84 rows and 6 variables:
#' \describe{
#'   \item{\code{BeginAge}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{EndAge}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{BeginPSA}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{EndPSA}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Grade}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Prob}}{double COLUMN_DESCRIPTION}
#'}
#' @export
"fhcrcData"
fhcrcData$biopsyOpportunisticComplianceTable <- swedenOpportunisticBiopsyCompliance
fhcrcData$biopsyFormalComplianceTable <- swedenFormalBiopsyCompliance
fhcrcData$secularTrendTreatment2008OR <- secularTrendTreatment2008OR
## https://www.socialstyrelsen.se/Lists/Artikelkatalog/Attachments/20008/2015-12-26.pdf
#' @title DATASET_TITLE
#' @description DATASET_DESCRIPTION
#' @format A data frame with 18 rows and 3 variables:
#' \describe{
#'   \item{\code{Age}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Sweden2000}}{double COLUMN_DESCRIPTION}
#'   \item{\code{World}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FUNCTION_NAME
#' @export
"ageStandards"
ageStandards <- data.frame(Age = cut(seq(0, 85, 5),
                                    breaks = c(seq(0,85,5), Inf),
                                    right = FALSE),
                          Sweden2000 = c(5.3, 6.9, 6.4, 5.7, 5.9, 6.7, 7.2, 6.9,
                                         6.6, 6.6, 7.6, 6.3, 4.9, 4.3, 4.1, 3.9,
                                         2.6, 2.3),
                          World = c(12.0, 10.0, 9.0, 9.0, 8.0, 8.0, 6.0, 6.0, 6.0,
                                    6.0, 5.0, 4.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.5))

#' @title Define and run simulation
#' @description Function to specify and run the microsimulation. A large number
#'     of simulated men, \code{n}, will cause the simulation to take longer time
#'     but will reduce the stochastic variation particularly for rare events.
#' @param n Integer number of men to simulate. Default: 10
#' @param screen String with one of the following simulated screening scenarios:
#'    \describe{
#'      \item{\code{noScreening}}{no screening test, only diagnosis from symptoms}
#'      \item{\code{randomScreen50to70}}{TBA}
#'      \item{\code{twoYearlyScreen50to70}}{two-yearly screening from age 50 to 70}
#'      \item{\code{fourYearlyScreen50to70}}{four-yearly screening from age 50 to 70}
#'      \item{\code{screen50}}{one screen at age 50}
#'      \item{\code{screen60}}{one screen at age 60}
#'      \item{\code{screen70}}{one screen at age 70}
#'      \item{\code{screenUptake}}{current testing pattern in Sweden}
#'      \item{\code{stockholm3_goteborg}}{TBA}
#'      \item{\code{stockholm3_risk_stratified}}{TBA}
#'      \item{\code{goteborg}}{risk stratified re-screening 2+4 from age 50 to 70}
#'      \item{\code{risk_stratified}}{risk stratified re-screening 4+8 from age 50 to 70}
#'      \item{\code{mixed_screening}}{risk stratified re-screening 2+4 from age 50 to 70 & opportunistic testing for other ages}
#'      \item{\code{regular_screen}}{TBA}
#'      \item{\code{single_screen}}{TBA}
#'      \item{\code{introduced_screening}}{TBA}
#'      \item{\code{stopped_screening}}{TBA}
#'    } . Default: 'noScreening'
#'
#' @param nLifeHistories Integer with number of men for all events should be
#'     recorded, Default: 10
#' @param seed Integer random number seed, Default: 12345
#' @param panel Boolean to use the Stockholm3 biomarker panel test
#'     characteristics, otherwise PSA is used, Default: FALSE
#' @param includePSArecords Boolean for extra reporting at the time of a
#'     screening test, Default: FALSE
#' @param includeDiagnoses Boolean for extra reporting at the time of diagnosis,
#'     Default: FALSE
#' @param flatPop Boolean to set all birth cohorts of equal size, Default: FALSE
#' @param pop Data.frame with two integer columns \code{cohort} with year of the
#'     birth cohorts and \code{pop} with the size of the corresponding birth
#'     cohorts, Default: pop1
#' @param tables List of data.frames to update the tables in fhcrcData, Default:
#'     IHE
#' @param debug Boolean to print debugging, Default: FALSE
#' @param parms List to update FhcrcParameters, Default: NULL
#' @param mc.cores Integer with the number of cores to use for the computation,
#'     Default: 1
#' @param print.timing Boolean should the required time be printed after the
#'     simulation run, Default: TRUE
#' @param ... TBA
#' @return A fhcrc object
#' @details TBA
#' @examples
#' \dontrun{
#' if(interactive()){
#'  library(prostata)
#'  sim1 <- callFhcrc(n=1e6, screen="screenUptake", mc.cores=3)
#'  }
#' }
#' @seealso
#'  \code{\link[parallel]{mclapply}}
#' @rdname callFhcrc
#' @export
#' @importFrom parallel mclapply
callFhcrc <- function(n=10, screen= "noScreening", nLifeHistories=10,
                      seed=12345, panel=FALSE, includePSArecords=FALSE,
                      includeDiagnoses=FALSE, flatPop = FALSE, pop = pop1,
                      tables = IHE, debug=FALSE, parms = NULL, mc.cores = 1,
                      print.timing = TRUE,...) {
  ## save the random number state for resetting later
  state <- RNGstate(); on.exit(state$reset())
  ## yes, we use the user-defined RNG
  RNGkind("user")
  set.user.Random.seed(seed)
  ## birth cohorts that should give approximately the number of men alive in Stockholm in 2012
  ## check the input arguments
  screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70",
               "fourYearlyScreen50to70", "screen50", "screen60", "screen70",
               "screenUptake", "stockholm3_goteborg",
               "stockholm3_risk_stratified", "goteborg", "risk_stratified",
               "mixed_screening","regular_screen", "single_screen",
               "introduced_screening", "stopped_screening")
  screen <- match.arg(screen, screenT)
  stopifnot(is.na(n) || is.integer(as.integer(n)))
  stopifnot(is.integer(as.integer(nLifeHistories)))
  ## these enum strings should be moved to C++
  stateT <- c("Healthy","Localised","Metastatic")
  ext_stateT <- c("Healthy","T1_T2","T3plus","Metastatic")
  gradeT <- c("Gleason_le_6","Gleason_7","Gleason_ge_8","Healthy")
  eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis", "toCancerDeath",
            "toOtherDeath","toScreen","toBiopsyFollowUpScreen",
            "toScreenInitiatedBiopsy","toClinicalDiagnosticBiopsy",
            "toScreenDiagnosis", "toOverDiagnosis", "toOrganised","toTreatment",
            "toCM","toRP", "toRT","toADT","toUtilityChange","toUtilityRemove",
            "toSTHLM3", "toOpportunistic","toT3plus", "toCancelScreens")
  diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
  treatmentT <- c("no_treatment","CM","RP","RT")
  psaT <- c("PSA<3","PSA>=3") # not sure where to put this...
  screenIndex <- which(screen == screenT) - 1
  timingfunction <- if (print.timing) function(x) print(system.time(x)) else function(x) x
  ## NB: sample() calls the random number generator (!)
  if (is.vector(pop)) {
      flatPop <- TRUE
      pop <- data.frame(cohort=pop,pop=1)
  }
  if (is.na(n)) {
    cohort <- pop$cohort[rep.int(1:nrow(pop),times=pop$pop)]
    n <- length(cohort)
  } else {
      if (flatPop) {
          cohort <- rep(pop$cohort,each=ceiling(n/nrow(pop))) #Need ceiling so int n=!0
          n <- ceiling(n/nrow(pop)) * nrow(pop) #To get the chunks right
      } else
          cohort <- sample(pop$cohort,n,prob=pop$pop/sum(pop$pop),replace=TRUE)
  }
  cohort <- sort(cohort)
  ## now separate the data into chunks
  chunks <- tapply(cohort, sort((0:(n-1)) %% mc.cores), I)
  ## set the initial random numbers
  currentSeed <- user.Random.seed()
  powerFun <- function(obj,FUN,n,...) {
    for(i in 1:n)
      obj <- FUN(obj,...)
    obj
  }
  initialSeeds <- Reduce(function(seed,i) powerFun(seed,parallel::nextRNGStream,10),
                         1:mc.cores, currentSeed, accumulate=TRUE)[-1]
  ns <- cumsum(sapply(chunks,length))
  ns <- c(0,ns[-length(ns)])
  ## Minor changes to fhcrcData
  fhcrcData$biopsyOpportunisticComplianceTable <- swedenOpportunisticBiopsyCompliance
  fhcrcData$biopsyFormalComplianceTable <- swedenFormalBiopsyCompliance
  fhcrcData$secularTrendTreatment2008OR <- secularTrendTreatment2008OR
  if (!is.null(tables))
      for (name  in names(tables))
          fhcrcData[[name]] <- tables[[name]]
  fhcrcData$rescreening <- rescreening
  fhcrcData$rescreening$total <- fhcrcData$rescreening$total_cat
  fhcrcData$prtx$Age <- as.double(fhcrcData$prtx$Age)
  fhcrcData$prtx$DxY <- as.double(fhcrcData$prtx$DxY)
  fhcrcData$prtx$G <- fhcrcData$prtx$G - 1L
  fhcrcData$pradt$Grade <- fhcrcData$pradt$Grade - 1L
  fhcrcData$biopsy_sensitivity$Year <- as.double(fhcrcData$biopsy_sensitivity$Year)
  fhcrcData$pradt$Age <- as.double(fhcrcData$pradt$Age)
  fhcrcData$pradt$DxY <- as.double(fhcrcData$pradt$DxY)
  ## fhcrcData$biopsyComplianceTable <-
  ##     data.frame(expand.grid(psa=c(4,7,10),age=seq(55,75,by=5)),
  ##                compliance=unlist(fhcrcData$biopsy_frequency[,-(1:2),]))
  fhcrcData$survival_local <-
      with(fhcrcData$survival_local,
           data.frame(Age=as.double(AgeLow),Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  fhcrcData$survival_dist <-
      with(fhcrcData$survival_dist,
           data.frame(Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  updateParameters <- parms
  updateParameters$nLifeHistories <- as.integer(nLifeHistories)
  updateParameters$screen <- as.integer(screenIndex)
  parameter <- FhcrcParameters
  for (name in names(updateParameters)){
      if(!(name %in% names(parameter)))
          warning("Name in parms argument not in FhcrcParameters: ",name,".")
      parameter[[name]] <- updateParameters[[name]]
  }
  pind <- sapply(parameter,class)=="numeric" & sapply(parameter,length)==1
  bInd <- sapply(parameter,class)=="logical" & sapply(parameter,length)==1
  if (parameter$stockholmTreatment)
      fhcrcData$prtx <- stockholmTreatment
  ## check some parameters for sanity
  if (panel && parameter["rTPF"]>1) stop("Panel: rTPF>1 (not currently implemented)")
  if (panel && parameter["rFPF"]>1) stop("Panel: rFPF>1 (not currently implemented)")
  ## now run the chunks separately
  timingfunction(out <- parallel::mclapply(1:mc.cores,
                function(i) {
                  chunk <- chunks[[i]]
                  set.user.Random.seed(initialSeeds[[i]])
                  .Call("callFhcrc",
                        parms=list(n=as.integer(length(chunk)),
                            firstId=ns[i],
                            panel=panel, # bool
                            debug=debug, # bool
                            cohort=as.double(chunk),
                            parameter=unlist(parameter[pind]),
                            bparameter=unlist(parameter[bInd]),
                            otherParameters=parameter[!pind & !bInd],
                            tables=fhcrcData,
                            includePSArecords=includePSArecords,
                            includeDiagnoses=includeDiagnoses),
                        PACKAGE="prostata")
                }, mc.cores = mc.cores))
  ## Apologies: we now need to massage the chunks from C++
  ## reader <- function(obj) {
  ##   out <- cbind(data.frame(state=enum(obj$state[[1]],stateT),
  ##                           dx=enum(obj$state[[2]],diagnosisT),
  ##                           psa=enum(obj$state[[3]],psaT),
  ##                           cohort=obj$state[[4]]),
  ##                data.frame(obj[-1]))
  ##   out$year <- out$cohort + out$age
  ##   out
  ## }
  ext_state2state <- function(obj) # collapses T-stages to localised
      `levels<-`(factor(obj),list(Healthy="Healthy",Localised=list("T1_T2","T3plus"),Metastatic="Metastatic"))
  cbindList <- function(obj) # recursive
    if (is.list(obj)) do.call("cbind",lapply(obj,cbindList)) else data.frame(obj)
  rbindList <- function(obj) # recursive
      if (is.list(obj)) do.call("rbind",lapply(obj,rbindList)) else data.frame(obj)
  rbindExtract <- function(obj,name)
      do.call("rbind",lapply(obj, function(obji) data.frame(obji[[name]])))
  reader <- function(obj) {
    obj <- cbindList(obj)
    out <- cbind(data.frame(state=ext_state2state(enum(obj[[1]],ext_stateT)),
                            ext_state=enum(obj[[1]],ext_stateT),
                            grade=enum(obj[[2]],gradeT),
                            dx=enum(obj[[3]],diagnosisT),
                            psa=enum(obj[[4]],psaT),
                            cohort=obj[[5]]),
                 data.frame(obj[,-(1:5)]))
    out
  }
  ## grab all of the pt, prev, ut, events from summary
  ## pt <- lapply(out, function(obj) obj$summary$pt)
  if (length(out[[1]]$summary) == 0) summary <- list()
  else {
      summary <- lapply(seq_along(out[[1]]$summary),
                        function(i) do.call("rbind",
                                            lapply(out, function(obj) reader(obj$summary[[i]]))))
      names(summary) <- names(out[[1]]$summary)
      states <- c("state","ext_state","grade","dx","psa","cohort")
      names(summary$prev) <- c(states,"age","count")
      names(summary$pt) <- c(states,"age","pt")
      names(summary$ut) <- c(states,"age","ut")
      names(summary$events) <- c(states,"event","age","n")
      if(FALSE) age <- NULL # To pass false-positive check note
      summary <- lapply(summary,function(obj) within(obj,year <- cohort+age))
      enum(summary$events$event) <- eventT
}

  ## lifeHistories <- do.call("rbind",lapply(out,function(obj) data.frame(obj$lifeHistories)))
  ## psarecord <- do.call("rbind",lapply(out,function(obj) data.frame(obj$psarecord)))
  ## diagnoses <- do.call("rbind",lapply(out,function(obj) data.frame(obj$diagnoses)))
  ## falsePositives <- do.call("rbind",lapply(out,function(obj) data.frame(obj$falsePositives)))
  ## parameters <- do.call("rbind",lapply(out,function(obj) data.frame(obj$parameters)))
  lifeHistories <- rbindExtract(out,"lifeHistories")
  psarecord <- rbindExtract(out,"psarecord")
  diagnoses <- rbindExtract(out,"diagnoses")
  falsePositives <- rbindExtract(out,"falsePositives")
  parameters <- rbindExtract(out,"parameters")

  appendMeans <- function(x) c(x,
                              mean.sum = x[["sum"]] / x[["n"]],
                              mean.sumsq = x[["sumsq"]] / x[["n"]])
  natural.history.summary <- data.frame(tmc_minus_t0 = appendMeans(sapply(rbindExtract(out,"tmc_minus_t0"), sum)))

  ## Identifying elements without name which also need to be rbind:ed
  societal.costs <- do.call("rbind",lapply(out,function(obj) data.frame(obj$costs))) #split in sociatal and healthcare perspective
  ## names(costs) <- c("type","item","cohort","age","costs")
  names(societal.costs) <- c("type","item","age","costs")
  societal.costs$type <- factor(ifelse(societal.costs$type,
                                       "Productivity loss",
                                       "Health sector cost")) # societal perspective
  healthsector.costs <- societal.costs[societal.costs["type"] == "Health sector cost", c("item", "age", "costs")] # healthcare perspective
  names(lifeHistories) <- c("id", "ext_state", "ext_grade", "dx", "event", "begin", "end", "year", "psa", "utility")
  enum(lifeHistories$ext_state) <- ext_stateT
  lifeHistories$state <- ext_state2state(lifeHistories$ext_state)
  lifeHistories <- lifeHistories[c(names(lifeHistories)[1], "state", names(lifeHistories)[-1])] # shift col order
  enum(lifeHistories$dx) <- diagnosisT
  enum(lifeHistories$event) <- eventT
  enum(diagnoses$ext_state) <- ext_stateT
  diagnoses$state <- ext_state2state(diagnoses$ext_state)
  enum(diagnoses$ext_grade) <- gradeT
  enum(diagnoses$dx) <- diagnosisT
  enum(diagnoses$tx) <- treatmentT
  enum <- list(stateT = stateT, ext_stateT = ext_stateT, eventT = eventT, screenT = screenT,
              diagnosisT = diagnosisT, psaT = psaT)
  out <- list(n=n,screen=screen,enum=enum,lifeHistories=lifeHistories,
              parameters=parameters, summary=summary,
              healthsector.costs=healthsector.costs, societal.costs=societal.costs,
              psarecord=psarecord, diagnoses=diagnoses,
              cohort=data.frame(table(cohort)),simulation.parameters=parameter,
              falsePositives=falsePositives,
              natural.history.summary=natural.history.summary)
  class(out) <- "fhcrc"
  out
}

## R --slave -e "options(width=200); require(microsimulation); callFhcrc(100,nLifeHistories=1e5,screen=\"screen50\")[[\"parameters\"]]"

#' @title Summarise simulation results
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname summary.fhcrc
#' @export
summary.fhcrc <- function(object, ...) {
    newobj <- object[c("n","screen")]
    with(object,
         structure(.Data=c(newobj,
                       with(object, list(
                           discountRate.costs = simulation.parameters$discountRate.costs,
                           discountRate.effectiveness = simulation.parameters$discountRate.effectiveness,
                           screening.tests = with(summary$events,
                                                sum(n[event == "toScreen"])) / n,
                           biopsies = with(summary$events,
                                           sum(n[event %in% c("toScreenInitiatedBiopsy",
                                                              "toClinicalDiagnosticBiopsy")])) / n,
                           negative.biopsies = with(summary$events,
                                                   sum(n[event %in% c("toScreenInitiatedBiopsy",
                                                                      "toClinicalDiagnosticBiopsy")])
                                                   - sum(n[event %in% c("toScreenDiagnosis",
                                                                        "toClinicalDiagnosis")])) / n,
                           screen.diagnosis = with(summary$events,
                                                   sum(n[event == "toScreenDiagnosis"])) / n,
                           over.diagnosis = with(summary$events,
                                                 sum(n[event == "toOverDiagnosis"])) / n,
                           cancer.deaths = with(summary$events,
                                                sum(n[event == "toCancerDeath"])) / n,
                           LE = sum(summary$pt$pt) / n,
                           QALE = sum(summary$ut$ut) / n,
                           healthsector.costs = sum(healthsector.costs$costs) / n,
                           societal.costs = sum(societal.costs$costs) / n))),
                   class="summary.fhcrc"))
}

#' @title Minus operator method for \code{summary.fhcrc} objects
#' @description Calculates the difference between two \code{summary.fhcrc}
#'     objects.
#' @param object1 The positive reference \code{summary.fhcrc} object.
#' @param object2 The negative \code{summary.fhcrc} object which is subtracted
#'     from \code{object1}.
#' @return A \code{summary.fhcrc} object representing the difference between the
#'     two compared \code{summary.fhcrc}.
#' @details Calculates the difference between two \code{summary.fhcrc} objects,
#'     e.g. the difference in life expectancy, quality of life etc. The
#'     exceptions are \code{discountRate.costs} and
#'     \code{discountRate.effectiveness} which are required two be the same for
#'     the two scenarios and the original value is returned. If \code{n} has the
#'     same value the original value is returned if they differ a list with both
#'     is returned.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scenario1 <- summary(callFhcrc(screen = "screenUptake"))
#'  scenario2 <-summary(callFhcrc(screen = "noScreening"))
#'  scenario1 - scenario2
#'  }
#' }
#' @rdname summary.fhcrc
#' @export
'-.summary.fhcrc' <- function(object1, object2) {
    if(object1$discountRate.costs != object2$discountRate.costs)
        stop("The two scenarios have different discount rates for the costs.")
    if(object1$discountRate.effectiveness != object2$discountRate.effectiveness)
        stop("The two scenarios have different discount rates for the effectiveness.")
    if(object1$n != object2$n)
        warning("The two scenarios have different number of men.")
    structure(.Data= c(list(n = if(object1$n == object2$n) {object1$n} else {list(object1$n, object2$n)},
                            screen = sprintf("%s-%s", object1$screen, object2$screen),
                            discountRate.costs = object1$discountRate.costs,
                            discountRate.effectiveness = object1$discountRate.effectiveness),
                       ## exempt special elements calc difference on all others
                       sapply(names(object1)[!names(object1) %in% c("n", "screen",
                                                                    "discountRate.costs",
                                                                    "discountRate.effectiveness")],
                              function(x) object1[[x]] - object2[[x]],
                              simplify = FALSE, USE.NAMES = TRUE)),
              class="summary.fhcrc")
}

#' @title Multiplaction with scalar operator method for \code{summary.fhcrc} objects
#' @description Calculates the results per a number of persons of a
#'     \code{summary.fhcrc} object.
#' @param obj The \code{summary.fhcrc} object.
#' @param perpersons a scalar to multiply the summary results with.
#' @return A \code{summary.fhcrc} object scaled for a number of persons.
#' @details N.b. this is unfortunatly not a commutative
#'     operator. Therefor the first parameter must be the
#'     \code{summary.fhcrc} object and the second parameter a nummeric
#'     scalar.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ## Example 1
#'  scenario1 <- summary(callFhcrc(screen = "screenUptake"))
#'  scenario1 * 1000 # the summery per thousand persons
#'
#'  ## Example 2
#'  scenario2 <-summary(callFhcrc(screen = "noScreening"))
#'  (scenario1 - scenario2) * 1000 # the difference per thousand persons
#'  }
#' }
#' @rdname summary.fhcrc
#' @export
'*.summary.fhcrc' <- function(obj, perpersons) {
    if(!is.numeric(perpersons))
        stop("The 'perpersons' variable needs to be a numeric.")
    structure(.Data= c(list(n = ifelse(is.list(obj$n),
                                       list(lapply(obj$n, function(x) x * perpersons)),
                                       obj$n * perpersons),
                            screen = sprintf("(%s)*%d", obj$screen, perpersons),
                            discountRate.costs = obj$discountRate.costs,
                            discountRate.effectiveness = obj$discountRate.effectiveness),
                       ## exempt special elements scale all others
                       lapply(obj[!names(obj) %in% c("n", "screen",
                                                     "discountRate.costs",
                                                     "discountRate.effectiveness")],
                              function(x) x * perpersons)),
              class="summary.fhcrc")
}

#' @title Divide with scalar operator method for \code{summary.fhcrc} objects
#' @description Calculates the results divided by a denominator \code{summary.fhcrc} object.
#' @param obj The \code{summary.fhcrc} object.
#' @param denominator a scalar to divide the summary results with.
#' @return A \code{summary.fhcrc} object scaled by a denominator.
#' @details N.b. this is unfortunatly not a commutative
#'     operator. Therefor the first parameter must be the
#'     \code{summary.fhcrc} object and the second parameter a nummeric
#'     scalar.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  }
#' }
#' @rdname summary.fhcrc
#' @export
'/.summary.fhcrc' <- function(obj, denominator) {
    if(!is.numeric(denominator))
        stop("The 'denominator' variable needs to be a numeric.")
    structure(.Data= c(list(n = ifelse(is.list(obj$n),
                                       list(lapply(obj$n, function(x) x / denominator)),
                                       obj$n * denominator),
                            screen = sprintf("(%s)/%d", obj$screen, denominator),
                            discountRate.costs = obj$discountRate.costs,
                            discountRate.effectiveness = obj$discountRate.effectiveness),
                       ## exempt special elements scale all others
                       lapply(obj[!names(obj) %in% c("n", "screen",
                                                     "discountRate.costs",
                                                     "discountRate.effectiveness")],
                              function(x) x / denominator)),
              class="summary.fhcrc")
}

#' @title Print a selection of the summarised simulation
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname print.summary.fhcrc
#' @export
print.summary.fhcrc <- function(x, ...) {
    obj <- x
    cat(sprintf(
"Screening scenario:            %s
Life expectancy:                %f
Discounted QALE:                %f
Discounted health sector costs: %f
Discounted societal costs:      %f
Discounted rate (effect.):      %f
Discounted rate (costs):        %f
", obj$screen, obj$LE, obj$QALE,
obj$healthsector.costs, obj$societal.costs,
obj$discountRate.effectiveness,
obj$discountRate.costs))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object1 PARAM_DESCRIPTION
#' @param object2 PARAM_DESCRIPTION
#' @param perspective PARAM_DESCRIPTION, Default: c("societal.costs", "healthsector.costs")
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ICER.fhcrc
#' @export
ICER.fhcrc <- function(object1, object2,
                       perspective = c("societal.costs",
                                       "healthsector.costs"), ...) {
    perspective <- match.arg(perspective)
    p1 <- object1$simulation.parameters
    p2 <- object2$simulation.parameters
    stopifnot(p1$discountRate.costs == p2$discountRate.costs)
    stopifnot(p1$discountRate.effectiveness == p2$discountRate.effectiveness)
    summary1 <- summary(object1, ...)
    summary2 <- summary(object2, ...)
    out <- list(ICER.QALE=(summary1[[perspective]] - summary2[[perspective]]) /
                    (summary1$QALE - summary2$QALE),
                delta.QALE=summary1$QALE - summary2$QALE,
                delta.costs=summary1[[perspective]] - summary2[[perspective]])
    if (p1$discountRate.costs == 0 && p2$discountRate.costs == 0 &&
        p1$discountRate.effectiveness == 0 && p2$discountRate.effectiveness == 0) {
        out <- c(out,
                 list(ICER.LE = (summary1[[perspective]] - summary2[[perspective]]) /
                          (summary1$LE - summary2$LE),
                      delta.LE = summary1$LE - summary2$LE))
    }
    out
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname print.fhcrc
#' @export
print.fhcrc <- function(x, ...)
    cat(sprintf("FHCRC prostate cancer model with %i individual(s) under scenario '%s'.\n",
                x$n, x$screen),
        ...)

## TODO: better solve issue below with testing.rate for noScreening scenario
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param scenarios PARAM_DESCRIPTION, Default: NULL
#' @param type PARAM_DESCRIPTION, Default: 'incidence.rate'
#' @param group PARAM_DESCRIPTION, Default: 'age'
#' @param age.breaks PARAM_DESCRIPTION, Default: NULL
#' @param year.breaks PARAM_DESCRIPTION, Default: NULL
#' @param cohort.breaks PARAM_DESCRIPTION, Default: NULL
#' @param age.weights PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @importFrom stats predict
#' @rdname predict.fhcrc
#' @export
predict.fhcrc <- function(object, scenarios=NULL, type = "incidence.rate",
                         group = "age", age.breaks = NULL, year.breaks = NULL,
                         cohort.breaks = NULL, age.weights = NULL, ...) {
    if(!inherits(object,"fhcrc")) stop("Expecting object to be an fhcrc object")
    if(!(is.null(scenarios) || all(sapply(scenarios,inherits,"fhcrc")) || inherits(object,"fhcrc")))
        stop("Expecting scenarios is NULL, a fhcrc object or a list of fhcrc objects")
    if(is.null(scenarios) && (grepl(".?rate.?ratio$|.?rr$", type, ignore.case = TRUE)))
        stop("More than one simulation object is needed to calculate a rate-ratio")

    ## Stripping of potential rate ratio option before matching
    abbr_type <- match.arg(sub(".?rr$|.?rate.?ratio$",
                               "", type, ignore.case = TRUE),
                           c("incidence.rate", "symptomatic.incidence.rate",
                             "screen.incidence.rate", "overdiagnosis.rate",
                             "testing.rate", "biopsy.rate", "metastasis.rate",
                             "pc.mortality.rate", "allcause.mortality.rate",
                             "prevalence", "ly", "life.years"))

    ## Allowing for several groups
    group <- match.arg(group,
                       c("state", "ext_state", "grade", "dx", "psa", "age", "year", "cohort"),
                       several.ok = TRUE)

    event_types <- switch(abbr_type,
                         incidence.rate = c("toClinicalDiagnosis", "toScreenDiagnosis"),
                         symptomatic.incidence.rate = "toClinicalDiagnosis",
                         screen.incidence.rate = "toScreenDiagnosis",
                         overdiagnosis.rate = "toOverDiagnosis",
                         testing.rate = "toScreen",
                         biopsy.rate = c("toClinicalDiagnosticBiopsy", "toScreenInitiatedBiopsy"),
                         metastasis.rate = "toMetastatic",
                         pc.mortality.rate = "toCancerDeath",
                         allcause.mortality.rate = c("toCancerDeath", "toOtherDeath"))

    interval2break <- function(interval_vector) {
        as.numeric(unique(unlist(regmatches(interval_vector,
                                            gregexpr("[0-9]+|Inf", interval_vector)))))
    }

    if(FALSE) {age <- year <- cohort <- Freq <- event <- Freq.y <- Freq.x <-
    scenario.x <- rate.x <- rate.y <- NULL} # To pass false-positive check note

    ## Manipulate input, make sure age.breaks exists from age standardisation
    if(!is.null(age.weights)) age.breaks <- interval2break(age.weights[,1])
    ## Manipulate input, make sure time group exists for specified time interval
    if(!is.null(age.breaks) && !("age" %in% group)) group <- c(group, "age")
    if(!is.null(year.breaks) && !("year" %in% group)) group <- c(group, "year")
    if(!is.null(cohort.breaks) && !("cohort" %in% group)) group <- c(group, "cohort")

    ## fast operations by group using base-R
    ## http://stackoverflow.com/questions/3685492/r-speeding-up-group-by-operations
    grp_apply = function(XS, INDEX, FUN, ..., simplify=T) {
        FUN = match.fun(FUN)
        if (!is.list(XS))
            XS = list(XS)
        as.data.frame(as.table(tapply(1:length(XS[[1L]]), INDEX, function(s, ...)
            do.call(FUN, c(lapply(XS, `[`, s), list(...))), ..., simplify=simplify)))
    }


    ## Fixes colnames after group operation
    name_grp <- function(x, grp = group) {names(x)[grep("^Var[0-9]+$", names(x))] <- grp; x}

    ## After a group operation the group will become a factor, if the time scale
    ## is not reported as an interval we convert it back to numeric. N.b. This
    ## needs to be called after name_grp().
    numeric_time_scale <- function(df, group) {
        within(df, {
            if("age" %in% group && is.null(age.breaks)) age <- as.numeric(levels(age))[age] #important factor conversion
            if("year" %in% group && is.null(year.breaks)) year <- as.numeric(levels(year))[year] #important factor conversion
            if("cohort" %in% group && is.null(cohort.breaks)) cohort <- as.numeric(levels(cohort))[cohort] #important factor conversion
        })
    }
    categorise_time <- function(df, age.breaks, year.breaks, cohort.breaks) {
        cut_t <- function(time, time.breaks) {
            ## using first column or single vector as breaks
            as.factor(cut(time, breaks = data.frame(time.breaks)[,1],
                          right = FALSE, dig.lab=4, ordered_result = TRUE))
        }
        if(!is.null(age.breaks)) df <- transform(df, age = cut_t(age, age.breaks))
        if(!is.null(year.breaks)) df <- transform(df, year = cut_t(year, year.breaks))
        if(!is.null(cohort.breaks)) df <- transform(df, cohort = cut_t(cohort, cohort.breaks))
        df
    }

    ## For standardising over a time scale e.g. ages
    standardise_time <- function(df, timestr = "age", parstr = "rate",
                                time.weights = age.weights) {
        if(is.null(time.weights)) {
            return(df)
        } else {
            collapsed_grps <- c(group[!group == timestr], "scenario")
            ## Scale with weight
            weight <- function(time, par, time.weights) {
                par * sapply(t(time), {function(x) time.weights[x == time.weights[,1],2]})
            }
            ## Sum over the groups (at least scenario)
            collapse <- function(df, collapsed_grps, parstr) {
                within(with(df, grp_apply(eval(parse(text = parstr)),
                                          lapply(as.list(collapsed_grps),
                                          {function(x) eval(parse(text = x))}),
                                          sum, na.rm = TRUE)) ,{
                                              assign(eval(substitute(parstr)),
                                                     Freq / sum(time.weights[,2]))
                                              rm(Freq)
                                          })
            }
            ## Standardise and rename
            ## browser()
            return(numeric_time_scale(name_grp(collapse(
                within(df,
                {assign(eval(substitute(parstr)), weight(eval(parse(text =timestr)),
                                                         eval(parse(text = parstr)),
                                                         time.weights))}),
                collapsed_grps, parstr), collapsed_grps), collapsed_grps))
        }
    }

    ## Calculates rates of specific events by specified groups
    calc_rate <- function(object, event_types, group){
        pt <- with(categorise_time(object$summary$pt, age.breaks, year.breaks,
                                  cohort.breaks),
                  numeric_time_scale(name_grp(grp_apply(pt,
                                      lapply(as.list(group), function(x) eval(parse(text = x))),
                                      sum)), group))
        ## TODO: if subset has no dim replace with zeros, temp fix below:
        if(!any(object$summary$event$event %in% event_types)) {
            stop(paste("The event(s)", paste(event_types, collapse = ", "),
                       "was not found in the", object$screen, "scenario"))
        }
        events <- with(subset(
            categorise_time(object$summary$events, age.breaks, year.breaks,
                            cohort.breaks), event %in% event_types),
            numeric_time_scale(name_grp(grp_apply(n,
                                                  lapply(as.list(group),
                                                  {function(x) eval(parse(text = x))}),
                                                  sum)), group))
        within(merge(pt, events, by = group, all = TRUE),{
            rate <- ifelse(is.na(Freq.y) & !is.na(Freq.x), 0, Freq.y/Freq.x) #no events but some pt -> 0
            n <- Freq.y
            pt <- Freq.x
            rm(Freq.x,Freq.y)})
    }

    ## Calculate prevalences by specified groups
    calc_prev <- function(object, group){
        within(with(categorise_time(object$summary$prev, age.breaks,
                                    year.breaks, cohort.breaks),
                    numeric_time_scale(name_grp(grp_apply(count,
                                       lapply(as.list(group),
                                              function(x) eval(parse(text = x))), sum)), group)), {
                                                  prevalence <- Freq/object$n
                                                  rm("Freq")})
    }

    ## Calculate prevalences by specified groups
    calc_ly <- function(object, group){
        within(with(categorise_time(object$summary$pt, age.breaks,
                                    year.breaks, cohort.breaks),
                    numeric_time_scale(name_grp(grp_apply(pt,
                                       lapply(as.list(group),
                                              function(x) eval(parse(text = x))), sum)), group)), {
                                                  ly <- Freq/object$n
                                                  rm("Freq")})
    }

    ## Calculates the outcome in the passed function for all
    ## simulation objects in the 'scenarios' list. Then the object
    ## outcome (e.g. rates or prev) are for the scenarios are added as
    ## rows and the scenario name as a column.
    predict_scenarios <- function(scenarios, calc_outcome, ...) {
        do.call(rbind, lapply(scenarios,
        {function(object, ...)
            cbind(calc_outcome(object, ...), scenario = object$screen)}, ...))
    }

    ## Input checks allow for scenarios to be a single fhcrc object or
    ## list of fhcrc objects. Now make sure scenarios is a list.
    if(inherits(scenarios, "fhcrc")) {scenarios <- list(scenarios)}

    ## Rate-ratio if type ends with rate.ratio or RR
    if(grepl(".?rate.?ratio$|.?rr$", type, ignore.case = TRUE)){
        scenario_rates <- standardise_time(predict_scenarios(unique(scenarios),
                                                            calc_rate,
                                                            event_types, group),
                                          timestr = "age",
                                          parstr = "rate",
                                          time.weights = age.weights)
        reference_rate <- standardise_time(predict_scenarios(list(object),
                                                            calc_rate,
                                                            event_types, group),
                                          timestr = "age",
                                          parstr = "rate",
                                          time.weights = age.weights)

        ## Alt. check if it is part of the rates
        if(!is.null(age.weights)) group <- group[!group == "age"]
        within(merge(scenario_rates, reference_rate, by = group),{
            scenario <- scenario.x
            rate.ratio <- rate.x/rate.y
            rate.ratio[!is.finite(rate.ratio)] <- NaN
            rm(list=ls(pattern=".x$|.y$"))})


        ## Prevalence if type ends with rate.ratio or RR
    } else if(grepl(".?prev$|.?prevalence$", type, ignore.case = TRUE)){
        standardise_time(predict_scenarios(unique(c(list(object),scenarios)), calc_prev, group),
                         timestr = "age", parstr = "prevalence", time.weights = age.weights)

        ## Life-years if type is ly or life*years
    } else if(grepl("^ly$|^life.?years$", type, ignore.case = TRUE)){
        standardise_time(predict_scenarios(unique(c(list(object),scenarios)), calc_ly, group),
                         timestr = "age", parstr = "LY", time.weights = age.weights)

        ## Defaults to plain rates. If reference object exist add it to
        ## scenario list and remove duplicates.
    } else {
        standardise_time(predict_scenarios(unique(c(list(object),scenarios)), calc_rate, event_types, group),
                         timestr = "age", parstr = "rate", time.weights = age.weights)
    }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: c("incidence.rate", "testing.rate", "biopsy.rate", "metastasis.rate",
#'    "pc.mortality.rate", "allcause.mortality.rate")
#' @param plot.type PARAM_DESCRIPTION, Default: 'l'
#' @param add PARAM_DESCRIPTION, Default: FALSE
#' @param xlab PARAM_DESCRIPTION, Default: 'Age (years)'
#' @param ylab PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @importFrom graphics plot
#' @rdname plot.fhcrc
#' @export
plot.fhcrc <- function(x, type=c("incidence.rate", "testing.rate",
                                 "biopsy.rate", "metastasis.rate",
                                 "pc.mortality.rate",
                                 "allcause.mortality.rate"),
                       plot.type="l", add=FALSE, xlab="Age (years)",
                       ylab=NULL, ...) {
    type <- match.arg(type)
    if (is.null(ylab)) {ylab <- switch(type,
                                       incidence.rate="Prostate cancer incidence rates per 100,000",
                                       testing.rate="PSA rates per 1000",
                                       biopsy.rate="Biopsies per 1000",
                                       metastasis.rate="Metastatic onset per 100,000",
                                       pc.mortality.rate="Cancer mortality rates per 100,000",
                                       allcause.mortality.rate="All cause mortality rates per 100,000")}
    rates <- predict(object = x, type = type)
    rates$rate = rates$rate*switch(type, testing.rate=1000, biopsy.rate=1000, incidence.rate=1e5, metastasis.rate=1e5,pc.mortality.rate=1e5,allcause.mortality.rate=1e5)
    if (!add) plot(rate~age, data=rates, type=plot.type, xlab=xlab, ylab=ylab, ...) else lines(rate~age, data=rates,  ...)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname lines.fhcrc
#' @importFrom graphics lines
#' @export
lines.fhcrc <- function(x,...) {
    plot(x, ..., add=TRUE)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object1 PARAM_DESCRIPTION
#' @param object2 PARAM_DESCRIPTION
#' @param startAge PARAM_DESCRIPTION, Default: 50
#' @param stopAge PARAM_DESCRIPTION, Default: Inf
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname nn.fhcrc
#' @export
nn <- function(object1, object2, startAge = 50, stopAge = Inf, ...){
  UseMethod("nn")
}

#' @export
nn.fhcrc <- function(object1, object2, startAge = 50, stopAge = Inf, ...) {
    pNNS <- function(thisScenario) {
        with(thisScenario$summary$events,
             sum(n[event=="toCancerDeath" & age>=startAge
                   & age<stopAge])) / # divided by
            with(thisScenario$summary$prev,
                        sum(count[age==round(startAge)]))
    }
    pNND <- function(thisScenario) {
        with(thisScenario$summary$events,
             sum(n[event=="toCancerDeath" & age>=startAge & age<stopAge])) / # divided by
            with(thisScenario$summary$events,
                 sum(n[event %in% c("toScreenDiagnosis","toClinicalDiagnosis")
                       & age>=startAge & age<stopAge]))
    }
    NNS <- 1 / (pNNS(object2) - pNNS(object1)) #number needed to screen to prevent 1 PCa death
    NND <- 1 / (pNND(object2) - pNND(object1)) #number needed to detect to prevent 1 PCa death
    ## Include additional number needed to treat (NNT) [Gulati 2011] to show overdiagnosis?
    structure(.Data = list(NNS=NNS,NND=NND), class = "nn.fhcrc")
}
