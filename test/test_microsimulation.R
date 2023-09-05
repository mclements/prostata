## try(detach("package:microsimulation", unload=TRUE))
## require(microsimulation)
## microsimulation:::.testPackage()

library(prostata)
parms=modifyList(prostata:::XiaoyangParameters(2022),
                 list(AI_assisted_pathology=TRUE,
                      start_screening=50,
                      stop_screening=70,
                      screening_interval=4))
set.seed(12345)
sim = callFhcrc(1e4, screen="regular_screen",
                nLifeHistories = 1e8,
                mc.cores=6,
                parms=parms)

AIcost.Fun <- function(object,ratio=1) {
    index = object$healthsector.costs$item == "AI pathology"
    object$healthsector.costs$costs[index] = object$healthsector.costs$costs[index]*ratio
    object
}
summary(AIcost.Fun(sim,ratio=1))
summary(AIcost.Fun(sim,ratio=10))
summary(AIcost.Fun(sim,ratio=20))


AIcost.Fun <- function(object,per=1e4,subset=TRUE,from=50,ratio=1) {
    n = sum(subset(object$summary$prev,age==from)$count) # number alive at 'from'
    TotalCosts <- sum(subset(object$healthsector.costs,age>=from)$costs)
    AIcosts <- sum(subset(object$healthsector.costs,age>=from & item == "AI pathology")$costs)
    (TotalCosts + AIcosts*(ratio-1))/n
}
AIcost.Fun(sim,ratio=1)
AIcost.Fun(sim,ratio=10)
AIcost.Fun(sim,ratio=20)

table(sim$healthsector.costs$item)
table(sim$societal.costs$item)

## table(sim$lifeHistories$event)
## tmp = subset(sim$lifeHistories, event=="toScreenInitiatedBiopsy")
tmp=subset(sim$summary$events, event=="toScreenInitiatedBiopsy")
pReduced = unlist(sim$simulation.parameters[c("pReducedBxCostG1","pReducedBxCostG2","pReducedBxCostG4plus","pReducedBxCostG0")])
tmp$pReduced = pReduced[tmp$grade]
sum(tmp$n)
table(sim$summary$events$grade)
myreport3 <- function(object,names="toScreenInitiatedBiopsy",per=1e4,subset=TRUE,from=50) {
    n = sum(subset(object$summary$prev,age==from)$count) # number alive at 'from'
    pReduced = unlist(object$simulation.parameters[c("pReducedBxCostG1","pReducedBxCostG2","pReducedBxCostG4plus","pReducedBxCostG0")])
    data = subset(object$summary$events, event %in% names & subset & age>=from)
    list(men_biopsied=sum(data$n)/n*per, reduced_cores=12*sum(data$n*pReduced[data$grade])/n*per,
         p_reduced_cores = sum(data$n*pReduced[data$grade])/sum(data$n))
}
myreport3(sim)

## Trust: probabilistic analysis
4/19
set.seed(12345)
sims1 = rbinom(1e5,19,4/19)/19
##
library(boot)
y=c(rep(1,4),rep(0,19-4))
boot1 = boot(y, function(y,i) sum(y[i])/length(y), R=100000)
##
par(mfrow=1:2)
plot(table(boot1$t[,1]))
plot(table(sims1))

## Base strategies
library(prostata)
if (FALSE)
    save.image("~/Downloads/trust_20230512.RData")
load("~/Downloads/trust_20230512.RData")
n = 1e6
## Probase
parm = modifyList(prostata:::TrustParameters(MRI_diagnostics=FALSE), # or TRUE
                  list(start_screening=45, # or 50
                       stop_screening=60,
                       psaThreshold=3, # for referral
                       risk_psa_threshold=1.5, # rescreening
                       risk_lower_interval=5,  # rescreening
                       risk_upper_interval=2)) # rescreening
set.seed(12345)
strategy1 = callFhcrc(n, screen="probase",
                      pop=1950,
                      parm=parm,
                      mc.cores=5)
set.seed(12345)
strategy2 = callFhcrc(n, screen="probase",
                      pop=1950,
                      parm=modifyList(parm,list(start_screening=50)),
                      mc.cores=5)
parm = modifyList(prostata:::TrustParameters(MRI_diagnostics=TRUE),
                  list(start_screening=45, # or 50
                       stop_screening=60,
                       psaThreshold=3, # for referral
                       risk_psa_threshold=1.5, # rescreening
                       risk_lower_interval=5,  # rescreening
                       risk_upper_interval=2)) # rescreening
set.seed(12345)
strategy3 = callFhcrc(n, screen="probase",
                      pop=1950,
                      parm=parm,
                      mc.cores=5)
set.seed(12345)
strategy4 = callFhcrc(n, screen="probase",
                      pop=1950,
                      parm=modifyList(parm,list(start_screening=50)),
                      mc.cores=5)
## Germany 2021 guidelines
parm = modifyList(prostata:::TrustParameters(MRI_diagnostics=FALSE, DRE=TRUE), # or FALSE
                  list(start_screening=45,
                       stop_screening=75,
                       risk_psa_threshold_lower=1.5, # rescreening
                       risk_psa_threshold_moderate=2, # rescreening
                       psaThreshold=4, # for referral
                       risk_lower_interval=5, # rescreening
                       risk_upper_interval=2)) # rescreening
set.seed(12345)
strategy5 = callFhcrc(n, screen="germany_2021",
                      pop=1950,
                      parm=parm,
                      mc.cores=5)
parm = modifyList(prostata:::TrustParameters(MRI_diagnostics=TRUE, DRE=TRUE),
                  list(start_screening=45,
                       stop_screening=75,
                       risk_psa_threshold_lower=1.5, # rescreening
                       risk_psa_threshold_moderate=2, # rescreening
                       psaThreshold=4, # for referral
                       risk_lower_interval=5, # rescreening
                       risk_upper_interval=2)) # rescreening
set.seed(12345)
strategy6 = callFhcrc(n, screen="germany_2021",
                      pop=1950,
                      parm=parm,
                      mc.cores=5)
set.seed(12345)
strategy7 = callFhcrc(n, screen="noScreening",
                      pop=1950,
                      parm=parm,
                      mc.cores=5)

## Quick plot of the cost-efficiency frontier
strategies = list(strategy1,
                  strategy2,
                  strategy3,
                  strategy4,
                  strategy5,
                  strategy6,
                  strategy7)

d = t(sapply(lapply(strategies, summary),
             function(object) c(object$QALE, object$healthsector.costs)))
d = data.frame(QALE=d[,1], Costs=d[,2],
               labels=c("ProBase SBx, 45 years",
                        "ProBase SBx, 50 years",
                        "ProBase MRI, 45 years",
                        "ProBase MRI, 50 years",
                        "2021 Guidelines, SBx",
                        "2021 Guidelines, MRI",
                        "No screening"),
               pos=c(4,4,3,2,2,2,1))
with(d, {
    plot(QALE, Costs,
         xlim=c(min(QALE)-0.0002,max(QALE)+0.0009),
         ## ylim=c(240,325),
         ylab="Life-time discounted costs (€)",
         xlab="Life-time discounted utilities (QALYs)")
    text(QALE, Costs, labels=labels,pos=pos)
    })
with(d[c(7,2,5), ], lines(QALE, Costs))

ICER(strategy2, strategy7, perspective="healthsector.costs")
ICER(strategy2, strategy5, perspective="healthsector.costs")
ICER(strategy2, strategy7, perspective="healthsector.costs", from=35)
ICER(strategy6, strategy2, perspective="healthsector.costs", from=35)


## Trust: model extensions to include DRE
library(prostata)
mortality_1990_2021 =
    data.frame(age=rep(seq(40,85,by=5),2021-1990+1),
               year=rep(1990:2021,each=10),
               rate = read.table(text="0.167662722959985	0.233291198350787	0.258712991568544	0.285536029650061	0.175569996773023	0.208606477439731	0.172109811567294	0.203127486195795	0.500093850946028	0.294173597723489	0.12733127650878	0.216672449757528	0.090273402025374	0.293283514240381	0.25679757446131	0.139944772195101	0.332710978963239	0.24930879137591	0.167655882947133	0.0855974085669881	0.0294501278282798	0.184113575982752	0.289002584646449	0.101694673944841	0.144674303782039	0.114892169868839	0.0794026380732474	0	0.124418034642958	0.0410011160503789	0	0
0.832281999335561	0.962394622600798	1.25197028824112	1.05635528693145	1.06710375236363	1.0942037920895	1.0159865483381	1.1928339012341	1.41579915690975	1.32340517158127	0.886825148232823	1.05610970022678	0.901032652730233	0.612865685037429	1.26621251844339	0.875325450865551	0.788241455383799	0.79769160322206	0.598646818730941	0.961356106402894	0.973003168556201	0.764811206353712	0.620679051094582	0.594053524222532	0.803296960238235	0.440632678024415	0.726712277955334	0.379417486656361	0.500128866537945	0.462268736729776	0.641766140418431	0.514591030757502
3.31231297141909	4.46330676695523	4.20459160515062	4.33401928149781	3.63728704189548	3.90605303283712	3.52181228503095	3.99954241794272	3.67352989925344	3.48566267521522	4.21111308401518	3.83108079225138	4.07277511211658	3.51088448891233	2.75678739435072	3.937549014282	4.16166211718589	2.81640229918249	2.87305477549004	3.73653262798567	3.3360856039566	3.17784175256027	3.17877359767301	2.99569689538812	2.60596119546101	3.15502702902972	2.70591515902094	2.66571986346977	2.16032511756006	2.67994381915624	2.86662424850854	2.17449162070373
9.81706499074708	11.7080752267421	13.154357809971	11.6807293107016	11.7486783643374	11.9628248263391	12.5561154004627	11.0176089266526	12.1268936348797	11.5304504386331	10.634658575402	10.949849688427	12.524389897773	11.2889462934093	11.2701453848755	11.1148512033182	10.8238638508958	11.5950790484518	10.1472044451754	9.6419175478308	9.87728611686653	9.79887064247049	10.1580277707133	9.14875881838697	8.84651995597972	8.9171402182573	8.79670677854307	7.69788971532946	8.22827937605239	7.39915467690257	6.94281670667193	8.25980544350469
29.4632462655976	30.4893341424481	32.3284025281776	30.6172996408244	32.6537016108669	30.1474213682629	31.5843736196245	27.3825687131334	29.3115481271955	28.4975714859442	30.1320453185962	28.870487997535	28.519706905166	27.7960625389635	26.9143537965005	25.4906555373304	24.8616052030085	26.7176268223868	26.500641364146	25.3265189677537	26.6771200033778	25.8589994659071	23.9796469568938	26.0220114151534	24.26119415939	24.9787281151371	22.9711638749388	21.916045226324	21.9939954109414	20.1059084290868	21.047306790221	19.6751438396043
68.9927674414609	69.5976844440987	74.1020437060462	73.2864169977584	75.2539821898909	75.1600593707095	74.5602763583261	69.3304734513035	67.9478294834072	63.5018196878773	63.0605043877595	58.2326468768182	61.123266551947	60.0795745027968	59.7327121473333	53.5469665623051	56.3004136565683	54.6078886741615	57.514531343833	54.1482258456866	53.3088680361238	52.7255185721177	54.2849230823995	52.8808413082477	51.8115461181155	50.0527267053623	50.424886147433	49.7574325164822	50.4406676106531	48.7861499698469	48.1291677891761	48.5048342008786
159.280451532216	160.395383433096	154.887498160389	154.857837011207	157.318404640256	146.578499231165	145.310133809478	133.758468357085	132.468368534429	134.007435716842	124.29077256349	126.124321244305	126.880100535433	121.742762240448	109.132560327127	107.127522141656	114.05151405756	106.359208115476	102.725661194949	102.880207610604	103.91412973073	99.8880043587493	97.2314776871651	92.870512226908	96.7984054343742	98.7445748740206	97.2192375002098	92.5624840789006	97.7742844727075	93.04202785189	97.1902134692832	90.5358492845026
285.681850495369	295.200880456509	296.967193309275	295.259234640585	300.49056911453	320.708453519073	298.750085802897	281.03254647846	274.95098383739	257.418814710988	242.247222174071	234.567364715684	222.508176505281	224.618856954913	216.866809233975	210.52451339215	204.081298657515	192.599575018773	201.016496669277	199.990092380561	187.719241834749	194.565222118963	170.462488073733	179.038917573188	171.622372306523	164.172883872802	168.971391285987	171.856547952017	172.836951128588	168.676324062533	162.086057367732	159.313488456058
492.278382781842	511.330244405365	531.422444421839	531.553114289592	534.393510681469	515.624010210147	509.718431880322	489.796710117939	470.020039312774	462.945240777575	483.195119282599	443.558058073315	439.132806427169	424.991272435846	404.858807792068	380.985003071229	365.530942723605	353.878820733979	371.003580792062	343.552804300173	352.328247686054	354.823813585864	324.033969247384	330.585205515462	321.761284250258	309.393752577648	310.747569587151	289.848218483895	299.072445644068	297.251857629879	293.076648008816	287.99470508644
765.081359599392	792.235618390735	809.614876518796	851.349400452916	892.17007055141	926.436550190985	920.441264318856	909.250709507018	872.911809989451	822.274506498673	789.711854853391	793.194178499023	796.556517459719	813.45638197511	744.125924941651	762.014639925632	718.107091232786	650.25563174523	668.225927868333	655.213054808512	663.79224162755	696.339180905129	660.114560044299	648.815176275806	653.950644029284	629.666196816807	636.562416490574	608.971936564805	618.425655814135	618.813241295458	612.639109203823	592.447564012749", sep="\t") |> unlist() |> "names<-"(NULL))
incidence_1999_2019 =
    data.frame(age=rep(seq(40,85,by=5),21),
           year=rep(1999:2019,each=10),
           rate=read.table(text="1.8	2	1.9	2.3	2.6	2.2	2.2	2.6	2.9	2.4	2.5	3	3	3	2.3	2.5	2.1	2.3	1.7	2.3	2
13.3	13.4	13.8	16	18.4	17.1	16.9	18.6	19.1	17.9	17.8	20	19.9	18.1	17.9	16.3	14.9	15.4	15.5	16	16.3
47.2	48.1	49.2	57.3	65	61.7	61.3	66.8	70.3	68.3	64.9	67.1	68	63.7	60.2	56.7	54.2	52.3	54.6	55.6	57.3
135.6	137.1	140.2	162.9	190.9	177	171.3	187	200	187.8	184.1	187.9	184.1	174.6	158.5	153.8	144.7	144.6	150.9	152.7	157.6
261.3	264.1	276.1	341.6	421.8	401.4	395	432.1	448.9	403.6	400.7	389.6	387.7	364.5	323	307.6	296.5	299.7	310.7	311.9	334.9
416.9	441	472.2	527.6	652.6	619.2	580	639	695.1	657.1	655.7	657.3	650.2	615.5	534.8	517.7	493	506.6	518.3	535.7	566.7
587.1	598.4	600.2	725.7	834.7	738.4	758.2	800.5	840.5	789.2	763.9	788.7	780.1	731.8	673.3	640.5	653.5	656.9	698.2	701.8	720.4
743.1	717.8	738.2	788.2	934	793.9	803.5	857.7	881.7	821.6	775.1	767.7	782.1	737.4	669.3	644.1	663.1	679.2	716.7	733.1	784.7
842.4	863	756.7	894	911.5	846.6	791.8	806.1	853.9	813.4	747.7	725.8	749.5	671.2	623.1	591.2	590.5	584.8	580.3	603.9	619.8
1151	1090.5	1043.1	1117.7	1094.5	901	914.1	874.4	903.9	840.9	839.3	789.2	844.2	739.7	720.9	698.3	663.5	610.4	591.7	548.2	550.2", sep="\t") |> unlist()) |> "rownames<-"(NULL)

if (FALSE) {
    n = 1e5
    set.seed(12345)
    parm = prostata:::TrustParameters(MRI_diagnostics=FALSE)
    parm$prtx = prostata:::germany_observed_tables$prtx
    sim1 = callFhcrc(n, screen="germany_observed",
                     pop=1950,
                     parm=parm,
                     mc.cores=5)
    set.seed(12345)
    sim2 = callFhcrc(n, screen="germany_2021",
                     pop=1950,
                     parm=prostata:::TrustParameters(),
                     mc.cores=5)
    set.seed(12345)
    sim0 = callFhcrc(n, screen="noScreening",
                     pop=1950,
                     parm=parm,
                     mc.cores=5)
    set.seed(12345)
    sim3 = callFhcrc(n, screen="regular_screen",
                     pop=1950,
                     parm=modifyList(prostata:::TrustParameters(),
                                     list(start_screening=50, stop_screening=75, screening_interval=2)),
                     mc.cores=5)
}

## The fit for mortality is okay...
library(dplyr)
pred1 = predict(sim1, type="pc.mortality.rate") |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred1, plot(age+0.5, rate, lty=1, type="l", xlab="Age (years)", ylab="Rate per 100,000",
                 log="y"))
pred2 = predict(sim2, type="pc.mortality.rate") |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred2, lines(age+0.5, rate, lty=2))
pred0 = predict(sim0, type="pc.mortality.rate") |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred0, lines(age+0.5, rate, lty=3))
pred3 = predict(sim3, type="pc.mortality.rate") |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred3, lines(age+0.5, rate, lty=4))
years = unique(mortality_1990_2021$year)
years = seq(1990,2020,by=5)
for (i in 1:length(years))
    with(subset(mortality_1990_2021, year==years[i]),
         lines(age+2.5, rate, col=i+1))

## Incidence has been calibrated
library(dplyr)
pred1 = predict(sim1) |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred1, plot(age+0.5, rate, lty=1, type="l", xlab="Age (years)", ylab="Rate per 100,000",
                 log="y"))
pred2 = predict(sim2) |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred2, lines(age+0.5, rate, lty=2))
pred0 = predict(sim0) |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred0, lines(age+0.5, rate, lty=3))
pred3 = predict(sim3) |> mutate(age=pmin(85,age)) |>
    group_by(age) |>
    summarise(rate = sum(n)/sum(pt)*1e5)
with(pred3, lines(age+0.5, rate, lty=4))
years = unique(incidence_1999_2019$year)
for (i in 1:length(years))
    with(subset(incidence_1999_2019, year==years[i]),
         lines(age+2.5, rate, col=i+1))



library(prostata)
fitted = list(weibull_onset_shape=exp(-0.0754039582652198),
              weibull_onset_scale=exp(5.2319848747848),
              beta7=exp(-2.35969408469582),
              beta8=exp(-1.30383055358671))
germany_2021_tables =
    list(prtx =
             data.frame(DxY=2008,
                        Age=c(50,50,50,75,75,75),
                        G=as.integer(c(6,7,8,6,7,8)-6),
                        CM=c(1,0,0,1,1,1),
                        RP=c(0,0.5,0.5,0,0,0),
                        RT=c(0,0.5,0.5,0,0,0)))
## Germany cancer incidence 2014-2019
rates=read.table(text="Age	Cases	Rate	Pop
40	323	2.15	14992758
45	2944	15.72	18731382
50	11492	55.10	20856325
55	28090	150.84	18622558
60	48318	310.57	15557976
65	67263	524.08	12834517
70	73839	676.33	10917643
75	77804	703.38	11061437
80	39854	596.25	6684070
85	25127	606.40	4143644", sep="\t", header=TRUE)
## Saarland cancer incidence by Gleason, 2014-2019
gleason=read.table(text="Age	GS6	GS7	GS8plus
40	2	0	1
45	8	19	4
50	28	80	25
55	103	208	74
60	170	359	129
65	191	444	223
70	211	430	260
75	255	462	334
80	114	193	238
85	26	44	102", sep="\t", header=TRUE) |> subset(Age>=45)
## Germany prostate cancer mortality rates (for validation)
mortality_2014_2019_mean =
    data.frame(age=seq(40,85,by=5),
               rate=c(0.0840647104029104,0.552076167690344,2.66214869728296,
                      8.31428178684408,22.704505869303,50.2122348446488,96.0235023686838,
                      169.689411768075,304.679188028816,627.731681835177))

set.seed(12345)
sim1 = callFhcrc(1e5, screen="germany_observed",
                 pop=1950,
                 parm=modifyList(prostata:::TrustParameters(MRI_diagnostics=FALSE),
                                            prostata:::germany_2021_tables) |> modifyList(fitted),
                 mc.cores=6)




set.seed(12345)
sim1 = callFhcrc(1e5, screen="germany_2021",
                 pop=1950,
                 parm=modifyList(modifyList(prostata:::TrustParameters(MRI_diagnostics=TRUE),
                                            germany_2021_tables),
                                 list(stop_screening=75,
                                      weibull_onset=TRUE)) |> modifyList(fitted),
                 mc.cores=6)
set.seed(12345)
sim2 = callFhcrc(1e5, screen="germany_observed",
                 pop=1950,
                 parm=modifyList(modifyList(prostata:::TrustParameters(MRI_diagnostics=FALSE),
                                            prostata:::germany_observed_tables),
                                 list(stop_screening=75,
                                      weibull_onset=TRUE)) |> modifyList(fitted),
                 mc.cores=6)
plot(sim1)
plot(sim2, add=TRUE, col=2)
with(rates, lines(Age+2.5-0.5, Rate, lty=2))
legend("topleft", legend=c("Germany 2018 guidelines","Germany observed"),
       lty=1, col=1:2, bty="n")

## Optimisation code:)
library(minqa)
library(dplyr)
Base = modifyList(modifyList(prostata:::TrustParameters(MRI_diagnostics=FALSE),
                             germany_observed_tables),
                  list(stop_screening=75, weibull_onset=TRUE))
mc.cores=3
theta=log(c(1,80,0.06840404,0.1989443))
nsim=1e4
negll = function(theta=log(c(1,80,0.06840404,0.1989443))) {
    set.seed(12345)
    base = data.frame(age=seq(45,85,by=5))
    sim2 <- callFhcrc(n=nsim,screen="germany_observed",pop=1950,
                             parms=modifyList(Base, list(weibull_onset_shape = exp(theta[1]),
                                                         weibull_onset_scale= exp(theta[2]),
                                                         beta7 = exp(theta[3]),
                                                         beta8 = exp(theta[4]))),
                             mc.cores=mc.cores, print.timing=FALSE)
    pred = predict(sim2) |>
        subset(age>=40) |>
        transform(age5=pmin(85,floor(age/5)*5)) |>
        group_by(age5) |>
        summarise(mu = sum(rate*pt)/sum(pt)) |> data.frame() |> mutate(age=age5, age5=NULL)
    pred = merge(pred, base, all=TRUE, by="age") |> mutate(mu=ifelse(is.na(mu), 1e-5, mu))
    ons2 = subset(rates, Age>=35)
    library(dplyr)
    gs6 = sim2$summary$events |>
        filter(grade=="Gleason_le_6" & event %in% c("toClinicalDiagnosis","toScreenDiagnosis")) |>
        mutate(age = pmin(85,floor(age/5)*5)) |>
        group_by(age) |>
        summarise(Freq6 = sum(n))
    gs7 = sim2$summary$events |>
        filter(grade=="Gleason_7" & event %in% c("toClinicalDiagnosis","toScreenDiagnosis")) |>
        mutate(age = pmin(85,floor(age/5)*5)) |>
        group_by(age) |>
        summarise(Freq7 = sum(n))
    gs8 = sim2$summary$events |>
        filter(grade=="Gleason_ge_8" & event %in% c("toClinicalDiagnosis","toScreenDiagnosis")) |>
        mutate(age = pmin(85,floor(age/5)*5)) |>
        group_by(age) |>
        summarise(Freq8 = sum(n))
    mu = merge(merge(merge(base,gs6,all=TRUE,by="age"),
                     gs7,all=TRUE,by="age"),
               gs8, all=TRUE, by="age") |>
        filter(age>=45) |>
        mutate(Freq6=ifelse(is.na(Freq6),1,Freq6),
               Freq7=ifelse(is.na(Freq7),1,Freq7),
               Freq8=ifelse(is.na(Freq8),1,Freq8))
    ## log-likelihood           
    ll1 = sum(sapply(1:nrow(gleason), function(i)
        dmultinom(unlist(gleason[i,-1]), prob=unlist(mu[i,-1]), log=TRUE)))
    val = -sum(dpois(x=ons2$Cases, lambda=ons2$Pop*pred$mu, log=TRUE)) - ll1
    print(theta)
    print(val)
    val
}
## optim(log(c(1,80)), negll)
nsim=1e5
negll()
bobyqa(c(0.143657727284436, 5.29666758041499, -2.9240240511791, -2.54198318725383), negll) # n=1e3
bobyqa(c(-0.0851374969515032, 5.27047615268318, -3.21044004054539, -1.35280299939307), negll) # n=1e4
bobyqa(c(-0.0754039582652198, 5.2319848747848, -2.35969408469582, -1.30383055358671), negll) # n=1e5


## set.seed(12345)
## sim3 = callFhcrc(1e4, screen="regular_screen",
##                  parm=modifyList(c(germany_2021_tables,
##                                    prostata:::TrustParameters()),
##                                  list(MRI_screen=TRUE, screening_interval=2)),
##                  mc.cores=6)
## plot(sim3, add=TRUE, col=3)

## Edna: calibration to ONS data
## 2017
library(prostata)
library(dplyr)
library(minqa)
nsim=1e6
mc.cores=6
ons = data.frame(age = seq(35,90,by=5),
                 rate = c(0.4,4,21.3,73.8,188.7,339.1,579.4,713.9,810,677.5,666.4,663.2), # correct
                 pop = c(3642643,3442758,3850108,3907196,3479034,2982920,2890646,2604535,1813420,1369854,856812,495244), # correct
                 cases = c(15,138,820,2884,6565,10115,16748,18594,14689,9281,5710,3284)) # approximations
Base <-
    modifyList(prostata:::ShuangParameters(),
               list(includeEventHistories=TRUE, includeBxrecords = TRUE,
                    ##gc=0.001615, beta7=0.069848, beta8=0.224633,
                    ##g0=0.000422, gc= 0.003450, beta7=0.089526, beta8=0.158090,
                    ##beta7=0.055989, beta8=0.257808, #with frailty
                    ## beta7=0.059567, beta8=0.249589, #with frailty, new onset
                    beta7 = 0.039052, beta8 = 0.245344, # with frailty 2021-07-28
                    ## beta7=0.022520, beta8=0.223489, #without frailty
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
                                         "Formal PSA" = 21 #from NICE guideline inflated to 2021 prices                     # test sampling, primary care #callendar 2020
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
                                         + 0 * 0,                          # GP primary care #why is a GP appointment 1493?
                                         "Biopsy" = 581                           # Callendar 2021
                                         + 545,                                # Pathology of biopsy #Callender 2021
                                         "MRI" = 339,                             # MRI cost  #Callendar 2021
                                         "Combined biopsy" = 581              # Biopsy cost (TBx), still called Combined biopsy
                                         + 545,                                # Pathology of biopsy
                                         "Assessment" = 0,                      # Urologist and nurse consultation
                                         "Prostatectomy" = 9808               # Robot assisted surgery #Callendar 2021
                                         +  0*0*0                         # Radiation therapy
                                         + 0*1,                                 # Urology and nurse visit
                                         "Radiation therapy" = 6462*1          # Radiation therapy  #Callendar 2021
                                         + 0*1                               # Oncologist new visit
                                         + 0*1                               # Oncologist further visit
                                         + 0*20                                  # Nurse visit
                                         + 0*0.2,                           # Hormone therapy
                                         "Active surveillance - yearly - w/o MRI" = 5052    # Callender 2021
                                         + 0*3                                # PSA sampling
                                         + 0*3                                  # PSA analysis
                                         + 0*0.33                               # Systematic biopsy
                                         + 0*0.33,                           # Pathology of biopsy
                                         "Active surveillance - yearly - with MRI" = 5052   # Urology visit and nurse visit
                                         + 0*3                                # PSA sampling
                                         + 0*3                                  # PSA analysis
                                         + 0*0.33                               # MRI cost
                                         + 0*1.5*0.33                           # Biopsy cost (SBx|TBx)
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
                    utility_estimates = c("Invitation" = 1,                    # Heijnsdijk 2012
                                          "Formal PSA" = 1,                 # Heijnsdijk 2012
                                          "Formal panel" = 1,               # Heijnsdijk 2012
                                          "Opportunistic PSA" = 1,          # Heijnsdijk 2012
                                          "Opportunistic panel" = 1,        # Heijnsdijk 2012
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
                    ## Latest review 2019 based on Heijnsdijk 2012, Magnus 2019 and extended review by Shuang - EQ-5D
                                        # c("Invitation" = 1,                   # Heijnsdijk 2012
                                        #  "Formal PSA" = 0.99,                 # Heijnsdijk 2012
                                        #  "Formal panel" = 0.99,               # Heijnsdijk 2012
                                        #  "Opportunistic PSA" = 0.99,          # Heijnsdijk 2012
                                        #  "Opportunistic panel" = 0.99,        # Heijnsdijk 2012
                                        #  "Biopsy" = 0.90,                     # Heijnsdijk 2012
                                        #  "Combined biopsy" = 0.90,            # Heijnsdijk 2012
                                        #  "Cancer diagnosis" = 0.80,           # Heijnsdijk 2012
                                        #  "Prostatectomy part 1" = 0.829,      # Hall 2015
                                        #  "Prostatectomy part 2" = 0.893,      # Extended review by SH (Glazener 2011, Korfarge 2005, Hall 2015)
                                        #  "Radiation therapy part 1" = 0.818,  # Hall 2015
                                        #  "Radiation therapy part 2" = 0.828,  # Extended review by SH (Korfarge 2005, Hall 2015)
                                        #  "Active surveillance" = 0.9,         # Loeb 2018
                                        #  "Postrecovery period" = 0.861,       # Extended review by SH (Torvinen 2013, Waston 2016)
                                        #  "ADT+chemo" = 0.727,                 # Extended review by SH (Hall 2019, Diels 2015, Loriot 2015, Skaltsa 2014, Wu 2007, Chi 2018)
                                        #  "Palliative therapy" = 0.62,         # Magnus 2019
                                        #  "Terminal illness" = 0.40,           # Heijnsdijk 2012
                                        #  "Death" = 0.00),
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
negll = function(theta=log(c(1,80))) {
    set.seed(12345)
    sim_control <- callFhcrc(n=nsim,screen="cap_control",pop=1950,
                             parms=modifyList(Base, list(weibull_onset_shape = exp(theta[1]),
                                                         weibull_onset_scale= exp(theta[2]))),
                             mc.cores=mc.cores, print.timing=FALSE)
    pred = predict(sim_control) |>
        subset(age>=40) |>
        transform(age5=pmin(90,floor(age/5)*5)) |>
        group_by(age5) |>
        summarise(mu = sum(rate*pt)/sum(pt)) |> data.frame()
    ons2 = subset(ons, age>=40)
    val = -sum(dpois(x=ons2$cases, lambda=ons2$pop*pred$mu, log=TRUE))
    print(theta)
    print(val)
    val
}
## optim(log(c(1,80)), negll)
bobyqa(c(-0.221027947626028, 4.52299277671647), negll) # 1e4
bobyqa(c(-0.21630531253576, 4.52280852804834), negll) # 1e5
bobyqa(c(-0.206920553407683, 4.5231809059259), negll) # 1e6

plotting = function(theta=log(c(1,80)), nsim=1e4, add=FALSE, ...) {
    set.seed(12345)
    sim_control <- callFhcrc(n=nsim,screen="cap_control",pop=1950,
                             parms=modifyList(Base, list(weibull_onset_shape = exp(theta[1]),
                                                         weibull_onset_scale= exp(theta[2]))),
                             mc.cores=mc.cores, print.timing=FALSE)
    pred = predict(sim_control) |>
        subset(age>=40) |>
        transform(age5=pmin(90,floor(age/5)*5)) |>
        group_by(age5) |>
        summarise(mu = 1e5*sum(rate*pt)/sum(pt)) |> data.frame()
    pred$ons = subset(ons, age>=40)$rate
    with(pred,
         if (!add) matplot(age5, cbind(ons, mu), type="l", ...) else lines(age5, mu, ...))
    invisible(pred)
}
## df = plotting(c(-0.202809536492609, 4.52584849024777),nsim=1e6)
df = plotting(c(-0.206920553407683, 4.5231809059259),nsim=1e6)
exp(c(-0.206920553407683, 4.5231809059259)) # 0.8130842 92.1281835
data.frame(age5 = c(40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
90), mu = c(8.54696749276783, 26.1416717925314, 78.5317722700056, 
190.428065941549, 356.148652489121, 540.590847820778, 686.624474581162, 
767.779307157743, 765.177184511065, 727.733159764657, 664.553391185488
), ons = c(4, 21.3, 73.8, 188.7, 339.1, 579.4, 713.9, 810, 677.5, 
666.4, 663.2))

## Check mortality rate ratio
## Base: four-yearly screening 55-69 years (NOTE: now includes new.table)
library(prostata)
nsim=1e6
mc.cores=6
new.table <- list( background_utilities =
                       data.frame(lower=c(0, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
                                  upper=c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 1.0e55),
                                  utility=c(1, 0.934, 0.922, 0.922, 0.905, 0.905, 0.849, 0.849, 0.804, 0.804, 0.785, 0.785,
                                            0.734, 0.734))) #Updated to UK
Base <-modifyList(prostata:::ShuangParameters(),
                  list(includeEventHistories=FALSE, includeBxrecords = FALSE, indiv_reports=FALSE,
                       beta7=0.055989, beta8=0.257808, # with frailty
                       ## beta7=0.022520, beta8=0.223489, # without frailty
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
                       shapeT = 4.8, # = shapeA
                       scaleT = 15.0, # = scaleA
                       cap_pScreened = c(0.3235438, # 50-54
                                         0.3531057, # 55-59
                                         0.3582405, # 60-64
                                         0.3236445),# 65-59
                       screeningParticipation = 0.85,
                       biopsyCompliance = 0.85,
                       eol=1,
                       MRI_screen = TRUE, MRI_clinical = FALSE, MRInegSBx= FALSE, # No SBx for MRI- (by default)
                       RP_mortHR = 0.63, #(95% CI, 0.21 to 1.93) 
                       weibull_onset = TRUE,
                       weibull_onset_shape = 1,
                       weibull_onset_scale= 80,
                       discountRate.effectiveness = 0.035,
                       discountRate.costs = 0.035,
                       frailty = TRUE,
                                        #other_variance = 0,
                       cost_parameters = c("Invitation" = 0 # Invitation letter
                                           + 0, # Results letter
                                           "Formal PSA" = 21 # test sampling, primary care #NICE guideline
                                           + 0 # PSA analysis
                                           + 0 * 0, # No GP primary care
                                           "Formal panel" = 0 # test sampling, primary care + 0 # PSA analysis not included in panel price
                                           + 0 # From BergusMedical (official lab for Sthlm3)
                                           + 0 * 0, # No GP for formal
                                           "Opportunistic PSA" = 21 # PSA analysis
                                           + 0 * 0, # GP primary care
                                           "Opportunistic panel" = 0 # PSA analysis not included in panel price
                                           + 0 # From BergusMedical (official lab for Sthlm3)
                                           + 0 * 0, # GP primary care #why is a GP appointment 1493?
                                           "Biopsy" = 581 # Callendar 2021
                                           + 545, # Pathology of biopsy #Callender 2021
                                           "MRI" = 339, # MRI cost #Callendar 2021
                                           "Combined biopsy" = 581 # Biopsy cost (TBx), still called Combined biopsy
                                           + 545, # Pathology of biopsy
                                           "Assessment" = 0, # Urologist and nurse consultation
                                           "Prostatectomy" = 9808 # Robot assisted surgery #Callendar 2021
                                           + 0*0*0 # Radiation therapy
                                           + 0*1, # Urology and nurse visit
                                           "Radiation therapy" = 6462*1 # Radiation therapy #Callendar 2021
                                           + 0*1 # Oncologist new visit
                                           + 0*1 # Oncologist further visit
                                           + 0*20 # Nurse visit
                                           + 0*0.2, # Hormone therapy
                                           "Active surveillance - yearly - w/o MRI" = 5052 # Callender 2021
                                           + 0*3 # PSA sampling
                                           + 0*3 # PSA analysis
                                           + 0*0.33 # Systematic biopsy
                                           + 0*0.33, # Pathology of biopsy
                                           "Active surveillance - yearly - with MRI" = 5052 # Urology visit and nurse visit
                                           + 0*3 # PSA sampling
                                           + 0*3 # PSA analysis
                                           + 0*0.33 # MRI cost
                                           + 0*1.5*0.33 # Biopsy cost (SBx|TBx)
                                           + 0*0.33, # Pathology of biopsy
                                           "ADT+chemo" = 0, # NEW: Chemo and hormone therapy #NICE model LHRH treatment:Decapeptyl 11.25 injection (3-month dose)
                                           "Post-Tx follow-up - yearly first" = 0 # Urologist and nurse consultation
                                           + 0 # PSA test sampling
                                           + 0, # PSA analysis
                                           "Post-Tx follow-up - yearly after" = 0 # PSA test sampling
                                           + 0 # PSA analysis,
                                           + 0, # Telefollow-up by urologist
                                           "Palliative therapy - yearly" = 7383, # Palliative care cost #Callendar 2021
                                           "Terminal illness" = 7383/2, # Terminal illness cost #
                                           "Polygenic risk stratification" = 25),# Callender et al (2021) with exchange rate of approximately 12
                       background_utilities <-
                           data.frame(lower=c(0, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
                                      upper=c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 1.0e55),
                                      utility=c(1, 0.934, 0.922, 0.922, 0.905, 0.905, 0.849, 0.849, 0.804, 0.804, 0.785, 0.785,
                                                0.734, 0.734)),
                       ## Latest review 2019 based on Heijnsdijk 2012, Magnus 2019 and extended review by Shuang - PORPUS-U
                       utility_estimates = c("Invitation" = 1,                    # Heijnsdijk 2012
                                             "Formal PSA" = 1,                 # Heijnsdijk 2012
                                             "Formal panel" = 1,               # Heijnsdijk 2012
                                             "Opportunistic PSA" = 1,          # Heijnsdijk 2012
                                             "Opportunistic panel" = 1,        # Heijnsdijk 2012
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
                       pMRIposG0=0.4515496, # Pr(MRI+ | ISUP 0 || undetectable) 2020-03-23
                       pMRIposG1=0.7145305, # Pr(MRI+ | ISUP 1 && detectable) 2020-03-23
                       pMRIposG2=0.9305352, # Pr(MRI+ | ISUP 2+ && detectable) 2020-03-23
                       pSBxG0ifG1=0.1402583, # Pr(SBx gives ISUP 0 | ISUP 1) 2020-03-23
                       pSBxG0ifG2=0.1032593, # Pr(SBx gives ISUP 0 | ISUP 2) 2020-03-23
                       pSBxG1ifG2=0.119, # Pr(SBx gives ISUP 1 | ISUP 2) (not used)
                       pTBxG0ifG1_MRIpos=0.2474775, # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) #Updated 2020-03-24 from Bx -> TBx
                       pTBxG0ifG2_MRIpos=0.06570613, # Pr(TBx gives ISUP 0 | ISUP 2, MRI+) #Updated 2020-03-24 from Bx -> TBx
                       pTBxG1ifG2_MRIpos=0, # Pr(TBx gives ISUP 1 | ISUP 2, MRI+) (not used) #Updated 2020-03-24 from Bx -> TBx
                       currency_rate = 1/1, # Riksbanken 2018 #What is this?
                       risk_psa_threshold=1.5, # PSA threshold for risk-stratified screening
                       risk_lower_interval=6, # re-screening interval for lower risk
                       risk_upper_interval=4, # re-screening interval for higher risk
                       start_screening = 55.0, # start of organised screening
                       stop_screening = 69.0, # end of organised screening
                       screening_interval = 4.0,
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
                             ) # 2017-2019 death, rates from UK life tables,
                       ))
set.seed(12345)
sim_control <- callFhcrc(n=nsim,screen="cap_control",pop=1950, parms=Base, mc.cores=mc.cores,
                         tables=c(new.table,list(rescreening=uk_rescreening,stockholmTreatment = stockholmTreatment)))
noScreen <- callFhcrc(n=nsim,screen="noScreen",pop=1950, parms=Base, mc.cores=mc.cores,
                      tables=c(new.table,list(rescreening=uk_rescreening,stockholmTreatment = stockholmTreatment)))
sim4_5570 = callFhcrc(n=nsim,screen="regular_screen", pop=1950, mc.cores=mc.cores,
                      parms=Base,
                      tables=c(new.table,list(rescreening=uk_rescreening, stockholmTreatment = stockholmTreatment)))
sim4_6070 = callFhcrc(n=nsim,screen="regular_screen", pop=1950, mc.cores=mc.cores,
                      parms=modifyList(Base,list(start_screening=60)),
                      tables=c(new.table,list(rescreening=uk_rescreening, stockholmTreatment = stockholmTreatment)))
sim50 = callFhcrc(n=nsim,screen="screen50", pop=1950, mc.cores=mc.cores,
                                parms=Base,
                                tables=c(new.table,list(rescreening=uk_rescreening, stockholmTreatment = stockholmTreatment)))
sim60 = callFhcrc(n=nsim,screen="screen60", pop=1950, mc.cores=mc.cores,
                                parms=Base,
                                tables=c(new.table,list(rescreening=uk_rescreening, stockholmTreatment = stockholmTreatment)))
##
## Mortality rates
## plot(rate~age, predict(noScreen, type="pc.mortality.rate"), type="l")
## lines(rate~age, predict(sim_control, type="pc.mortality.rate"), col="blue") ## ok - // noScreen
## lines(rate~age, predict(sim4_5570, type="pc.mortality.rate"), col="red")
cumhaz = function(object,ages)
    sum(subset(predict(object, type="pc.mortality.rate"), age %in% ages)$rate)
mrr = function(object1,object2,ages) cumhaz(object2,ages) /cumhaz(object1,ages)
mrr(noScreen,sim4_5570,1+55:(55+15)) # 0.80
mrr(noScreen,sim4_6070,1+60:(60+15)) # 0.83
mrr(noScreen,sim50,1+50:(50+9))      # 0.90
mrr(noScreen,sim60,1+60:(60+9))      # 0.91

## Edna extensions
## Check coding for Weibull distribution with a frailty
library(prostata)
frailty = 0.1
rweibullFrailty1 = function(n,shape,scale,frailty) {
    scale*rexp(n,frailty)^(1/shape)
}
set.seed(123456)
rweibullFrailty2 = function(n,shape,scale,frailty) {
    rweibull(n,shape,scale/frailty^(1/shape))
}
set.seed(123456)
y1=rweibullFrailty1(1e5,2,3,frailty)
set.seed(123456)
y2=rweibullFrailty2(1e5,2,3,frailty)
plot(density(y1,from=0))
lines(density(y2,from=0),lty=2)
## Compare frailty distributions with the previous onset distribution
g0 = 0.0005
function(age) pmax(0,age-35)*g0
H = function(age) pmax(0,age-35)^2*g0/2
H2 = Vectorize(function(age) integrate(h,0,age)$value)
par(mfrow=c(2,2))
x = seq(0,100,length=301)
## survival plot
plot(x,exp(-H(x)),type="l",ylab="Survival",ylim=0:1)
lines(x,exp(-H2(x)),col=2)
lines(35+x,pweibull(x,shape=2,scale=sqrt(2/g0),lower.tail=FALSE),col=3)
## density plot
plot(x,h(x)*exp(-H(x)),type="l",ylab="Density")
lines(x,h(x)*exp(-H2(x)),col=2)
lines(35+x,dweibull(x,shape=2,scale=sqrt(2/g0)),col=3)
## hazard plot
plot(x,h(x),type="l",ylab="Hazard")
lines(x,h(x),col=2)
lines(35+x,
      dweibull(x,shape=2,scale=sqrt(2/g0))/pweibull(x,shape=2,scale=sqrt(2/g0),lower.tail=FALSE),
      col=3)
##
set.seed(12345)
n = 1e5
x = seq(0,100,length=301)
g0 = 0.0005
grs_variance = 0.68           # Callender et al (2019)
other_variance = 1.14          # total variance = 1.82 from Kicinski et al (2011)
grs_log_mean = -grs_variance/2
grs_log_sd = sqrt(grs_variance)
other_log_mean = -other_variance/2
other_log_sd = sqrt(other_variance)
grs_frailty = rlnorm(n, grs_log_mean, grs_log_sd)
other_frailty = rlnorm(n, other_log_mean, other_log_sd)
plot(density(pmin(grs_frailty*other_frailty,10),from=0))
rweibullFrailty = function(n,shape,scale,frailty) {
    rweibull(n,shape,scale/frailty^(1/shape))
}
## base shape=1
plot(35+x,dweibull(x,shape=2,scale=sqrt(2/g0)), type="l")
y = pmin(200,35+rweibullFrailty(n, shape=1, scale=sqrt(2/g0)*1, frailty=grs_frailty*other_frailty))
lines(density(y,from=35),col=3)
## base shape=1.5
plot(35+x,dweibull(x,shape=2,scale=sqrt(2/g0)), type="l")
y = pmin(200,35+rweibullFrailty(n, shape=1.5, scale=sqrt(2/g0)*0.4, frailty=grs_frailty*other_frailty))
lines(density(y,from=35),col=4)
## base shape=2
plot(35+x,dweibull(x,shape=2,scale=sqrt(2/g0)), type="l")
y = pmin(200,35+rweibullFrailty(n, shape=2, scale=sqrt(2/g0)*0.5, frailty=grs_frailty*other_frailty))
lines(density(y,from=35),col=4)
## base shape=3
plot(35+x,dweibull(x,shape=2,scale=sqrt(2/g0)), type="l")
y = pmin(200,35+rweibullFrailty(n, shape=3, scale=sqrt(2/g0)*0.75, frailty=grs_frailty*other_frailty))
lines(density(y,from=35),col=4)

plot(35+x,dweibull(x,shape=1,scale=120), type="l")
y = pmin(200,35+rweibullFrailty(n, shape=1.5, scale=120*0.5, frailty=grs_frailty*other_frailty))
lines(density(y,from=35),col=4)


## STHLM3-MRI simulation
library(prostata)
library(dplyr)
## prostata::sthlm3_mri_arm
## thresholds1.5 = c(5.322, 8.933, 3.427)
s3m_parms1.5 = modifyList(prostata:::ShuangParameters(year=2019),
                       list(includePSArecords = FALSE,
                            includeEventHistories = FALSE,
                            formal_compliance=1,
                            pMRIposG0=0.167,          # Pr(MRI+ | ISUP 0 || undetectable)
                            pMRIposG1=0.96,     # Pr(MRI+ | ISUP 1 && detectable)
                            pMRIposG2=0.96,     # Pr(MRI+ | ISUP 2+ && detectable)
                            MRI_screen = TRUE,
                            PSA_FP_threshold_nCa=5.322, # reduce FP in no cancers with PSA threshold
                            PSA_FP_threshold_GG6=8.933, # reduce FP in GG 6 with PSA threshold
                            PSA_FP_threshold_GG7plus=3.427, # reduce FP in GG >= 7 with PSA threshold
                            panelReflexThreshold = 1.5,
                            SplitS3M25plus = FALSE,
                            start_screening=55,
                            stop_screening=70,
                            screening_interval=4))
psa_parms = modifyList(s3m_parms1.5,
                   list(pMRIposG0=0.148,     # Pr(MRI+ | ISUP 0 || undetectable)
                        pMRIposG1=0.743,     # Pr(MRI+ | ISUP 1 && detectable)
                        pMRIposG2=0.948,     # Pr(MRI+ | ISUP 2+ && detectable)
                        SplitS3M25plus = FALSE,
                        PSA_FP_threshold_nCa=3, # reduce FP in no cancers with PSA threshold
                        PSA_FP_threshold_GG6=3, # reduce FP in GG 6 with PSA threshold
                        PSA_FP_threshold_GG7plus=3 # reduce FP in GG >= 7 with PSA threshold
                        ))
## PSA>=2.0
## c(7.40, 7.94, 4.74) # pseudo-thresholds
## c(0.5428571, 0.8000000, 0.9597701) # =1/adjBoth2.0 -- the split proportions
## adjBoth2.0 = c(1.84210526315789,1.25,1.04191616766467)
## PrMRIpos_s3m2.0 = c(0.159,0.768,0.950)
## adjBoth2.0*PrMRIpos_s3m2.0
## c(0.2928947, 0.96, 0.9898204)
thresholds2.0 = c(6.123, 9.901, 4.455)
s3m_parms2.0 = modifyList(s3m_parms1.5,
                       list(pMRIposG0=0.164,
                            pMRIposG1=0.96,
                            pMRIposG2=0.959,
                            PSA_FP_threshold_nCa=6.123,
                            PSA_FP_threshold_GG6=9.901,
                            PSA_FP_threshold_GG7plus=4.455,
                            panelReflexThreshold = 2.0))
thresholds11_2.0 = c(3.311,5.813,1.167)
PrMRIpos11_s3m2.0 = c(0.149,0.968,0.962)
s3m11_parms2.0 = modifyList(s3m_parms1.5,
                            list(pMRIposG0=0.149,
                                 pMRIposG1=0.968,
                                 pMRIposG2=0.962,
                                 PSA_FP_threshold_nCa=3.311,
                                 PSA_FP_threshold_GG6=5.813,
                                 PSA_FP_threshold_GG7plus=1.167,
                                 panelReflexThreshold = 2.0))
thresholds11_2.5 = c(4.230,8.612,3.165)
PrMRIpos11_s3m2.5 = c(0.161,0.968,0.961)
s3m11_parms2.5 = modifyList(s3m_parms1.5,
                            list(pMRIposG0=0.161,
                                 pMRIposG1=0.968,
                                 pMRIposG2=0.961,
                                 PSA_FP_threshold_nCa=4.230,
                                 PSA_FP_threshold_GG6=8.612,
                                 PSA_FP_threshold_GG7plus=3.165,
                                 panelReflexThreshold = 2.0))
s3m_parms2.0_nomri = modifyList(s3m_parms1.5,
                                list(PSA_FP_threshold_nCa = 4.56,
                                     PSA_FP_threshold_GG6 = 4.31,
                                     PSA_FP_threshold_GG7plus = 3.28,
                                     panelReflexThreshold = 2.0,
                                     MRI_screen = FALSE))
thresholds2.5 = c(6.696, 13.603, 5.250)
PrMRIpos_s3m2.5 = c(0.164,0.96,0.959)
s3m_parms2.5 = modifyList(s3m_parms1.5,
                       list(pMRIposG0=0.164,
                            pMRIposG1=0.96,
                            pMRIposG2=0.959,
                            PSA_FP_threshold_nCa=6.696,
                            PSA_FP_threshold_GG6=13.603,
                            PSA_FP_threshold_GG7plus=5.250,
                            panelReflexThreshold = 2.5))
## ## STHLM3-MRI simulation (OLD!)
## s3m_parms = modifyList(prostata:::ShuangParameters(),
##                        list(includePSArecords = FALSE,
##                             includeEventHistories = FALSE,
##                             formal_compliance=1,
##                             pMRIposG0=0.28,          # Pr(MRI+ or MRI-S3m>=25 | ISUP 0 || undetectable) 2020-03-23
##                             pMRIposG1=0.9595385,     # Pr(MRI+ or MRI-S3m>=25 | ISUP 1 && detectable) 2020-03-23
##                             pMRIposG2=0.9900455,     # Pr(MRI+ or MRI-S3m>=25 | ISUP 2+ && detectable) 2020-03-23
##                             MRI_screen = TRUE,
##                             PSA_FP_threshold_nCa=6.38, # reduce FP in no cancers with PSA threshold
##                             PSA_FP_threshold_GG6=7.04, # reduce FP in GG 6 with PSA threshold
##                             PSA_FP_threshold_GG7plus=3.55, # reduce FP in GG >= 7 with PSA threshold
##                             panelReflexThreshold = 1.5,
##                             SplitS3M25plus = TRUE,
##                             PrS3MposIfBx_nCa = 0.575,
##                             PrS3MposIfBx_GG6 = 0.7878788,
##                             PrS3MposIfBx_GG7plus = 0.9565217,
##                             start_screening=55,
##                             stop_screening=70,
##                             screening_interval=4))
## psa_parms = modifyList(s3m_parms,
##                    list(pMRIposG0=0.148,          # Pr(MRI+ or MRI-S3m>=25 | ISUP 0 || undetectable) 2020-03-23
##                         pMRIposG1=0.743,     # Pr(MRI+ or MRI-S3m>=25 | ISUP 1 && detectable) 2020-03-23
##                         pMRIposG2=0.948,     # Pr(MRI+ or MRI-S3m>=25 | ISUP 2+ && detectable) 2020-03-23
##                         SplitS3M25plus = FALSE,
##                         PSA_FP_threshold_nCa=3, # reduce FP in no cancers with PSA threshold
##                         PSA_FP_threshold_GG6=3, # reduce FP in GG 6 with PSA threshold
##                         PSA_FP_threshold_GG7plus=3 # reduce FP in GG >= 7 with PSA threshold
##                         ))
## ## PSA>=2.0
## ## c(7.40, 7.94, 4.74) # pseudo-thresholds
## ## c(0.5428571, 0.8000000, 0.9597701) # =1/adjBoth2.0 -- the split proportions
## ## adjBoth2.0 = c(1.84210526315789,1.25,1.04191616766467)
## ## PrMRIpos_s3m2.0 = c(0.159,0.768,0.950)
## ## adjBoth2.0*PrMRIpos_s3m2.0
## ## c(0.2928947, 0.96, 0.9898204)
## s3m_parms2.0 = modifyList(s3m_parms,
##                        list(pMRIposG0=0.2928947,
##                             pMRIposG1=0.96,
##                             pMRIposG2=0.9898204,
##                             PSA_FP_threshold_nCa=7.4,
##                             PSA_FP_threshold_GG6=7.94,
##                             PSA_FP_threshold_GG7plus=4.74,
##                             panelReflexThreshold = 2.0,
##                             PrS3MposIfBx_nCa = 0.5428571,
##                             PrS3MposIfBx_GG6 = 0.8,
##                             PrS3MposIfBx_GG7plus = 0.9597701))

n.sim <- 1e6
noscreen <- callFhcrc(n.sim, screen="noScreening", pop=1940, flatPop=TRUE, parms=s3m_parms1.5,
                   mc.cores=6, panel=TRUE)
s3m1.5 <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms1.5,
                   mc.cores=6, panel=TRUE)
psa <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=psa_parms,
                   mc.cores=6, panel=FALSE)
s3m2.0 <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms2.0,
                    mc.cores=6, panel=TRUE)
s3m11_2.0 <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m11_parms2.0,
                    mc.cores=6, panel=TRUE)
s3m11_2.5 <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m11_parms2.5,
                    mc.cores=6, panel=TRUE)
s3m_2.0_nomri <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms2.0_nomri,
                    mc.cores=6, panel=TRUE)
s3m2.5 <- callFhcrc(n.sim, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms2.5,
                    mc.cores=6, panel=TRUE)
ICER(s3m1.5,psa,from=55)
ICER(s3m2.0,psa,from=55)
ICER(s3m1.5,s3m2.0,from=55)
ICER(s3m11_2.0,s3m2.0,from=55)
ICER(psa,noscreen,from=55)
ICER(s3m11_2.0,noscreen,from=55)
ICER(s3m_2.0_nomri,noscreen,from=55)
## save(s3m1.5, psa, s3m2.0, file="~/src/R/prostata/test/sims-20210128.RData")
## load("~/src/R/prostata/test/sims-20210128.RData")

d = rbind(with(summary(psa,from=55), c(QALE,societal.costs)),
          with(summary(s3m1.5,from=55), c(QALE,societal.costs)),
          with(summary(s3m2.0,from=55), c(QALE,societal.costs)),
          with(summary(s3m2.5,from=55), c(QALE,societal.costs)),
          with(summary(s3m11_2.0,from=55), c(QALE,societal.costs)),
          with(summary(s3m11_2.5,from=55), c(QALE,societal.costs)),
          with(summary(s3m_2.0_nomri,from=55), c(QALE,societal.costs)),
          with(summary(noscreen,from=55), c(QALE,societal.costs))
          )
rownames(d)=c("PSA+MRI","S3M+MRI (15%, 1.5)", "S3M+MRI (15%, 2.0)",
              "S3M+MRI (15%, 2.5)",
              "S3M+MRI (11%, 2.0)", "S3M+MRI (11%, 2.5)",
              "S3M+SBx (11%, 2.0)", "No screening")
eps=c(3e-4,1e-2)
plot(d, xlab="QALE",ylab="E(costs)",
     xlim=range(d[,1])*c(1-eps[1],1+eps[1]),
     ylim=range(d[,2])*c(1-eps[2],1+eps[2]))
lines(d[c("No screening","S3M+MRI (15%, 2.5)","PSA+MRI","S3M+MRI (11%, 2.0)"),], lty=2)
iord = order(d[,2])
text(d[iord,],labels=rownames(d[iord,]),pos=c(2,4,4,2,2,2,4,2,2))

s3m1.5 <- callFhcrc(1e6, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms,
                   mc.cores=6, nLifeHistories=1e8, panel=TRUE)
psa <- callFhcrc(1e6, screen="regular_screen", pop=1940, flatPop=TRUE, parms=psa_parms,
                   mc.cores=6, nLifeHistories=1e8, panel=FALSE)
s3m2.0 <- callFhcrc(1e6, screen="regular_screen", pop=1940, flatPop=TRUE, parms=s3m_parms2.0,
                    mc.cores=6, nLifeHistories=1e8, panel=TRUE)
ICER(s3m1.5,psa,from=55)
ICER(s3m2.0,psa,from=55)
ICER(s3m1.5,s3m2.0,from=55)

merge(as.data.frame(xtabs(costs~item,s3m1.5$healthsector.costs)),
      as.data.frame(xtabs(costs~item,psa$healthsector.costs)), by="item",all=TRUE)


## Code to compare the simulation with the trial
s3m <- callFhcrc(1e5, screen="sthlm3_mri_arm", pop=sthlm3_mri_arm,
                 parms=modifyList(s3m_parms1.5, list(includePSArecords=TRUE,
                                                  includeEventHistories=TRUE)),
                 mc.cores=6, nLifeHistories=1e8, panel=TRUE)
s3m2.0 <- callFhcrc(1e5, screen="sthlm3_mri_arm", pop=sthlm3_mri_arm,
                 parms=modifyList(s3m_parms2.0, list(includePSArecords=TRUE,
                                                  includeEventHistories=TRUE)),
                 mc.cores=6, nLifeHistories=1e8, panel=TRUE)
psa <- callFhcrc(1e5, screen="sthlm3_mri_arm", pop=sthlm3_mri_arm,
                 parms=modifyList(psa_parms, list(includePSArecords=TRUE,
                                                  includeEventHistories=TRUE)),
                 mc.cores=6, nLifeHistories=1e8, panel=FALSE)
table(s3m_entry$isup)
table(psa_entry$isup) # not identical??
table(s3m$lifeHistories$event)

## join entry with other events
ISUP = function(x) ifelse(x==0,1,
                   ifelse(x %in% c(1,2), 2,
                   ifelse(x==3,0,
                          NA)))
report <- function(data,type=c("parameters","entry","events")) {
    type <- match.arg(type)
    parameters <- filter(data$parameters, (age_pca==-1 | age_pca>ageEntry) & age_d>ageEntry)
    if (type %in% c("entry","events"))
        entry = inner_join(parameters, data$psarecord, by="id") %>% filter(age==ageEntry) %>%
            mutate(isup=ISUP(ifelse(detectable==1,future_ext_grade,3)))
    if (type=="events") {
        events = inner_join(entry, data$lifeHistories, by="id")
        eps=1e-8
        events = filter(events, (ageEntry+1/52-eps<=end & end<=ageEntry+1/52+eps) |
                              (ageEntry+4/52-eps<=end & end<=ageEntry+4/52+eps))
    }
    switch(type,parameters=parameters,entry=entry,events=events)
}
temp2 = report(s3m,type="entry")
mean(temp2$psa>=1.5)
length(unique(temp2$id))

BaseReport = function(data) {
    temp2 = report(data,type="events")
    temp2 %>%
        filter(event %in% c("toScreenInitiatedBiopsy","toScreenDiagnosis","toMRI")) %>%
        "[["("event") %>%
        table %>% as.data.frame %>% filter(Freq!=0) %>%
        mutate(p=Freq/length(unique(temp2$id)))
}
BaseReport(psa)
BaseReport(s3m)
BaseReport(s3m2.0)


BxReport = function(data) {
    events = report(data,type="events")
    tab = table(subset(events,event=="toScreenDiagnosis")$isup) %>%
        as.data.frame
    tab$Freq[1] = nrow(subset(events,event=="toScreenInitiatedBiopsy")) -
        sum(tab$Freq[2:3])
    tab = mutate(tab,p = Freq/length(unique(events$id)))
    tab
}
BxReport(psa)
BxReport(s3m)
BxReport(s3m2.0)

table(subset(report(s3m,type="events"),event=="toScreenInitiatedBiopsy")$isup) /
    table(subset(report(psa,type="events"),event=="toScreenDiagnosis")$isup)
table(subset(report(s3m2.0,type="events"),event=="toScreenDiagnosis")$isup) /
    table(subset(report(psa,type="events"),event=="toScreenDiagnosis")$isup)


## Phase II
## Denominator: Expected number of men in the MRI arm
## total*(men randomised to MRI)/(men randomised)
##
## Expected number of men S3M+
## Expected number of men S3M+ with MRI performed
## Expected number of men S3M+MRI+
## Expected number of men MRI-S3M>=25
## Expected number of men (S3M+MRI+ or MRI-S3M>=25) by GG (=Bx)
##
## Expected number of men PSA+
## Expected number of men PSA+MRI+
## Expected number of men PSA+MRI+ by GG (=Bx)

sapply(names(s3m$parameters), function(nm) sum(s3m$parameters[[nm]]!=psa$parameters[[nm]]))
head(which(s3m$parameters$age_pca != psa$parameters$age_pca)-1)

options(width=140)
a=subset(s3m$lifeHistories, id==65984)
b=subset(psa$lifeHistories, id==65984)
sapply(names(a), function(nm) sum(!peq(a[[nm]],b[[nm]])))
cbind(a$state,b$state)
cbind(a$event,b$event)
a[23,]
b[23,]

plot(table(table(subset(s3m$lifeHistories,event=="toMRI")$id)))
table(subset(s3m$lifeHistories,event=="toMRI")$id) %>% as.data.frame %>% filter(Freq>=10) %>% nrow
subset(s3m$lifeHistories,id==65984)

dim(subset(s3m$lifeHistories, id==30))
dim(subset(psa$lifeHistories, id==30))


## look at opportunistic rescreening
myplot = function(data,age) {
    x=seq(0,20,length=1001)
    myline = function(x, row, ...) lines(x, pweibull(x,row$shape,row$scale,lower.tail=TRUE)*(1-row$cure), ...)
    plot(c(0,20),0:1,type="n",xlab="Time (years)",ylab="Pr(Rescreened)",main=age)
    for(i in 1:nrow(data))
        myline(x, data[i,],col=i,lty=1)
}
par(mfrow=c(4,4))
for (.age5 in unique(rescreening$age5))
    myplot(subset(rescreening,age5==.age5),.age5)
plot(0:1,0:1,type="n", axes=FALSE)
legend("topleft",legend=c("0-","1-","3-","10-"),lty=1,col=1:4)
## rescreening is steep 

oldpar = options()
options(width=140)
subset(s3m$lifeHistories,id==4)
options(width=80)


## STHLM3-MRI test calibration
library(prostata)
library(dplyr)
parms = modifyList(prostata:::ShuangParameters(),
                   list(includePSArecords = TRUE, includeLifeHistories=TRUE))
model <- callFhcrc(1e6, screen="sthlm3_mri_arm", pop=sthlm3_mri_arm, parms=parms,
                   mc.cores=6, nLifeHistories=1e8)
## remove men who got cancer or died before "ageEntry"
parameters <- filter(model$parameters, (age_pca==-1 | age_pca>ageEntry) & age_d>ageEntry)
## get PSA records at age of study entry
ISUP = function(x) ifelse(x==0,1,
                   ifelse(x %in% c(1,2), 2,
                   ifelse(x==3,0,
                          NA)))
temp = inner_join(parameters, model$psarecord, by="id") %>% filter(age==ageEntry) %>%
    mutate(isup=ISUP(ifelse(detectable==1,future_ext_grade,3)))
## ## ratio of (MRI-S3M>=25 or MRI+S3M+) to MRI+S3M+ --> number of biopsies
## adjBoth25.1.5 = c(1.73913043478261,1.26923076923077,1.04545454545455) 
## adjBoth25.2.0 = c(1.84210526315789,1.25,1.04191616766467)
## see Table 1A
PrMRIpos_s3m1.5 = c(0.167,0.96,0.96)
PrMRIpos_s3m2.0 = c(0.164,0.96,0.959)
PrMRIpos11_s3m1.5 = c(0.145,0.968,0.962)
PrMRIpos11_s3m2.0 = c(0.149,0.968,0.962)
PrMRIpos11_s3m2.5 = c(0.161,0.968,0.961)
PrMRIpos_s3m2.5 = c(0.164,0.96,0.959)
PrMRIpos_psa = c(0.148,0.743,0.948)
##
PrMRIpos_s3m1.5_tbx = c(0.224,0.96,0.96)
PrMRIpos_s3m2.0_tbx = c(0.223,0.96,0.959)
PrMRIpos_psa_tbx = c(0.175,0.743,0.948)
PrMRIpos_s3m1.5_ITT = c(0.175,0.960,0.960)
PrMRIpos_s3m2.0_ITT = c(0.174,0.960,0.959)
PrMRIpos_psa_ITT = c(0.150,0.743,0.948)
##
rpf1.5 = c(0.597, 0.743, 0.994)
rpf2.0 = c(0.494,0.686,0.944)
rpf11_2.0 = c(0.909,1,1.09)
rpf11_2.5 = c(0.753,0.771,1.006)
rpf2.5 = c(0.442,0.514,0.904)
rpf11_1.5=c(1.078,1.171,1.192)
##
rpf1.5_tbx = c(0.63,1.174,0.975)
rpf2.0_tbx = c(0.528,1.087,0.938)
rpf1.5_ITT = c(0.81,0.829,1)
rpf2.0_ITT = c(0.714,0.756,0.948)

## PrMRIpos_s3m1.5 = c(0.161,0.756,0.947)
## PrMRIpos_s3m2.0 = c(0.159,0.768,0.950)
## PrMRIpos_psa = c(0.148,0.743,0.948)
## rpf25.1.5 = c(0.8, 0.825, 1.005)
## rpf25.2.0 = c(0.7,0.75,0.951)
## rpf1.5 = c(0.597, 0.743, 0.994)
## rpf2.0 = c(0.494,0.686,0.944)
## expected number of MRI+ (== number of biopsies if 100% compliance == number of diagnoses if 100% accuracy) for PSA>=3
temp1 = temp %>% filter(psa>=3) %>% mutate(mripos=PrMRIpos_psa[isup+1])
EmriPos = temp1 %>% group_by(isup) %>%
    summarise(Emri=sum(mripos), .groups="drop_last") %>% data.frame

## assumes four-yearly rescreening
do = function(PrMRIpos_s3m,rpf,adjBoth=rep(1,3)) {
    temp2 <- temp %>% mutate(mripos=PrMRIpos_s3m[isup+1], bx=mripos*adjBoth[isup+1]) 
    f = function(data,.isup,.psa,EmriPosData)
        sum(filter(data,isup==.isup & psa>=.psa)$bx)/EmriPosData$Emri[.isup+1]
    list(uniroot(function(psa) f(temp2,0,psa,EmriPos)-rpf[1], c(0.1,100)),
         uniroot(function(psa) f(temp2,1,psa,EmriPos)-rpf[2], c(0.1,100)),
         uniroot(function(psa) f(temp2,2,psa,EmriPos)-rpf[3], c(0.1,100)))
}
## Do.call = function(list,f) do.call(f,list)
do.call(rbind,do(PrMRIpos_s3m1.5,rpf1.5))
do.call(rbind,do(PrMRIpos_s3m2.0,rpf2.0))
do.call(rbind,do(PrMRIpos11_s3m1.5,rpf11_1.5))
do.call(rbind,do(PrMRIpos11_s3m2.0,rpf11_2.0))
do.call(rbind,do(PrMRIpos11_s3m2.5,rpf11_2.5))
do.call(rbind,do(PrMRIpos_s3m2.5,rpf2.5))
## do.call(rbind,do(PrMRIpos_s3m1.5,rpf25.1.5,adjBoth25.1.5))
## do.call(rbind,do(PrMRIpos_s3m2.0,rpf25.2.0,adjBoth25.2.0))
do.call(rbind,do(PrMRIpos_s3m1.5_tbx,rpf1.5_tbx))
do.call(rbind,do(PrMRIpos_s3m2.0_tbx,rpf2.0_tbx))
do.call(rbind,do(PrMRIpos_s3m1.5_ITT,rpf1.5_ITT))
do.call(rbind,do(PrMRIpos_s3m2.0_ITT,rpf2.0_ITT))

## given this pseudo threshold: (i) test if MRI+; (ii) if MRI-, also test if S3M>=25
## Alternatively: test if bx; if bx, split by whether MRI+ or MRI-/S3M>=25
getRoots = function(x) dput(round(sapply(x, "[[", "root"),3))
thresholds1.5 = getRoots(do(PrMRIpos_s3m1.5,rpf1.5)) # c(5.322, 8.933, 3.427)
thresholds2.0 = sapply(do(PrMRIpos_s3m2.0,rpf2.0), "[[", "root") # c(6.123, 9.901, 4.455)
thresholds11_2.0 = sapply(do(PrMRIpos11_s3m2.0,rpf11_2.0), "[[", "root") # c(3.311,5.813,1.167)
thresholds11_2.5 = sapply(do(PrMRIpos11_s3m2.5,rpf11_2.5), "[[", "root") # c(4.230,8.612,3.165)
thresholds2.5 = sapply(do(PrMRIpos_s3m2.5,rpf2.5), "[[", "root") # c(6.696, 13.603, 5.250)
getRoots(do(PrMRIpos_s3m1.5_tbx,rpf1.5_tbx)) # c(6.475, 4.021, 3.836)
getRoots(do(PrMRIpos_s3m2.0_tbx,rpf2.0_tbx)) # c(7.414, 4.841, 4.593)
getRoots(do(PrMRIpos_s3m1.5_ITT,rpf1.5_ITT)) # c(4.268, 7.715, 3.282)
getRoots(do(PrMRIpos_s3m2.0_ITT,rpf2.0_ITT)) # c(4.74, 8.759, 4.391)
## thresholds25.1.5 = sapply(do(PrMRIpos_s3m1.5,rpf25.1.5,adjBoth25.1.5), "[[", "root")
## thresholds25.2.0 = sapply(do(PrMRIpos_s3m2.0,rpf25.2.0,adjBoth25.2.0), "[[", "root")

do = function(thresholds,PrMRIpos_s3m,rpf,adjBoth=rep(1,3)) {
    temp2 <- temp %>% mutate(mripos=PrMRIpos_s3m[isup+1], bx=mripos*adjBoth[isup+1]) 
    ## ratio of MRI tests (S3M+ vs PSA+)
    out = data.frame(ratioMRI=c(nrow(subset(temp2,psa>=thresholds[1] & isup==0))/nrow(subset(temp2,psa>=3 & isup==0)),
                                nrow(subset(temp,psa>=thresholds[2] & isup==1))/nrow(subset(temp,psa>=3 & isup==1)),
                                nrow(subset(temp,psa>=thresholds[3] & isup==2))/nrow(subset(temp,psa>=3 & isup==2))),
                     ## ratio of biopsies
                     ratioBiopsies=c(sum(subset(temp2,psa>=thresholds[1] & isup==0)$bx)/sum(subset(temp1,psa>=3 & isup==0)$mripos),
                                     sum(subset(temp2,psa>=thresholds[2] & isup==1)$bx)/sum(subset(temp1,psa>=3 & isup==1)$mripos),
                                     sum(subset(temp2,psa>=thresholds[3] & isup==2)$bx)/sum(subset(temp1,psa>=3 & isup==2)$mripos)))
    rownames(out)=c("GG=0","GG=1","GG>=2")
    out
}
do(thresholds1.5,PrMRIpos_s3m1.5,rpf1.5)
do(thresholds2.0,PrMRIpos_s3m2.0,rpf2.0)
do(thresholds25.1.5,PrMRIpos_s3m1.5,rpf25.1.5,adjBoth25.1.5)
do(thresholds25.2.0,PrMRIpos_s3m2.0,rpf25.2.0,adjBoth25.2.0)

## simulate for thresholds
library(prostata)
library(dplyr)
parms = modifyList(prostata:::ShuangParameters(),
                   list(includePSArecords = TRUE, includeLifeHistories=TRUE))
model <- callFhcrc(1e6, screen="sthlm3_mri_arm", pop=sthlm3_mri_arm, parms=parms,
                   mc.cores=6, nLifeHistories=1e8)
## remove men who got cancer or died before "ageEntry"
parameters <- filter(model$parameters, (age_pca==-1 | age_pca>ageEntry) & age_d>ageEntry)
## get PSA records at age of study entry
ISUP = function(x) ifelse(x==0,1,
                   ifelse(x %in% c(1,2), 2,
                   ifelse(x==3,0,
                          NA)))
temp = inner_join(parameters, model$psarecord, by="id") %>% filter(age==ageEntry) %>%
    mutate(isup=ISUP(ifelse(detectable==1,future_ext_grade,3)))
logitInputs=list(PrMRIpos_psaGG0 = c(p=0.148, lower=0.126, upper=0.192),
                 PrMRIpos_psaGG1 = c(0.743, 0.676, 0.816),
                 PrMRIpos_psaGG2 = c(0.948, 0.925, 0.971),
                 PrMRIpos_s3m1.5GG0 = c(0.167, 0.124, 0.224),
                 PrMRIpos_s3m1.5GG1 = c(0.96, 0.796, 0.999),
                 PrMRIpos_s3m1.5GG2 = c(0.96, 0.9, 0.989),
                 PrMRIpos_s3mGG0 = c(0.164, 0.119, 0.226),
                 PrMRIpos_s3mGG1 = c(0.960, 0.796, 0.999),
                 PrMRIpos_s3mGG2 = c(0.959, 0.899, 0.989))
logInputs=list(rpf1.5GG0 = c(0.597, 0.454, 0.785),
               rpf1.5GG1 = c(0.743, 0.524, 1.054),
               rpf1.5GG2 = c(0.994, 0.91, 1.086),
               rpf2.0GG0 = c(0.494, 0.372, 0.655),
               rpf2.0GG1 = c(0.686, 0.483, 0.974),
               rpf2.0GG2 = c(0.944, 0.868, 1.026))
logit=binomial()$linkfun
expit=binomial()$linkinv
sapply(logitInputs, function(tuple) tuple[1] - expit(mean(logit(tuple[2:3])))) # check
sapply(logInputs, function(tuple) tuple[1] - exp(mean(log(tuple[2:3])))) # check
set.seed(1234567)
n = 1000
logitSims = lapply(logitInputs, function(tuple) {
    se = diff(logit(tuple[2:3]))/2/1.96
    expit(rnorm(n, logit(tuple[1]), se))
})
logSims = lapply(logInputs, function(tuple) {
    se = diff(log(tuple[2:3]))/2/1.96
    exp(rnorm(n, log(tuple[1]), se))
})
## t(sapply(logitSims, function(x) c(mean=mean(x), quantile(x,c(0.025, 0.975)))))
## t(sapply(logSims, function(x) c(mean=mean(x), quantile(x,c(0.025, 0.975)))))
PrMRIpos_psa = with(logitSims, cbind(PrMRIpos_psaGG0, PrMRIpos_psaGG1, PrMRIpos_psaGG2))
PrMRIpos_s3m = with(logitSims, cbind(PrMRIpos_s3mGG0, PrMRIpos_s3mGG1, PrMRIpos_s3mGG2))
PrMRIpos_s3m1.5 = with(logitSims, cbind(PrMRIpos_s3m1.5GG0, PrMRIpos_s3m1.5GG1, PrMRIpos_s3m1.5GG2))
rpf = with(logSims, cbind(rpf2.0GG0, rpf2.0GG1, rpf2.0GG2))
rpf1.5 = with(logSims, cbind(rpf1.5GG0, rpf1.5GG1, rpf1.5GG2))
nMRI = sapply(0:2,function(.isup) nrow(filter(temp,psa>=3 & isup==.isup)))
PSAs = sapply(0:2, function(.isup) sort(filter(temp, isup==.isup)$psa))
do = function(.isup,PrMRIpos_psa,PrMRIpos_s3m,rpf) {
    psa = PSAs[[.isup+1]]
    index = findInterval(PrMRIpos_psa*nMRI[.isup+1]*rpf, (1:length(psa))*PrMRIpos_s3m)
    rev(psa)[index]
}
do2 = function(.isup,PrMRIpos_psa,PrMRIpos_s3m,rpf) {
    EmriPos = PrMRIpos_psa*nMRI[.isup+1]
    temp2 = data.frame(psa=PSAs[[.isup+1]]) %>% mutate(bx=PrMRIpos_s3m)
    f = function(.psa)
        sum(filter(temp2,psa>=.psa)$bx)/EmriPos
    uniroot(function(psa) f(psa)-rpf, c(0.1,100))$root
}
## Do.call = function(list,f) do.call(f,list)
do(0,PrMRIpos_psa[1,1],PrMRIpos_s3m1.5[1,1],rpf1.5[1,1])
do2(0,PrMRIpos_psa[1,1],PrMRIpos_s3m1.5[1,1],rpf1.5[1,1])
threshold1.5GG0 = mapply(function(a,b,c) do(0,a,b,c), PrMRIpos_psa[,1],PrMRIpos_s3m1.5[,1], rpf1.5[,1])
threshold1.5GG1 = mapply(function(a,b,c) do(1,a,b,c), PrMRIpos_psa[,2],PrMRIpos_s3m1.5[,2], rpf1.5[,2])
threshold1.5GG2 = mapply(function(a,b,c) do(2,a,b,c), PrMRIpos_psa[,3],PrMRIpos_s3m1.5[,3], rpf1.5[,3])
thresholdGG0 = mapply(function(a,b,c) do(0,a,b,c), PrMRIpos_psa[,1],PrMRIpos_s3m[,1], rpf[,1])
thresholdGG1 = mapply(function(a,b,c) do(1,a,b,c), PrMRIpos_psa[,2],PrMRIpos_s3m[,2], rpf[,2])
thresholdGG2 = mapply(function(a,b,c) do(2,a,b,c), PrMRIpos_psa[,3],PrMRIpos_s3m[,3], rpf[,3])
threshold = cbind(thresholdGG0,thresholdGG1,thresholdGG2)
threshold1.5 = cbind(threshold1.5GG0,threshold1.5GG1,threshold1.5GG2)
## save(threshold, threshold1.5, PrMRIpos_psa, PrMRIpos_s3m, PrMRIpos_s3m1.5, file="~/Documents/clients/shuang/Study3/Data/sensitivity_inputs.RData")
apply(threshold,2,mean)
apply(threshold1.5,2,mean)



## To split by S3M+MRI+ vs MRI-/S3M>=25; Pr(MRI+S3M+ | S3M+MRI+ vs MRI-/S3M>=25) =
PrS3MposIfBx = 1/adjBoth
## 
PrBxIfS3M1.5pos = PrMRIpos_s3m1.5(0:2)*adjBoth



## Pr(Survival to age 55 years)
library(prostata)
fit <- callFhcrc(1e3,"noScreening",pop=1960,flatPop=TRUE,mc.cores=2)
with(fit, binom.test(x=sum(subset(summary$prev,age==55)$count), n)) # 951/1e2
fit <- callFhcrc(1e4,"noScreening",pop=1960,flatPop=TRUE,mc.cores=2)
with(fit, binom.test(x=sum(subset(summary$prev,age==55)$count), n)) # 9525/1e4
fit <- callFhcrc(1e5,"noScreening",pop=1960,flatPop=TRUE,mc.cores=2)
with(fit, binom.test(x=sum(subset(summary$prev,age==55)$count), n)) # 95189/1e5
fit <- callFhcrc(1e6,"noScreening",pop=1960,flatPop=TRUE,mc.cores=2)
with(fit, binom.test(x=sum(subset(summary$prev,age==55)$count), n)) # 951355/1e6
fit <- callFhcrc(1e7,"noScreening",pop=1960,flatPop=TRUE,mc.cores=2)
with(fit, binom.test(x=sum(subset(summary$prev,age==55)$count), n)) # 9512712/1e7

## Assess Monte Carlo errors - Shuang
library(prostata)
n.sim <- 1e6
mc.cores <- 2
ShuangParametersPorpusU <- prostata:::ShuangParameters()
ShuangParameters_BaseTBx <- prostata:::ShuangParameters()
ShuangParameters_BaseTBx$pTBxG0ifG1_MRIpos=0.2474775    # Pr(TBx gives ISUP 0 | ISUP 1, MRI+) 2020-03-23
ShuangParameters_BaseTBx$pTBxG0ifG2_MRIpos=0.06570613   # Pr(TBx gives ISUP 0 | ISUP 2, MRI+) 2020-03-23
ShuangParameters_BaseTBx$cost_parameters["Combined biopsy"] =
    2990 +                  # Biopsy cost (TBx), still called Combined biopsy
    4238.25                 # Pathology of biopsy
ShuangParameters_BaseTBx$cost_parameters["Active surveillance - yearly - with MRI"] =
    1460+                   # Urology visit and nurse visit
    355.82*3+               # PSA sampling
    57.4*3+                 # PSA analysis
    3500*0.33+              # MRI cost
    2990*0.33+              # Biopsy cost (TBx)
    4238.25*0.33            # Pathology of biopsy
## fit3: MRI+TBx
fit3 <- callFhcrc(n.sim, screen="regular_screen",
                       flatPop=TRUE,pop=1995-55,
                       parms=modifyList(ShuangParameters_BaseTBx, #Using different parameters from fit1, fit2 and fit3
                                        list(MRI_screen=TRUE,
                                             start_screening=55,
                                             screening_interval=4,
                                             formal_compliance=1,
                                             indiv_reports=TRUE,
                                             startReportAge=55)),
                       mc.cores=mc.cores)
# fit4: MRI+TBx/SBx
fit4 <- callFhcrc(n.sim, screen="regular_screen",
                       flatPop=TRUE,pop=1995-55,
                       parms=modifyList(ShuangParametersPorpusU, # check: is this the correct set?
                                        list(MRI_screen=TRUE,
                                             start_screening=55,
                                             screening_interval=4,
                                             formal_compliance=1,
                                             indiv_reports=TRUE,
                                             startReportAge=55)),
                       mc.cores=mc.cores)
summary(fit3,from=55)
fit3$mean_utilities
fit3$mean_costs
summary(fit4,from=55)
fit4$mean_utilities
fit4$mean_costs
mean(fit3$indiv_utilities-fit4$indiv_utilities,na.rm=TRUE)
sd(fit3$indiv_utilities-fit4$indiv_utilities,na.rm=TRUE)/sqrt(fit3$mean_utilities$n)
mean(fit3$indiv_costs-fit4$indiv_costs,na.rm=TRUE)
sd(fit3$indiv_costs-fit4$indiv_costs,na.rm=TRUE)/sqrt(fit3$mean_costs$n)

## Monte Carlo uncertainty in the ICER
## NB: do *not* use do.boot=TRUE with n.sim=1e6 and R=1000 locally, as it segfaults.
mc.uncertainty <- function(title,fit1_costs,fit1_utilities,fit2_costs,fit2_utilities,R=1000,do.boot=FALSE) {
    out <- local({
        index <- !is.na(fit2_costs)
        X <- (fit2_costs-fit1_costs)[index]
        Y <- (fit2_utilities-fit1_utilities)[index]
        muX <- mean(X)
        muY <- mean(Y)
        n <- length(X)
        list(icer=muX/muY,
             d_costs=muX,
             d_utilities=muY,
             se_d_costs=sd(X)/sqrt(n-1),
             se_d_utilities=sd(Y)/sqrt(n-1),
             coefVar_costs=sd(X)/muX/sqrt(n-1),
             coefVar_utilities=sd(Y)/muY/sqrt(n-1),
             log_costs=log(muX),
             log_utilities=log(muY),
             se_log_costs=sd(X)/muX/sqrt(n-1),
             se_log_utilities=sd(Y)/muY/sqrt(n-1),
             se_log_ratio=sqrt((sd(X)/muX/sqrt(n-1))^2+(sd(Y)/muY/sqrt(n-1))^2))
    })
    out$icer.lower = with(out, icer*exp(-1.96*se_log_ratio))
    out$icer.upper = with(out, icer*exp(1.96*se_log_ratio))
    if (do.boot) {
        stopifnot(require(boot))
        index <- !is.na(fit2_costs)
        d <- data.frame(x=(fit2_costs-fit1_costs),
                        y=fit2_utilities-fit1_utilities)[index,]
        boot1 <- boot(d,
                      function(d, w) log(mean(d$x[w])/mean(d$y[w])),
                      R=R)
        out <- c(list(boot=boot1),out)
    }
    cat(title)
    out
}
mc.uncertainty("", fit3$indiv_costs,fit3$indiv_utilities,fit4$indiv_costs,fit4$indiv_utilities)

## Plot of costs over ages
## checks
ShuangParameters_BaseTBx$cost_parameters-
    ShuangParametersPorpusU$cost_parameters # only two costs have changed
all(names(ShuangParameters_BaseTBx$cost_parameters) ==
    names(ShuangParametersPorpusU$cost_parameters))
stacked <- function(x, y, col=1:ncol(y), ...) {
    cumy <- t(apply(cbind(0,y),1,cumsum))
    matplot(x,cumy,type="n",...)
    for (i in 1:ncol(y))
        polygon(c(x,rev(x)), c(cumy[,i+1],rev(cumy[,i])), border=col[i], col=col[i])
}
xtabs(costs~item,fit3$societal.costs, subset=age>=50)
xtabs(costs~item,fit4$societal.costs, subset=age>=50)
costs <- fit3$societal.costs
## xtabs(costs~item,costs,subset=age>=50)
dput(levels(costs$item))
## combine levels
levels(costs$item) <-
    c("Active surveillance", "ADT+chemo", "Diagnostic",
      "Diagnostic", "Diagnostic", "Productivity losses", "Diagnostic", "PSA testing",
      "Palliative therapy", "Post-Tx follow-up",
      "Post-Tx follow-up", "Post-Tx follow-up",
      "Prostatectomy", "Radiation therapy", "Terminal illness")
## xtabs(costs~item,costs,subset=age>=50)
## xtabs(costs~item,costs,subset=age>=50)
## reorder levels
costs$item <-
    factor(costs$item,
           levels=
               c("PSA testing","Diagnostic",
                 "Active surveillance", "Prostatectomy", "Radiation therapy", "ADT+chemo",
                 "Productivity losses", "Post-Tx follow-up",
                 "Palliative therapy",
                 "Terminal illness"))
## xtabs(costs~item,costs,subset=age>=50)
##
library(Cairo)
CairoPDF("~/Documents/clients/shuang/Study2/costsByAge.pdf")
library(RColorBrewer)
tab <- xtabs(costs~age+item,costs,subset=age>=54)
cols <- brewer.pal(ncol(tab), "Spectral")
stacked(as.numeric(rownames(tab)), tab/fit3$n, xlab="Age (years)", ylab="Age-specific annual cost (€) per man", col=cols)
legend("topright",legend=rev(levels(costs$item)),fill=rev(cols),bty="n")
dev.off()

## Changing the startReportAge
library(prostata)
## undebug(callFhcrc)
sim1 <- callFhcrc(1e3,screen="regular_screen",cohort=1960,
          parms=list(indiv_reports=TRUE,
                     screening_interval=2, start_screening=55,
                     startReportAge=55, full_report=1),mc.cores=2)
summary(sim1,from=55) # post hoc adjustment
sim1$mean_utilities # in-simulations calculations - including StdErrs
sim1$mean_costs # in-simulations calculations - including StdErrs

## Timing for full reporting vs brief reporting
library(prostata)
## callFhcrc(1e6,screen="regular_screen",cohort=1960,
##           parms=list(indiv_reports=TRUE,
##                      screening_interval=2, start_screening=55),mc.cores=2) # 350 s
## callFhcrc(1e5,screen="regular_screen",cohort=1960,
##           parms=list(indiv_reports=TRUE,
##                      screening_interval=2, start_screening=55),mc.cores=2) # 35 s
sim3 <- callFhcrc(1e4,screen="regular_screen",cohort=1960,
                  parms=list(indiv_reports=TRUE,
                             screening_interval=2, start_screening=55,
                             full_report=1), mc.cores=2) # 3.5 s
sim4 <- callFhcrc(1e4,screen="regular_screen",cohort=1960,
                  parms=list(indiv_reports=TRUE,
                             screening_interval=2, start_screening=55,
                             full_report=0), mc.cores=2) # 2.4 s
sim5 <- callFhcrc(1e4,screen="regular_screen",cohort=1960,
                  parms=list(indiv_reports=TRUE,
                             screening_interval=2, start_screening=55,
                             full_report=-1), mc.cores=2) # 1.8 s
summary(sim3)
summary(sim4)
## summary(sim5) # not currently available with full_report=-1
sim5$mean_costs
sim5$mean_utilities

## Assess Monte Carlo errors
library(prostata)
sim1 = callFhcrc(1e4,screen="regular_screen",cohort=1960,parms=list(indiv_reports=TRUE,
                                                                    screening_interval=2, start_screening=55),mc.cores=2)
sim2 = callFhcrc(1e4,screen="regular_screen",cohort=1960,parms=list(indiv_reports=TRUE,
                                                                    screening_interval=4, start_screening=55),mc.cores=2)
##
mean(sim1$indiv_costs-sim2$indiv_costs)
sd(sim1$indiv_costs-sim2$indiv_costs)/sqrt(sim1$n)
mean(sim1$indiv_utilities-sim2$indiv_utilities)
sd(sim1$indiv_utilities-sim2$indiv_utilities)/sqrt(sim1$n)
##
mean(sim1$indiv_utilities)
sd(sim1$indiv_utilities)/sqrt(sim1$n)
sim1$mean_utilities


## CAP reconstruction
library(prostata)
library(parallel)
parms = modifyList(prostata:::EdnaParameters(), list(includeEventHistories=FALSE))
## run simulations
cl <- makeCluster(detectCores(),type="PSOCK")
sim1 = callFhcrc(n=sum(cap_control$pop),screen="cap_control",pop=cap_control, cl=cl,
                 parms=parms,
                 tables=list(rescreening=uk_rescreening),
                 nLifeHistories=sum(cap_control$pop))
sim2 = callFhcrc(n=sum(cap_study$pop),screen="cap_study",pop=cap_study, cl=cl,
                 parms=parms,
                 tables=list(rescreening=uk_rescreening),
                 nLifeHistories=sum(cap_study$pop))
stopCluster(cl)
cap_study_by_year <- data.frame(year=seq(2002,2009),
                                p=c(9.89,13.79,17.23,17.756,15.944,21.348,4.06,0))
cap_control_by_year <- data.frame(year=seq(2002,2009),
                                  p=c(7.18,5.21,9.15,29.56,25.24,20.88,2.47,0.31))
## Incidence analysis
cap_incidence <- function(sim, year, p) {
    set.seed(12345) # NB: set the seed for re-producibility
    df <- transform(sim$parameters,
                    yearFup=2016.25,
                    yearEntry=sample(year,
                                     size=nrow(sim$parameters),
                                     replace=TRUE,
                                     prob=p)+
                        runif(nrow(sim$parameters)))
    df <- transform(df, cohort=yearEntry-ageEntry) # actual birth cohort (double)
    df <- subset(df, age_d>ageEntry & (age_pca==-1 | age_pca>ageEntry))
    df <- transform(df,
                    cap_pca=(age_pca!=-1 & age_pca<age_d & cohort+age_pca<yearFup))
    df <- transform(df,
                    pt=ifelse(cap_pca, age_pca, pmin(age_d, yearFup-cohort))-ageEntry)
    df
}
controls <- cap_incidence(sim1, cap_control_by_year$year, cap_control_by_year$p)
study <- cap_incidence(sim2, cap_study_by_year$year, cap_study_by_year$p)
## Plot of cumulative hazards curves for incidence // Figure 2 Panel B from Martin et al (2018)
library(survival)
plot(survfit(Surv(pt,cap_pca)~1, study), fun="cumhaz")
lines(survfit(Surv(pt,cap_pca)~1, controls),col=2,fun="cumhaz")
legend("topleft", c("Study arm","Control arm"), lty=1,col=1:2,bty="n")
## Table: Counts by age at study entry and Gleason score, simulated CAP control arm
xtabs(~I(floor(ageEntry/5)*5)+future_ext_grade, controls, subset=cap_pca)
## Cox regression for incidence
joint <- rbind(transform(controls,study=0),
               transform(study,study=1))
summary(coxph(Surv(pt,cap_pca)~study,joint)) ## HR=1.17 (1.14, 1.21) cf. 1.19 (1.14,1.25) from Martin et al
##
## Mortality analysis
cap_death <- function(sim, year, p) {
    set.seed(12345) # NB: set the seed for re-producibility
    df <- transform(sim$parameters,
                    yearFup=2016.25,
                    yearEntry=sample(year,
                                     size=nrow(sim$parameters),
                                     replace=TRUE,
                                     prob=p)+
                        runif(nrow(sim$parameters)))
    df <- transform(df, cohort=yearEntry-ageEntry) # actual birth cohort (double)
    df <- subset(df, age_d>ageEntry & (age_pca==-1 | age_pca>ageEntry))
    df <- transform(df,
                    cap_dth=(pca_death & cohort+age_d<yearFup),
                    pt=pmin(age_d, yearFup-cohort)-ageEntry)
    df
}
controls <- cap_death(sim1, cap_control_by_year$year, cap_control_by_year$p)
study <- cap_death(sim2, cap_study_by_year$year, cap_study_by_year$p)
## Plot of cumulative hazards curves for incidence // Figure 2 Panel A from Martin et al (2018)
library(survival)
plot(survfit(Surv(pt,cap_dth)~1, study), fun="cumhaz")
lines(survfit(Surv(pt,cap_dth)~1, controls),col=2,fun="cumhaz")
legend("topleft", c("Study arm","Control arm"), lty=1,col=1:2,bty="n")
## Table: Counts by age at study entry and Gleason score, simulated CAP control arm
xtabs(~I(floor(ageEntry/5)*5)+future_ext_grade, controls, subset=cap_dth)
## Cox regression for incidence
joint <- rbind(transform(controls,study=0),
               transform(study,study=1))
summary(coxph(Surv(pt,cap_dth)~study,joint)) ## HR=0.91 (0.80,1.03) cf. 0.96 (0.85, 1.08) from Martin et al
##
## p-values
local({
    a <- c(0.91, 0.80,1.03)
    b <- c(0.96, 0.85, 1.08)
    statistic <- (a[1]-b[1])/sqrt(((a[3]-a[2])/2/1.96)^2+((b[3]-b[2])/2/1.96)^2)
    (1-pnorm(abs(statistic)))*2
})
local({
    a <- c(1.17, 1.14, 1.21)
    b <- c(1.19, 1.14,1.25)
    statistic <- (a[1]-b[1])/sqrt(((a[3]-a[2])/2/1.96)^2+((b[3]-b[2])/2/1.96)^2)
    (1-pnorm(abs(statistic)))*2
})


## UK re-calibration for rescreening
## Assume: cure and shape related to Sweden; given p1, solve for scale
library(prostata) # rescreening data-frame
subset(transform(rescreening, ever=1-cure, p1 = (1-cure)*pweibull(1,shape,scale)),
       age5>=50 & age5<70)
cprd <- rbind(data.frame(age=50,psa=c(0,3,4,6,10,20),p1=c(.09,.27,.61,.67,.63,.49)),
              data.frame(age=55,psa=c(0,3,4,6,10,20),p1=c(.11,.20,.55,.62,.58,.47)),
              data.frame(age=60,psa=c(0,3,4,6,10,20),p1=c(.13,.18,.50,.58,.57,.42)),
              data.frame(age=65,psa=c(0,3,4,6,10,20),p1=c(.15,.20,.40,.57,.53,.38)))
solver <- function(age.,psa.) {
    row <- subset(rescreening,age5==floor(age./5)*5 &
                              total_cat==ifelse(psa.<3,1,ifelse(psa.<10,3,10)))
    shape <- row$shape[1]
    cure <- row$cure[1]
    p1 <- subset(cprd, age==age. & psa==psa.)$p1[1]
    objective <- function(scale) (1-cure)*pweibull(1,shape,scale) - p1
    data.frame(age5=age., total_cat=psa., shape=shape, cure=cure,
               scale=uniroot(objective, lower=0.01, upper=10)$root)
}
plotter <- function(row,t=seq(0,30,length=301))
    plot(t, 0.9*pweibull(t,row$shape,row$scale), type="l")
plotter(solver(50,0))
uk_rescreening1 <-
    do.call(rbind,
            lapply(c(50,55,60,65),
                   function(age) do.call(rbind,lapply(c(0,3,4,6,10,20),
                                                      function(psa) solver(age,psa)))))
uk_rescreening2 <- rbind(transform(subset(rescreening,age5>=80 & total_cat==1), total_cat=0),
                         subset(rescreening,age5>=80 & total_cat==3),
                         transform(subset(rescreening,age5>=80 & total_cat==3), total_cat=4),
                         transform(subset(rescreening,age5>=80 & total_cat==3), total_cat=6),
                         subset(rescreening,age5>=80 & total_cat==10),
                         transform(subset(rescreening,age5>=80 & total_cat==10), total_cat=20))
uk_rescreening <- rbind(transform(subset(uk_rescreening1,age5==50),age5=30),
                         uk_rescreening1,
                         uk_rescreening2)
uk_rescreening <- uk_rescreening[with(uk_rescreening, order(age5,total_cat)),]
dput(uk_rescreening)
##
library(prostata)
sim1 = callFhcrc(1e4,screen="cap_control",pop=cap, mc.cores=2)
sim2 = callFhcrc(1e4,screen="cap_control",pop=cap, mc.cores=2, tables=list(rescreening=uk_rescreening))
plot(sim1,type="testing.rate")
lines(sim2,type="testing.rate",lty=2)
plotter <- function(row,t=seq(0,30,length=301), ...)
    plot(t, 0.9*pweibull(t,row$shape,row$scale), type="l", ...)
liner <- function(row,t=seq(0,30,length=301), ...)
    lines(t, 0.9*pweibull(t,row$shape,row$scale), ...)
plotter(subset(rescreening, age5==30 & total_cat==10))
liner(subset(uk_rescreening, age5==30 & total_cat==10),lty=2, col="blue")


## test for CAP
library(prostata)
control = callFhcrc(1e4,screen="cap_control",pop=cap_control, mc.cores=2)
study = callFhcrc(1e4,screen="cap_study",pop=cap_study, mc.cores=2)
ICER(study, control)


## Bug with predictions
library(prostata)
library(dplyr)
library(ggplot2)
sim1 <- callFhcrc(1e4, mc.cores=2)
pred1 <- predict(sim1, group = c("age"), type = "incidence.rate")
pred2 <- predict(sim1, group = c("age", "grade"), type = "incidence.rate")
pred3 <- left_join(pred1[-(3:4)], pred2[-c(3,5)], by=c("age","scenario")) %>%
    filter(grade!='Healthy') %>%
    mutate(n=ifelse(is.na(n),0,n),
           rate=n/pt,
           grade=droplevels(grade,'Healthy'))
ggplot(pred3, aes(x=age, y=rate)) + facet_grid(~grade) + geom_line()
## checks
sum(sim1$summary$pt$pt)/1e4
table(sim1$summary$events$event)
sum(subset(sim1$summary$events,event %in% c("toScreenDiagnosis","toClinicalDiagnosis"))$n)
sum(pred1$pt)
sum(pred1$n, na.rm=TRUE)
sum(pred2$pt)
sum(pred2$n, na.rm=TRUE)





## Using parLapply (for Windows)
library(prostata)
library(parallel)
gc = 0.001
fit0 <- callFhcrc(1e4,pop=1960,parms=list(gc=gc))
fit1 <- callFhcrc(1e4,mc.cores=2,pop=1960,parms=list(gc=gc))
cl <- makeCluster(detectCores(),type="PSOCK")
fit2 <- callFhcrc(1e4,cl=cl,pop=1960, parms=list(gc=gc))
stopCluster(cl)
summary(fit0)
summary(fit1)
summary(fit2)

## Testing for using background utilities
library(prostata)
summary(test1 <- callFhcrc(1e3))
## # Old output:
## Screening scenario:            noScreening
## Life expectancy:                79.735905
## Discounted QALE:                28.019726
## Discounted health sector costs: 249.729055
## Discounted societal costs:      259.084014
## Discounted rate (effect.):      0.030000
## Discounted rate (costs):        0.030000

library(prostata)
new.table <- list(background_utilities=data.frame(lower=c(0,5,15,30,45,60,70,80),
                        upper=c(5,15,30,45,60,70,80,1.0e99),
                        utility=c(0.970, 0.983, 0.957, 0.941, 0.903, 0.826, 0.731, 0.642)))
## new.table$background_utilities <- transform(new.table$background_utilities,
##                                             lower=pmax(0,lower-1e-7),
##                                             upper=upper-1e-7)
## bg <- transform(prostata:::fhcrcData$background_utilities,
##                 lower=pmax(0,lower-1e-7),
##                 upper=upper-1e-7)
## summary(callFhcrc(1e4,screen="regular_screen"))
summary(test1 <- callFhcrc(1e4,screen="regular_screen",nLifeHistories=1e5))
summary(test2 <- callFhcrc(1e4,screen="regular_screen",tables=new.table,nLifeHistories=1e5))
## summary(test3 <- callFhcrc(1e4,screen="regular_screen",tables=list(background_utilities=bg),
##                            nLifeHistories=1e5))
merge(as.data.frame(table(test1$lifeHistories$event)),
      as.data.frame(table(test2$lifeHistories$event)), by="Var1", all=TRUE) # toRT, toCancerDeath, toYearlyActiveSurveillance, toYearlyPostTxFollowUp
apply(test1$parameters != test2$parameters, 2, sum) # age_d and pca_death => survival??

library(sqldf)
lh1 <- test1$lifeHistories
lh2 <- test2$lifeHistories
diff1 <- sqldf("select distinct id, event from lh1 except select id, event from lh2")
diff2 <- sqldf("select distinct id, event from lh2 except select id, event from lh1")
sqldf("select event, count(*) from diff1 group by event")
sqldf("select event, count(*) from diff2 group by event")
sqldf("select id from diff1 where event=\"toScreenInitiatedBiopsy\"")

options(width=200)
subset(test1$lifeHistories,id==138)
subset(test2$lifeHistories,id==138)


which(test1$parameters != test2$parameters)
which(test1$parameters$age_pca != test2$parameters$age_pca)

rbind(subset(test1$parameters, id==752),
      subset(test2$parameters, id==752))
options(width=200)
subset(test1$lifeHistories,id==752)
subset(test2$lifeHistories,id==752)

rbind(subset(test1$parameters, id==80),
      subset(test2$parameters, id==80))
options(width=200)
subset(test1$lifeHistories,id==80)
subset(test2$lifeHistories,id==80)



## Issue with small differences
library(prostata)
library(dplyr)
test1 = callFhcrc(1e4, nLifeHistories=1e4, mc.cores=2)
test2 = callFhcrc(1e4, screen="regular_screen",nLifeHistories=1e4, mc.cores=2,
                  parms=list(screening_interval=4))
summary(test1)
summary(test1,from=55)
ICER(test1, test2)
unlist(ICER(test1,test2))[2:3]*1000
summary(test1)
## approximate QALYs, costs and LE
QALE <- function(sim, discountRate=0)
    with(sim, summary$ut %>% mutate(ut=ut*((1+simulation.parameters$discountRate.effectiveness)/(1+discountRate))^(age+0.5)) %>% summarise(ut=sum(ut)/n))
costs <- function(sim, discountRate=0)
    with(sim, societal.costs %>% mutate(costs=costs*((1+simulation.parameters$discountRate.costs)/(1+discountRate))^(age+0.5)) %>% summarise(costs=sum(costs)/n))
LE <- function(sim, discountRate = 0)
    with(sim, summary$pt %>% mutate(pt=pt/(1+discountRate)^(age+0.5)) %>% summarise(pt=sum(pt)/n))
## differences in measures with and without discounting
(QALE(test1)-QALE(test2))*1000
(QALE(test1,0.03)-QALE(test2,0.03))*1000
(costs(test1)-costs(test2))*1000
(costs(test1,0.03)-costs(test2,0.03))*1000
(LE(test1)-LE(test2))*1000
(LE(test1,0.03)-LE(test2,0.03))*1000

QALE <- function(sim, discountRate=0, centre=0, minage=0)
    with(sim, filter(summary$ut, age>=minage) %>% mutate(ut=ut*(1+simulation.parameters$discountRate.effectiveness)^(age+0.5)/(1+discountRate)^((age-centre)+0.5)) %>% summarise(ut=sum(ut)/n))
(QALE(test1)-QALE(test2))*1000
(QALE(test1,0.03)-QALE(test2,0.03))*1000
(QALE(test1,0.03,55)-QALE(test2,0.03,55))*1000
(QALE(test1,minage=50)-QALE(test2,minage=50))*1000
(QALE(test1,0.03)-QALE(test2,0.03))*1000
(QALE(test1,0.03,55)-QALE(test2,0.03,55))*1000



## ERSPC replication
## Simple: compare no screening with four-yearly screening 60--69 years, with follow-up for fourteen years
library(prostata)
eform <- 
function(object, parm, level = 0.95, method=c("Profile","Wald"), name="exp(beta)") {
  method <- match.arg(method)
  if (missing(parm))
    parm <- TRUE
  estfun <- switch(method, Profile = MASS:::confint.glm, Wald = stats::confint.default)
  val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
  colnames(val) <- c(name,colnames(val)[-1])
  val[parm, ]
}
parms <- list(start_screening = 60,
              stop_screening = 70,
              screening_interval = 4,
              RR_T3plus = 2,
              formal_compliance=1,
              c_txlt_interaction = 0.75^(1/10))
parms2Int <- parms4Int <- parms2 <- parms4 <- parms
parms4Int$c_txlt_interaction <- parms2Int$
parms2Int$screening_interval <- parms2$screening_interval <- 2

noscreen <- callFhcrc(1e6,screen="noScreen",mc.cores=2,pop=1990-60,parms=parms)
screen4 <- callFhcrc(1e6,screen="regular_screen",mc.cores=2,pop=1990-60,parms=parms4)
screen2 <- callFhcrc(1e6,screen="regular_screen",mc.cores=2,pop=1990-60,parms=parms2)
screen4Int <- callFhcrc(1e6,screen="regular_screen",mc.cores=2,pop=1990-60,parms=parms4Int)
screen2Int <- callFhcrc(1e6,screen="regular_screen",mc.cores=2,pop=1990-60,parms=parms2Int)

summary(noscreen)
summary(screen4)
summary(screen4Int)
summary(screen2Int)
plot(screen, type="testing.rate")
plot(noscreen, type="pc.mortality.rate")
lines(screen4Int, type="pc.mortality.rate", col="red")
lines(screen2Int, type="pc.mortality.rate", col="green")

compare <- function(ref,exposed) {
    rates <- rbind(transform(subset(predict(ref, type="pc.mortality.rate"), age %in% 60:75), exposed=0),
                   transform(subset(predict(exposed, type="pc.mortality.rate"), age %in% 60:75), exposed=1))
    fit <- glm(n ~ exposed + offset(log(pt)), data=rates, family=poisson)
    eform(fit)
}
compare(noscreen, screen4)
compare(noscreen, screen2)
compare(noscreen, screen4Int)
compare(noscreen, screen2Int)


##
0.834 (4 interval)
0.782 (4 interval + int@0.95)
0.761 (2 interval + int)
0.737 (0.5 interval + int)

## In summary, varying RR_T3plus between 1 and 2 leads to approximately a 2% reduction in the prostate cancer mortality rate ratio; increasing RR_T3plus to 3 led to another 1% reduction in the rate ratio.



## re-check survival
library(prostata)
today <- callFhcrc(1e5,screen="screenUptake",mc.cores=2,pop=1990-60,
                   includeDiagnoses=TRUE, nLifeHistories=1e5)
dim(today$diagnoses)
subset(today$diagnoses, age_at_death<0)
## TODO competing risks
library(survival)
tmp <- subset(transform(today$diagnoses, age10 = floor(age/10)*10), age10>=50 & age10<90)
##
## Metastatic and loco-regional
par(mfrow=c(1,2))
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Metastatic"), col=1:4, main="Metastatic")
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised"), col=1:4, main="Loco-regional")
legend("topright",legend=levels(factor(tmp$age10)),col=1:4,lty=1)
##
## Loco-regional by PSA
par(mfrow=c(1,2))
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised" & psa<10), col=1:4, main="Loco-regional, PSA<10")
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised" & psa>=10), col=1:4, main="Loco-regional, PSA>=10")
legend("topright",legend=levels(factor(tmp$age10)),col=1:4,lty=1)
##
##
## Loco-regional by Gleason
par(mfrow=c(2,2))
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised" & ext_grade=="Gleason_le_6"), col=1:4, main="Loco-regional, Gleason 6")
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised" & ext_grade=="Gleason_7"), col=1:4, main="Loco-regional, Gleason 7")
plot(fit <- survfit(Surv(age_at_death - age, cancer_death) ~ age10, data=tmp, subset=state=="Localised" & ext_grade=="Gleason_ge_8"), col=1:4, main="Loco-regional, Gleason 8+")
legend("topright",legend=levels(factor(tmp$age10)),col=1:4,lty=1)

## Predicted risk to 10 and 15 years
tmp2 <- transform(droplevels(subset(tmp, state == "Localised")), psa10=(psa>=10), t3plus=ext_state=="T3plus", age10=factor(age10),time=age_at_death - age)
fit <- coxph(Surv(time, cancer_death) ~ age10+ext_grade+psa10+t3plus, data=tmp2) # assumes multiplicitive
newd <- unique(subset(tmp2,select=c(age10,ext_grade,t3plus,psa10)))
newd <- newd[with(newd, order(age10,ext_grade,t3plus,psa10)),]
pred <- data.frame(newd,
                   surv10=as.vector(summary(survfit(fit,newdata=newd), time=10)$surv),
                   surv15=as.vector(summary(survfit(fit,newdata=newd), time=15)$surv))
xtabs(surv10 ~ age10 + ext_grade+t3plus+psa10, data=pred)



## xtabs(Survival~ AgeLow+Grade,data=fhcrcData$survival_local, subset=Time==10)

## 9013 bug
require(microsimulation)
debug(callFhcrc)
out=callFhcrc(1e5,nLifeHistories=1e5,screen="screenUptake",mc.cores=4)

tmp <- out[[1]]$parameters
sapply(tmp,length)
with(c(tmp,list(i=9013+1)), data.frame(id=id[i],age_d=age_d[i],aoc=aoc[i],pca_death=pca_death[i]))

subset(lifeHistories, id==0)

with(lapply(tmp,function(var) var[9013+1]), tmc+35)
names(lifeHistories) <- c("id", "ext_state", "ext_grade", "dx", "event", "begin", "end", "year", "psa", "utility")
options(width=150)
subset(lifeHistories, id==9013)
`names<-`(0:(length(eventT)-1),eventT)


## unit tests
refresh
require(microsimulation)
temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp2=callFhcrc(1000)
stopifnot(abs(with(temp2,sum(summary$pt$pt)/n)-79.88935)<1e-3)
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)

temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)

refresh
require(microsimulation)
require(sqldf)
temp2=callFhcrc(1e7,screen="screenUptake",mc.cores=2)

events <- temp2$summary$events
pt <- temp2$summary$pt
m <- dbDriver("SQLite")
connection <- dbConnect(m, dbname = ":memory:")
init_extensions(connection)
dbWriteTable(connection, '`events`', events, row.names = FALSE)
dbWriteTable(connection, '`pt`', pt, row.names = FALSE)

## cancer mortality rates
t1 <- dbGetQuery(connection, "select *, n/pt as rate from (select year, min(85,floor(age/5)*5) as age5, sum(n*1.0) as n from events where event in ('toCancerDeath') and year>=1995 and year<=2010 group by year, age5) t1 natural join (select year, min(85,floor(age/5)*5) as age5, sum(pt) as pt from pt where year>=1995 and year<=2010 group by year, age5) as t2 order by t1.age5, t1.year")

## incidence rates
t1 <- dbGetQuery(connection, "select *, n/pt as rate from (select year, min(85,floor(age/5)*5) as age5, sum(n*1.0) as n from events where event in ('toClinicalDiagnosis','toScreenDiagnosis') and year>=1995 and year<=2010 group by year, age5) t1 natural join (select year, min(85,floor(age/5)*5) as age5, sum(pt) as pt from pt where year>=1995 and year<=2010 group by year, age5) as t2 order by t1.age5, t1.year")
require(lattice)
xyplot(rate ~ year | factor(age5), data=t1, type="l")

dbGetQuery(connection, 'select event,sum(n) from events group by event')

temp=callFhcrc(1e5,nLifeHistories=1e5)
life=temp$lifeHistories
deaths <- sqldf("select t1.id, t1.end as dx, t2.end as death, t1.state, t1.ext_grade from life as t1 inner join life as t2 on t1.id=t2.id where t1.event='toClinicalDiagnosis' and t2.event='toCancerDeath'")

dbDisconnect(connection)

with(subset(deaths,state=="Localised"),plot(density(death-dx,from=0)))

## what proportion of clinical diagnoses die due to cancer?
deaths <- sqldf("select t1.id, t1.end as dx, t2.end as death, t2.event, t1.state, t1.ext_grade from life as t1 inner join life as t2 on t1.id=t2.id where t1.event='toClinicalDiagnosis' and t2.event in ('toOtherDeath','toCancerDeath')")
xtabs(~event+pmin(85,floor(death/5)*5),deaths,subset=(state=="Localised"))

refresh
require(microsimulation)
model <- function(..., mc.cores=3) callFhcrc(..., mc.cores=mc.cores)
model0.0 <- model(1e6,discountRate=0)
model0 <- model(1e6)
model2 <- model(1e6,screen="twoYearlyScreen50to70")
model2.0 <- model(1e6,screen="twoYearlyScreen50to70",discountRate=0)
model4 <- model(1e6,screen="fourYearlyScreen50to70")
ICER(model0.0,model2.0)
ICER(model0,model2)
ICER(model0,model4)
ICER(model2,model4)

## ext_state
refresh
require(microsimulation)
report <- function(obj) {
    rowpct <- function(m) t(apply(m,1,function(row) row/sum(row)))
    print(xtabs(~I(floor(age/5)*5)+ext_state, data=obj$diagnoses))
    print(rowpct(xtabs(~I(floor(age/5)*5)+ext_state, data=obj$diagnoses)))
    rowpct(xtabs(~I(floor(age/5)*5)+ext_state, data=obj$diagnoses, subset=ext_state %in% c('T1_T2','T3plus')))
}
test <- callFhcrc(1e5,includeDiagnoses=TRUE,mc.cores=2)
report(test)
test2 <- callFhcrc(1e6,screen="screenUptake",includeDiagnoses=TRUE,mc.cores=2, parms=list(gamma_m_diff=-0.001))
test2a <- test2
test2a$diagnoses <- subset(test2$diagnoses,year>=2000)
(report2 <- report(test2a))
test3 <- callFhcrc(1e5,screen="twoYearlyScreen50to70",includeDiagnoses=TRUE,mc.cores=2)
report(test3)
s0 <- data.frame(age=c(40,45,50,55,60,65,70,75,80,85,90), n_M0=c(16,99,576,1231,2223,3170,1713,870,474,195,48), 
pT1_T2=c(0.9375,0.949495,0.939236,0.918765,0.891138,0.901262,0.841214,0.791954,0.71308,0.651282,0.416667))
plot(as.numeric(rownames(report2)),report2[,2], type="l")
with(s0, lines(age,pT1_T2,lty=2))


#### Calibrate for survival
require(microsimulation)
require(dplyr)
## competing risks - using FhcrcParameters$mu0 and fhcrcData$survival_local
CR <- function(agedx,times=c(10,15),HR=1,grade=0) {
    S0 <- data.frame(age=0:105,mu0=FhcrcParameters$mu0) %>%
        filter(age>=agedx) %>%
            mutate(Time=age-agedx,S0=exp(-cumsum(c(0,mu0[-length(mu0)])))) %>%
                select(Time,mu0,S0)
    S1 <- filter(fhcrcData$survival_local,AgeLow==agedx,Grade==grade) %>%
        mutate(S1=Survival^HR,mu1=-HR*log(c(Survival[-1],NA)/Survival)) %>%
            select(Time,mu1,S1)
    inner_join(S0,S1,by="Time") %>%
        mutate(x0=S0*S1*mu0, x1=S0*S1*mu1, CR0=cumsum(x0), CR1=cumsum(x1)) %>%
            filter(Time %in% (times-1)) %>%
                mutate(Time=times) %>%
                    select(Time,CR1,CR0)
}
CRsolve <- function(data) {
    data$grade <- ifelse(data$grade %in% c(0,6,7),0,1)
    optimize(function(hr)
             CR(unique(data$age),HR=hr,grade=data$grade) %>% inner_join(data, by="Time") %>%
             with(sum(c(CR1-cr1)^2)),
             c(0.1,10))
}
## PSA<10, M0/MX
CRsolve(data.frame(age=55,grade=6,Time=c(10,15),cr1=c(0.009,0.029))) # HR=0.137
## CRsolve(data.frame(age=50,grade=7,Time=c(10,15),cr1=c(NA,NA))) # too few at risk
CRsolve(data.frame(age=65,grade=6,Time=c(10,15),cr1=c(0.014,0.047))) # HR=0.181
CRsolve(data.frame(age=65,grade=7,Time=c(10,15),cr1=c(0.087,0.198))) # HR=0.874
CRsolve(data.frame(age=75,grade=6,Time=c(10,15),cr1=c(0.049,0.119))) # HR=0.480
CRsolve(data.frame(age=75,grade=7,Time=c(10,15),cr1=c(0.128,0.217))) # HR=1.020
##
## PSA>=10, M0/MX
CRsolve(data.frame(age=55,grade=6,Time=c(10,15),cr1=c(0.060,0.178))) # HR=0.902
CRsolve(data.frame(age=55,grade=7,Time=c(10,15),cr1=c(0.369,0.474))) # HR=3.631
CRsolve(data.frame(age=65,grade=6,Time=c(10,15),cr1=c(0.088,0.188))) # HR=0.840
CRsolve(data.frame(age=65,grade=7,Time=c(10,15),cr1=c(0.329,0.436))) # HR=2.644
CRsolve(data.frame(age=75,grade=6,Time=c(10,15),cr1=c(0.140,0.239))) # HR=1.134
CRsolve(data.frame(age=75,grade=7,Time=c(10,15),cr1=c(0.265,0.363))) # HR=2.035
CRsolve(data.frame(age=75,grade=8,Time=c(10,15),cr1=c(0.463,0.524))) # HR=0.797
## NB: the SEER survival data were NOT stratified by PSA value.

refresh
require(rstpm2)
require(foreign)
require(lattice)
d <- read.dta("~/Downloads/prostate-20141010.dta")
d2 <- subset(d, age>=60 & age<70 & diayear>=1980)
## debug(pstpm2)
fit <- pstpm2(Surv(time,pcdeath)~1,data=d2,smooth.formula=~s(time,diayear,k=30),sp=0.001)

xtabs(~diayear+addNA(m),data=d)
recent <- subset(d,diayear>=2004 & !is.na(m))
xtabs(~I(floor(age/5)*5)+m,recent)
require(sqldf)
sqldf("select min(90,floor(age/5)*5) as age5, count(*) as n, avg(m) as p_m from recent group by age5")

grid <- expand.grid(diayear=1980:2009,
                    time=seq(0.1,6000,length=50))
grid$haz <- predict(fit,newdata=grid,type="hazard")
grid$surv <- predict(fit,newdata=grid)
xyplot(haz~time | diayear, data=grid, type="l")
xyplot(surv~time | diayear, data=grid, type="l")

xyplot(surv~time | factor(diayear), data=grid,
       panel=function(x,y,subscripts) {
           panel.xyplot(x,y,type="l")
           d3 <- subset(d2,diayear == grid$diayear[subscripts][1])
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1)
           panel.lines(sfit$time,sfit$lower,col=1,lty=2)
           panel.lines(sfit$time,sfit$upper,col=1,lty=2)
       })

xyplot(pcdeath~time | factor(diayear), data=d,
       subset=age>=50 & age<70,
       panel=function(x,y,subscripts) {
           d3 <- d[subscripts,]
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1,type="S")
           panel.lines(sfit$time,sfit$lower,col=1,lty=2,type="S")
           panel.lines(sfit$time,sfit$upper,col=1,lty=2,type="S")
       })

xyplot(pcdeath~time | factor(diayear), data=d,
       subset=age>=85 & diayear>=1961,
       xlim=c(0,365*15),
       panel=function(x,y,subscripts) {
           d3 <- d[subscripts,]
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1,type="S")
           panel.lines(sfit$time,sfit$lower,col=1,lty=2,type="S")
           panel.lines(sfit$time,sfit$upper,col=1,lty=2,type="S")
       })




## all(c(callFhcrc(5)$lifeHistories == callFhcrc(5)$lifeHistories,
##       callFhcrc(5)$parameters == callFhcrc(5)$parameters,
##       callFhcrc(5)$summary$pt == callFhcrc(5)$summary$pt,
##       callFhcrc(5)$summary$prev == callFhcrc(5)$summary$prev))
## all(callFhcrc(10)$lifeHistories == callFhcrc(10)$lifeHistories) # fails
## all(callFhcrc(1e4)$parameters == callFhcrc(1e4)$parameters) # okay

## extract PSA values and calculate STHLM3 PSA "pseudo-thresholds" for the biomarker panel
refresh
require(microsimulation)
pos <- function(x) ifelse(x>0,x,0)
set.seed(12345)
temp <- callFhcrc(1e6,screen="stockholm3_risk_stratified",includePSArecords=TRUE,mc.cores=2)$psarecord
temp2 <- subset(temp,organised & age>=50 & age<70 & !dx) 
## first organised screen
i <- tapply(1:nrow(temp2),temp2$id,min)
temp3 <- temp2[i,]
cat("No cancer:\n")
with(subset(temp3,state==0 & psa>3), mean(psa<4.4)) # about 42% (from STHLM3)
cat("Loco-regional Cancer:\n")
with(subset(temp3,state>0 & ext_grade==0 & psa>3), mean(psa<3.6)) # about 17% (from STHLM3)
with(subset(temp3,state>0 & ext_grade==0 & psa>3),
     cumsum(table(cut(delta,c(0,5,10,15,Inf)))/length(delta)))
with(subset(temp3,state>0 & ext_grade==0 & psa>3.6), 
     plot(density(delta,from=0)))
with(subset(temp3,state>0 & ext_grade==0 & psa>3), 
    lines(density(delta,from=0),lty=2))

with(subset(temp3,state>0 & ext_grade==0 & psa>3), mean(delta))
with(subset(temp3,state>0 & ext_grade==0 & psa>3 & psa<3.6), mean(delta))
temp3 <- transform(temp3, delta=age-(t0+35))


temp2 <- transform(temp2,
                   advanced=(state>0 & ext_grade==2),
                   cancer=(state>0),
                   logpsa=log(psa),
                   logZ=log(Z),
                   logZstar=beta0+beta1*(age-35)+pos(beta2*(age-35-t0)),
                   grade=ifelse(ext_grade %in% 0:1,1,ext_grade))
temp2 <- within(temp2, {
    ext_grade <- ifelse(cancer, ext_grade, NA)
})
## ## rnormPos is now in the package...
## rnormPos <- function(n,mean=0,sd=1,lbound=0) {
##     if (length(mean)<n) mean <- rep(mean,length=n)
##     if (length(sd)<n) sd <- rep(sd,length=n)
##     x <- rnorm(n,mean,sd)
##     while(any(i <- which(x<lbound)))
##         x[i] <- rnorm(length(i),mean[i],sd[i])
##     x
## }
## onset ho(t) = g0 * t
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
n <- nrow(temp2)
correlatedValue <- function(y,mu,se=NULL,rho) {
    residual <- y - mu 
    if (is.null(se)) se <- sqrt(var(residual))
    u <- rnorm(length(y),0,se) # marginal error
    mu + rho*residual + sqrt(1-rho^2)*u # new value
}

rtpf <- function(marker1, threshold1, marker2, threshold2, disease) {
    n1 <- sum(marker1[disease] > threshold1)
    n2 <- sum(marker2[disease] > threshold2)
    list(n1=n1, n2=n2, rtpf=n1/n2)
}
rfpf <- function(marker1, threshold1, marker2, threshold2, disease) {
    n1 <- sum(marker1[!disease] > threshold1)
    n2 <- sum(marker2[!disease] > threshold2)
    list(n1=n1, n2=n2, rfpf=n1/n2)
}
variance <- tau2 <- 0.0829
optim1 <- optim(c(log(4.5),log(0.01)),
                function(par) {
                    set.seed(12345)
                    threshold <- par[1]
                    scale <- exp(par[2])
                    temp2$bbp <<- with(temp2, logpsa + scale*pos(age-t0) + rnorm(nrow(temp2), 0, sqrt(variance)))
                    with(temp2,
                         (rtpf(logpsa, log(3.0), bbp, threshold, advanced)$rtpf-1)^2 +
                         (rfpf(logpsa, log(3.0), bbp, threshold, advanced)$rfpf-1.25)^2)
                })
set.seed(12345)
threshold <- optim1$par[1]
scale <- exp(optim1$par[2])
temp2$bbp <- with(temp2, logpsa + scale*pos(age-t0) + rnorm(nrow(temp2), 0, sqrt(variance)))
with(temp2, rtpf(logpsa, log(3.0), bbp, threshold, advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, threshold, advanced))

with(temp2, rtpf(logpsa, log(3.0), bbp, log(10), advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, log(10), advanced))
root1 <- uniroot(function(threshold) with(temp2, rtpf(logpsa, log(3.0), bbp, threshold, advanced))$rtpf-1, c(log(2), log(100)))
with(temp2, rtpf(logpsa, log(3.0), bbp, root1$root, advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, root1$root, advanced))


## correlated biomarkers based on the mean
set.seed(12345+5)
biomarker2 <- exp(correlatedValue(log(temp2$psa),
                                  p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2[temp2$grade]*pos(temp2$age-35-temp2$t0),
                                  rho=0.25))
biomarker2 <- exp(log(temp2$psa) + 0.1*pos(temp2$age-35-temp2$t0)+rnorm(n,0,sqrt(p$tau2)))
if (FALSE) {
    plot(temp2$psa,biomarker2,log="xy")
    sqrt(var(log(temp2$psa) - p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2*pos(temp2$age-35-temp2$t0)))
    cor(log(biomarker2),log(temp2$psa))
    plot(density(log(temp2$psa)))
    lines(density(log(biomarker2)),lty=2)
    var(log(biomarker2))
    var(log(temp2$psa))
}
temp3 <- transform(temp2, BBP=biomarker2)
## STHLM3 simulation report
with(list(threshold=5),with(transform(temp3,BBPpos=(psa>=1 & BBP>=threshold),PSApos=(psa>=3)),
     cat(sprintf("
PSA+ & advanced:\t\t%i
PSA+ & cancer:\t\t\t%i
BBP+ & adv:\t\t\t%i
BBP+ & can:\t\t\t%i
BBP+ & PSA+:\t\t\t%i
BBP+ & PSA-:\t\t\t%i
BBP- & PSA+:\t\t\t%i
PSA+ | BBP+:\t\t\t%i
Prop reduction in biospies:\t%5.3f\n",
                 sum(PSApos & advanced),
                 sum(PSApos & cancer),
                 sum(BBPpos & advanced),
                 sum(BBPpos & cancer),
                 sum(BBPpos & PSApos),
                 sum(BBPpos & !PSApos),
                 sum(!BBPpos & PSApos),
                 sum(BBPpos | PSApos),
                 (sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
                       sum(PSApos)
                 ))))

report <- function(psa, BBP, advanced, threshold=3) {
    BBPpos <- (psa>=1 & BBP>=threshold)
    PSApos <- (psa>=3)
    c(PSAposAdvanced=sum(PSApos & advanced),
      BBPposAdvanced=sum(BBPpos & advanced),
      pBiopsy=(sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
      sum(PSApos))
}
report(temp2$psa,biomarker2,temp2$advanced,threshold=1.8)
reports <- sapply(1:200,function(i) {
    biomarker2 <- exp(correlatedValue(log(temp2$psa),
                                      p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2[temp2$grade]*pos(temp2$age-35-temp2$t0),
                                      rho=0.25))
    report(temp2$psa,biomarker2,temp2$advanced, threshold=1.8)
})
reports <- as.data.frame(t(reports))
with(reports, mean( PSAposAdvanced <= BBPposAdvanced))
with(reports, plot(table( PSAposAdvanced - BBPposAdvanced)))
with(reports, plot(density(pBiopsy)))
with(reports, plot(PSAposAdvanced - BBPposAdvanced,pBiopsy))

## Old natural history - ignoring diagnosis
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
## onset ho(t) = g0 * t, Ho(t) = g0/2*t*t = -logU => t=sqrt(-2*log(U)/g0)
set.seed(12345)
n <- 1e6
age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
## grade <- rep(1:2,c(0.9*n,0.1*n)) # this should depend on age of onset
grade <- ifelse(runif(n) < 0.006*(age_o-35), 1, 0)
beta0 <- with(p, rnorm(n,mubeta0,sebeta0))
beta1 <- with(p, rnormPos(n,mubeta1,sebeta1))
beta2 <- with(p, rnormPos(n,mubeta2[grade],sebeta2[grade]))
eps <- with(p, rnorm(n,0,tau2))
lpsa <- pmin(log(20),beta0+beta1*(50-35)+beta2*pmax(0,50-age_o)+eps)
psacut <- function(x) cut(x,c(0,1,3,10,Inf), right=FALSE)
table(psacut(exp(lpsa)))/length(lpsa)
plot(density(exp(lpsa)),xlim=c(0,20)) # density of PSA at age 50 years
i <- 1
for (age in seq(45,85,by=10)) {
    cat(age,"\n")
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    lines(density(exp(lpsa)),col=i)
    print(table(cut(exp(lpsa),psa_cuts))/length(lpsa))
    i <- i+1
}
tab <- sapply(ages <- seq(55,80,by=5), function(age) {
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o))
    tab <- table(psacut(exp(lpsa)))
    tab/sum(tab)
})
colnames(tab) <- ages
tab



## revised natural history?

p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
psaSim <- function(n,age=50,grade=NULL) {
    age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
    if (is.null(grade)) grade <- rep(1:2,c(0.9*n,0.1*n))
    grade <- rep(grade,length=n)
    beta0 <- with(p, rnorm(n,mubeta0,sebeta0))
    beta1 <- with(p, rnormPos(n,mubeta1,sebeta1))
    beta2 <- with(p, rnormPos(n,mubeta2[grade],sebeta2[grade]))
    eps <- with(p, rnorm(n,0,tau2))
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    psa <- exp(lpsa)
    as.data.frame(as.list(environment()))
}
set.seed(12345)
psaSim(11)

refresh
require(microsimulation)
y <- t(replicate(1000,.Call("rbinormPos_test",package="microsimulation")))
cor(y)
plot(y)


refresh
require(microsimulation)
require(sqldf)
require(dplyr)
load("~/Downloads/IHEdata.RData")
makeModel <- function(discountRate=0.03,
                      formal_compliance=1,
                      formal_costs=1,
                      panel=TRUE) {
    function(screen, ..., parms=NULL, n=1e6, mc.cores=3, pop=1960) {
        newparms <- list(formal_compliance=formal_compliance,
                         formal_costs=formal_costs)
        if (!is.null(parms))
            for (name in names(parms))
                newparms[[name]] <- parms[[name]]
        callFhcrc(n, screen=screen, mc.cores=mc.cores, pop=pop, discountRate=discountRate, parms=newparms, panel=panel, ...)
    }
}
modelSet <- function(model) {
    model0 <- model("noScreening")
    model2 <- model("twoYearlyScreen50to70")
    model4 <- model("fourYearlyScreen50to70")
    model50 <- model("screen50")
    model60 <- model("screen60")
    model70 <- model("screen70")
    modelUptake1930 <- model("screenUptake",pop=1930)
    modelUptake1960 <- model("screenUptake")
    modelGoteborg <- model("goteborg")
    modelRiskStratified <- model("risk_stratified")
    modelMixedScreening <- model("mixed_screening")
    cat("NOTE: Processing completed.\n")
    models <- list(model0,model2,model4,model50,model60,model70,
                   modelUptake1930,modelUptake1960,modelGoteborg,
                   modelRiskStratified,modelMixedScreening)
    names(models) <- c("No screening","2-yearly","4-yearly",
                       "50 only","60 only","70 only","Opportunistic 1930",
                       "Opportunistic 1960+",
                       "Göteborg","Risk stratified",
                         "Mixed screening")
    models
}
predict.fhcrc <-
function (obj, type = c("incidence", "cancerdeath"))
{
    type <- match.arg(type)
    event_types <- switch(type, incidence = c("toClinicalDiagnosis", 
        "toScreenDiagnosis"), cancerdeath = "toCancerDeath")
    if (require(dplyr)) {
        pt <- obj$summary$pt %>% group_by(age) %>% summarise(pt = sum(pt))
        events <- obj$summary$events %>% filter(event %in% event_types) %>% 
            group_by(age) %>% summarise(n = sum(n))
        out <- left_join(pt, events, by = "age") %>% mutate(rate = ifelse(is.na(n), 
            0, n/pt))
        return(with(out,data.frame(age=age,rate=rate,pt=pt,n=ifelse(is.na(n),0,n))))
    }
    else error("dplyr is not available for plotting")
}
plot.scenarios <- function(models,
                           costs="delta.costs",
                           effects="delta.QALE",
                           xlim=NULL, ylim=NULL,
                           ylab="Effectiveness (QALY)",
                           suffix="", prefix="",
                           textp=TRUE,
                           pos=rep(4,length(models)),
                           ...) {
    s <- data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=sprintf("%s%s%s",prefix,names(models),suffix),
                    pos=pos)
    costs <- s[[costs]]
    effects <- s[[effects]]
    plot(costs,
         effects,
         xlim=if (is.null(xlim)) c(0,max(costs)*1.3) else xlim,
         ylim=if (is.null(ylim)) c(0,max(effects)*1.1) else ylim,
         xlab="Costs (SEK)",
         ylab=ylab,
         pch=19, cex=1.5,
         ...)
    if (textp) text(costs,effects, labels=s$model, pos=pos)
    lines.frontier(costs,effects,type="c",lwd=2,col="grey")
}
points.scenarios <- function(models,
                             costs="delta.costs",
                             effects="delta.QALE",
                             suffix="",
                             prefix="",
                             textp = TRUE,
                             pos=rep(4,length(models)), ...) {
    s <- data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=sprintf("%s%s%s",prefix,names(models),suffix),
                    pos=pos)
    costs <- s[[costs]]
    effects <- s[[effects]]
    points(costs,
           effects,
           pch=19,cex=1.5,
           ...)
    if (textp) text(costs,effects, labels=s$model, pos=pos)
}
segments.scenarios <- function(modelsA,
                               modelsB,
                               costs="delta.costs",
                               textp=FALSE,
                               pos=rep(4,length(modelsA)),
                               effects="delta.QALE",
                               ...) {
    sA <- data.frame(t(sapply(modelsA,
                              function(obj) unlist(ICER(obj,modelsA[[1]])))))
    sB <- data.frame(t(sapply(modelsB,
                              function(obj) unlist(ICER(obj,modelsB[[1]])))))
    costsA <- sA[[costs]]
    effectsA <- sA[[effects]]
    costsB <- sB[[costs]]
    effectsB <- sB[[effects]]
    segments(costsA,effectsA,
             costsB,effectsB,
             lwd=2,
             ...)
    if (textp)
        text((costsA+costsB)/2,
             (effectsA+effectsB)/2,
             labels=names(modelsA),
             pos=pos)
}
summary.scenarios <- function(models) {
    data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=names(models))
}
post <- function(modelSet) {
    i <- c(1,4,5,6,8:11)
    names(modelSet) <- c("No screening","Göteborg","4-yearly",
                         "50 only","60 only","70 only","Opportunistic 1930",
                         "Opportunistic",
                         "Risk stratified (2+4)","Risk stratified (4+8)",
                         "Mixed screening")
    modelSet[i]
}
if (FALSE) {
    modelSetA <- modelSet(makeModel(discount=0,formal_compliance=0,formal_costs=0,panel=FALSE))
    modelSetB <- modelSet(makeModel(discount=0.03,formal_compliance=0,formal_costs=0,panel=FALSE))
    modelSetC <- modelSet(makeModel(discount=0.03,formal_compliance=1,formal_costs=1,panel=FALSE))
    modelSetD <- modelSet(makeModel(discount=0.03,formal_compliance=1,formal_costs=1,panel=TRUE))
    modelSetBD <- modelSet(makeModel(discount=0.03,formal_compliance=0,formal_costs=0,panel=TRUE))
    save(modelSetA,file="~/work/modelSetA-20150201.RData")
    save(modelSetB,file="~/work/modelSetB-20150201.RData")
    save(modelSetC,file="~/work/modelSetC-20150201.RData")
    save(modelSetD,file="~/work/modelSetD-20150201.RData")
    save(modelSetBD,file="~/work/modelSetBD-20150201.RData")
}
doOnce <- TRUE
if (doOnce) {
    load("~/work/modelSetA-20150201.RData")
    load("~/work/modelSetB-20150201.RData")
    load("~/work/modelSetC-20150201.RData")
    load("~/work/modelSetD-20150201.RData")
    load("~/work/modelSetBD-20150201.RData")
    modelSetA <- post(modelSetA)
    modelSetB <- post(modelSetB)
    modelSetC <- post(modelSetC)
    modelSetD <- post(modelSetD)
    modelSetBD <- post(modelSetBD)
    doOnce <- FALSE
}
## ## labels
##     c("1"="No screening",
##       "2"="50 only",
##       "3"="60 only",
##       "4"="70 only",
##       "5"="Opportunistic",
##       "6"="Risk stratified (2+4)",
##       "7"="Risk stratified (4+8)",
##       "8"="Mixed screening")

## modelSetBD0 <- modelSet(makeModel(discount=0,formal_compliance=0,formal_costs=0,panel=TRUE))
## modelSetBD0 <- post(modelSetBD0)
plot.scenarios(c(modelSetA,modelSetBD0),
               type="n",textp=FALSE)
points.scenarios(modelSetBD0,col="violet",textp=FALSE)
points.scenarios(modelSetA,col="red",textp=FALSE)
segments.scenarios(modelSetBD0, modelSetA,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("Panel + informal","PSA + informal"),col=c("violet","red"),bty="n",pch=19,pt.cex=1.5)

## Tables of costs and effectiveness
rbind(transform(summary.scenarios(modelSetB),set="B"),
      transform(summary.scenarios(modelSetC),set="C"),
      transform(summary.scenarios(modelSetD),set="D"))

## individual plots
pdf("~/Downloads/cea-A-LY.pdf")
plot.scenarios(modelSetA,effects="delta.LE",ylab="Effectiveness (LY)")
dev.off()
pdf("~/Downloads/cea-A-QALY.pdf")
plot.scenarios(modelSetA)
dev.off()
pdf("~/Downloads/cea-B.pdf")
plot.scenarios(modelSetB,col="red")
dev.off()
pdf("~/Downloads/cea-C.pdf")
plot.scenarios(modelSetC,col="orange")
dev.off()
pdf("~/Downloads/cea-D.pdf")
plot.scenarios(modelSetD,col="green")
dev.off()

## sets B, BD, C and D together
plot.scenarios(c(modelSetB,modelSetC,modelSetD,modelSetBD),type="n",textp=FALSE)
points.scenarios(modelSetC,col="orange")
points.scenarios(modelSetB,col="red")
points.scenarios(modelSetD,col="green",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
legend("bottomright",legend=c("Panel + formal","PSA + formal","Panel + informal","PSA + informal"),col=c("green","orange","violet","red"),bty="n",pch=19,pt.cex=1.5)

pdf("~/Downloads/cea-BC.pdf")
plot.scenarios(c(modelSetC,modelSetB),xlim=c(0,3000),
               type="n",textp=FALSE)
points.scenarios(modelSetC,col="orange",textp=FALSE)
points.scenarios(modelSetB,col="red",textp=FALSE)
segments.scenarios(modelSetB, modelSetC,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("PSA + formal","PSA + informal"),col=c("orange","red"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-BvCD.pdf")
plot.scenarios(c(modelSetBD,modelSetB),xlim=c(0,3000),
               type="n",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
points.scenarios(modelSetB,col="red",textp=FALSE)
segments.scenarios(modelSetB, modelSetBD,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("Panel + informal","PSA + informal"),col=c("violet","red"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-CD.pdf")
plot.scenarios(c(modelSetC,modelSetD),xlim=c(0,3000),col="orange",textp=FALSE,type="n")
points.scenarios(modelSetC,col="orange",textp=FALSE)
points.scenarios(modelSetD,col="green",textp=FALSE)
segments.scenarios(modelSetC, modelSetD,textp=TRUE)
legend("bottomright",legend=c("PSA + formal","Panel + formal"),col=c("orange","green"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-BDvD.pdf")
plot.scenarios(c(modelSetD,modelSetBD),xlim=c(0,3000),textp=FALSE,type="n")
points.scenarios(modelSetD,col="green",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
segments.scenarios(modelSetBD, modelSetD,textp=TRUE)
legend("bottomright",legend=c("Panel + informal","Panel + formal"),col=c("violet","green"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-incidence.pdf")
plot(modelSetA[[1]],type="incidence",xlab="Age (years)", ylab="Incidence rate", lwd=2)
for (i in 2:length(modelSetA))
    lines(modelSetA[[i]],col=i,type="incidence", lwd=2)
legend("bottomright",legend=names(modelSetA),lty=1,col=1:length(modelSetA),bty="n",lwd=2)
dev.off()


## Just compliance
model <- makeModel(discount=0,formal_compliance=1,formal_costs=0,panel=FALSE)
model0 <- model("noScreening")
model2 <- model("twoYearlyScreen50to70")
plot(model0,type="cancerdeath")
comparison1 <- rbind(predict(model2,type="cancerdeath") %>% mutate(screen=1),
                    predict(model0,type="cancerdeath") %>% mutate(screen=0)) %>%
    filter(age >= 50 & age<70)
require(splines)
exp(coef(glm(n ~ offset(log(pt)) + ns(age) + screen, data=comparison1, family=poisson)))
##
comparison <-
    inner_join(predict(model2,type="cancerdeath") %>% rename(rate2=rate) %>% select(age,rate2),
               predict(model0,type="cancerdeath") %>% rename(rate0=rate) %>% select(age,rate0)) %>%
    mutate(RR=rate2/rate0) %>% filter(age>=50 & age<90)
plot(RR ~ age, data=comparison, type="l")

## Treatment patterns
treat <- with(stockholmTreatment,
              data.frame(data.frame(year=DxY,Age=Age,G=factor(G+1,labels=c("Gleason 6","Gleason 7","Gleason 8+"))),
                         Treatment=factor(rep(1:3,each=nrow(stockholmTreatment)),labels=c("CM","RP","RT")),
                         Proportion=as.vector(cbind(CM,RP,RT))))
require(ggplot2)
pdf("~/work/treatment_patterns.pdf")
ggplot(treat, aes(x=Age,y=Proportion,group=Treatment,fill=Treatment)) + facet_wrap(~G) + geom_area(position="fill") + xlab("Age (years)")
dev.off()

if (FALSE) {
    plot(modelSetD[["No screening"]],type="incidence")
    lines(modelSetD[[""]],type="incidence",col="blue")
    lines(modelMixedScreening,type="incidence",col="red")
    lines(modelGoteborg,type="incidence",col="green")
    lines(modelRiskStratified,type="incidence",col="lightblue")
    lines(model1,type="incidence",col="orange")
    lines(modelUptake1960,type="incidence",col="pink")
}

## List of homogeneous elements
List <- function(...) {
    .Data <- list(...)
    class.element <- class(.Data[[1]])
    stopifnot(all(sapply(.Data, function(element) class(element)==class.element)))
    structure(.Data=.Data,
              element.class=class.element, # new attribute
              names=names(.Data),
              class = c("List","list"))
}
print.List <- function(obj,...) {
    i <- 1
    namess <- names(obj)
    for (obji in obj) {
        name <- if (is.null(namess)) sprintf("[[%i]]",i) else namess[i]
        cat(name,"\n")
        print(obji,...)
        i <- i+1
    }
    invisible(obj)
}
getListFunction <- function(fun,obj,...) {
    stopifnot(inherits(obj,"List"))
    VALUE <- sprintf("%s.List.%s",
                     deparse(substitute(fun)),
                     attr(obj,"element.class"))
    FUN <- tryCatch(get(VALUE,...))
    if (inherits(FUN,"try-error")) stop(sprintf("%s is not defined.\n",VALUE))
    FUN
}
plot.List <- function(obj,...) {
    getListFunction(plot,obj)(obj,...)
}
plot.List.fhcrc <- function(obj,...) {
    temp <- data.frame(t(sapply(obj,
                                function(obji) unlist(ICER(obji,obj[[1]])))))
    with(temp,
         plot(delta.costs,
              delta.QALE,
              xlim=c(0,max(delta.costs)*1.3),
              ylim=c(0,max(delta.QALE)*1.1),
              xlab="Costs (SEK)",
              ylab="Effectiveness (QALY)"))
    lines.frontier(temp$delta.costs,temp$delta.QALE)
    invisible(obj)
}
plot(List(model0,model1,model2,model50,model60,model70,
          modelUptake1930,modelUptake1960,modelGoteborg,
          modelRiskStratified,modelMixedScreening,modelFormalTestManagement))
List(model0,model1,model2,model50,model60,model70,
     modelUptake1930,modelUptake1960,modelGoteborg,
     modelRiskStratified,modelMixedScreening,modelFormalTestManagement)



## save(model0,model1,model1p,model2,model2p,model50,model60,model70,
##      modelUptake1930,modelUptake1960,modelGoteborg,
##      modelRiskStratified,modelMixedScreening,modelFormalTestManagement,
##      file="~/work/icer_20150111.RData")


mubeta2 <- c(0.0397,0.1678)
p <- 0.9 
fun <- function(par,a,b) {
    alpha <- par[1]
    beta <- par[2]
    (exp(alpha+beta*b)-exp(alpha+beta*a))/(b-a)/beta
}
objective <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    (mubeta2[1]-fun(par,0,p))^2+
    (mubeta2[2]-fun(par,p,1))^2
}
fit <- optim(c(1,1),objective,control=list(abstol=1e-16,reltol=1e-16))
fun(fit$par, p, 1)
fun(fit$par, 0, p)
p1 <- 0.6
fun(fit$par, 0, p1)
fun(fit$par, p1, p)
with(list(x=seq(0,1,length=301)),
     plot(x,sapply(x,function(xi) exp(fit$par[1]+fit$par[2]*xi)), type="l"))
abline(v=p,lty=2)
abline(v=p1, lty=2)
segments(0,mubeta2[1],p,mubeta2[1],lty=3)
segments(p,mubeta2[2],1,mubeta2[2],lty=3)

## Even width
mubeta2 <- c(0.0397,0.1678)
fun <- function(par,a,b) {
    alpha <- par[1]
    beta <- par[2]
    (exp(alpha+beta*b)-exp(alpha+beta*a))/(b-a)/beta
}
objective <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    (mubeta2[1]-fun(par,6,8))^2+
    (mubeta2[2]-fun(par,8,9))^2
}
fit <- optim(c(1,1),objective,control=list(abstol=1e-16,reltol=1e-16))
fun(fit$par, 6, 8)
fun(fit$par, 8, 9)
fun(fit$par, 6, 7)
fun(fit$par, 7, 8)
with(list(x=seq(6,9,length=301)),
     plot(x,sapply(x,function(xi) exp(fit$par[1]+fit$par[2]*xi)), type="l"))


## check cost calculations
model0 <- callFhcrc(1e5,screen="twoYearlyScreen50to70",mc.cores=3,pop=1995-50,discountRate=0)
model1 <- callFhcrc(1e5,screen="fourYearlyScreen50to70",mc.cores=3,pop=1995-50,discountRate=0)
costs <- model1$costs
pt <- model1$summary$pt
pop1 <- sqldf("select age, sum(pt) as pop from pt group by age")
costs1 <- sqldf("select age, item, sum(costs) as costs from costs group by age, item")
sqldf("select item, sum(costs/pop*IHE/1e6) as adj from pop1 natural join costs1 natural join IHEpop group by item")


## Correlated PSA values
refresh
require(microsimulation)
require(mvtnorm)
require(dplyr)
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
rmvnormPos <- function(n,mean=0,sigma=matrix(1),lbound=0) {
    x <- rmvnorm(n,mean,sigma)
    while(any(i <- which(apply(x,1,min) < lbound)))
        x[i,] <- rmvnorm(length(i),mean,sigma)
    x
}
## rmvnormPos(10,c(0,0),matrix(c(1,0,0,1),2))
prob_grade7 <- fhcrcData$prob_grade7 %>% "names<-"(c("x","y")) %>% approxfun()
psaSimCor <- function(n,age=50,rho=0.62,max.psa=50,mubeta2.scale) {
    age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
    grade <- ifelse(runif(n)>=1+FhcrcParameters$c_low_grade_slope*(age_o-35),8,7)
    cor0 <- cor1 <- cor2 <- matrix(c(1,rho,rho,1),2)
    beta0 <- with(p, rmvnorm(n,c(mubeta0,mubeta0),sebeta0^2*cor0))
    beta1 <- with(p, rmvnorm(n,c(mubeta1,mubeta1),sebeta1^2*cor1))
    beta2 <- matrix(NA,n,2)
    for (gradei in 7:8) {
        i <- which(grade == gradei)
        index <- if(gradei==7) 1 else 2
        if (any(i)) {
            x <- with(p, rmvnorm(length(i),c(mubeta2[index],mubeta2.scale*mubeta2[index]),sebeta2[index]^2*cor2))
            beta2[i,1] <- x[,1]
            beta2[i,2] <- x[,2]
        }
    }
    ext_grade <- ifelse(grade==7,
                        ifelse(runif(n)<prob_grade7(beta2),7,6),
                        8)
    eps <- with(p, cbind(rnorm(n,0,tau2),rnorm(n,0,tau2)))
    lpsa <- t(apply(beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps, 1, pmin, log(max.psa)))
    ## psa <- exp(lpsa)
    data.frame(age=age,
               cancer=age_o<age,
               advCancer=age_o<age, ### !!!!!!!
               lpsa=lpsa[,1],
               lbp=lpsa[,2],
               psa=exp(lpsa[,1]),
               bp=exp(lpsa[,2]),
               age_o=age_o,
               grade=grade,
               ext_grade=ext_grade)
}
dAgg <- function(data,threshold1,threshold2)
    mutate(data,posPSA=psa>=threshold1,posBP=bp>=threshold2) %>% group_by(advCancer,posPSA,posBP) %>% summarize(freq=n())
rTPF <- function(data) {
    a <- filter(data,advCancer & posPSA & posBP)$freq
    b <- filter(data,advCancer & !posPSA & posBP)$freq
    c <- filter(data,advCancer & posPSA & !posBP)$freq
    (a+b)/(a+c)
}
rFPF <- function(data) {
    e <- filter(data,!advCancer & posPSA & posBP)$freq
    f <- filter(data,!advCancer & !posPSA & posBP)$freq
    g <- filter(data,!advCancer & posPSA & !posBP)$freq
    (e+f)/(e+g)
}
rBiopsy <- function(data) 
    sum(filter(data,posBP)$freq)/sum(filter(data,posPSA)$freq)
RNGkind("Mersenne-Twister")
set.seed(12345)
d <- psaSimCor(10000,age=70,mubeta2.scale=2.1,rho=0.62)
plot(log(bp) ~ log(psa), data=d)
cor(subset(d,select=c(lpsa,lbp)))
## dAgg(d,3,3)
dAgg(d,3,3) %>% rTPF()
dAgg(d,3,3) %>% rFPF()

## uniroot1 <- uniroot(function(x) dAgg(d,3,x) %>% rBiopsy()-1, interval=c(1,20))
uniroot1 <- uniroot(function(x) dAgg(d,3,x) %>% rTPF()-1, interval=c(1,20))
uniroot1
## dAgg(d,3,uniroot1$root)
dAgg(d,3,uniroot1$root) %>% rTPF()
dAgg(d,3,uniroot1$root) %>% rFPF()
dAgg(d,3,uniroot1$root) %>% rBiopsy()

## Random draw from a bivariate normal distribution
rho <- 0.62
Sigma <- matrix(c(1,rho,rho,1),2)
A <- chol(Sigma)
z <- matrix(rnorm(2*1e5),nrow=1e5)
y <- cbind(z[,1],z[,1]*rho+z[,2]*sqrt(1-rho*rho))
## y <- z %*% A
cor(y)
apply(y,2,mean)
apply(y,2,sd)




lpsa1 <- psaSimCor(1e5,age=70,rho=0.0)
lpsa2 <- apply(lpsa1,1,mean)
## lpsa2 <- apply(lpsa1,1,function(x) sum(x*c(0.3,0.7)))
cor(lpsa1[,1],lpsa2) # cor>=0.71


plot(density(psaSim(1e4)$psa),xlim=c(0,20))
i <- 1
for (age in seq(55,80,by=5)) {
    lines(density(psaSim(1e5,age=age)$psa),col=i)
    i <- i+1
}





## The baseline FHCRC model assumes that PSA is an unbiased measure of the underlying diease process. The results here suggest that imprecision in the measure is less important than bias - and that PSA would need to be relatively biased to get the predicted change in biopsies from STHLM3.
## The main challenge now is that the FHCRC model was based on PCPT trial data which will not be available - nor, probably, will the bias be estimable from observed data.

## correlated betas
rho <- 0.5 # correlation
set.seed(12345+1)
beta0 <- correlatedValue(temp2$beta0,p$mubeta0,p$sebeta0, rho)
beta1 <- correlatedValue(temp2$beta1,p$mubeta1,p$sebeta1, rho)
beta2 <- pmax(0,correlatedValue(temp2$beta2,p$mubeta2[temp2$grade],p$sebeta2[temp2$grade], rho)) # should be a conditional distribution
biomarker2 <- exp(beta0+beta1*(temp2$age-35)+beta2*pos(temp2$age-35-temp2$t0)+rnorm(n,0,sqrt(p$tau2)))

## completely independent biomarker with more measurement error
set.seed(12345+1)
beta0 <- rnorm(n,p$mubeta0,p$sebeta0)
beta1 <- rnorm(n,p$mubeta1,p$sebeta1)
beta2 <- rnormPos(n,p$mubeta2[temp2$grade],p$sebeta2[temp2$grade])
biomarker2 <- exp(beta0+beta1*(temp2$age-35)+beta2*pos(temp2$age-35-temp2$t0)+rnorm(n,0,2*sqrt(p$tau2)))
plot(temp2$psa,biomarker2,log="xy")

set.seed(12345+2)
temp2 <- transform(temp2,
                   BBP=Z)
## STHLM3 simulation report
with(list(threshold=3.11),with(transform(temp2,BBPpos=(psa>=1 & BBP>=threshold),PSApos=(psa>=3)),
     cat(sprintf("
PSA+ & advanced:\t\t%i
PSA+ & cancer:\t\t\t%i
BBP+ & adv:\t\t\t%i
BBP+ & can:\t\t\t%i
BBP+ & PSA+:\t\t\t%i
BBP+ & PSA-:\t\t\t%i
BBP- & PSA+:\t\t\t%i
PSA+ | BBP+:\t\t\t%i
Prop reduction in biospies:\t%5.3f\n",
                 sum(PSApos & advanced),
                 sum(PSApos & cancer),
                 sum(BBPpos & advanced),
                 sum(BBPpos & cancer),
                 sum(BBPpos & PSApos),
                 sum(BBPpos & !PSApos),
                 sum(!BBPpos & PSApos),
                 sum(BBPpos | PSApos),
                 (sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
                       sum(BBPpos | PSApos))))))

logZ=rnorm(100000,0,0.1)
logpsa=logZ+rnorm(100000,0,sqrt(p$tau2))
Z=exp(logZ)
psa=exp(logpsa)
plot(density(logZ))
lines(density(logpsa),lty=2)
mean(logZ)
mean(logpsa)
mean(Z)
mean(psa)



## simulating sequentially from a bivariate normal distribution
set.seed(12345)
n <- 1e5
rho <- 0.6
sigma <- 2
u1 <- rnorm(n)
u2 <- rnorm(n)
x1 <- u1*sigma
x2 <- rho*u1*sigma+sqrt(1-rho^2)*u2*sigma
cor(x1,x2)
var(x1)
var(x2)

    
## testing the user-defined random number generator
init.seed <- as.integer(c(407,rep(12345,6)))
RNGkind("user")
set.user.Random.seed(init.seed)
testA <- runif(2)
next.user.Random.substream()
testB <- runif(2)
set.user.Random.seed(parallel::nextRNGStream(init.seed))
newSeed <- user.Random.seed()
testC <- runif(2)
set.user.Random.seed(parallel::nextRNGStream(newSeed))
testD <- runif(2)
##
RNGkind("L'Ecuyer-CMRG")
init.seed <- as.integer(c(407,rep(12345,6)))
.Random.seed <- init.seed
all(testA == runif(2))
.Random.seed <- parallel::nextRNGSubStream(init.seed)
all(testB == runif(2))
newSeed <- .Random.seed <- parallel::nextRNGStream(init.seed)
all(testC == runif(2))
.Random.seed <- parallel::nextRNGStream(newSeed)
all(testD == runif(2))

## Simulated likelihood
simulate <- function(i,mu=0,ntarget=100,sdsim=1,common.seed=12345,nsim=1000) {
    set.seed(i)
    y <- rnorm(ntarget,mu)
    sim <- stats::rnorm
    objective <- function(theta,common=TRUE) {
        if (common) 
            set.seed(common.seed)
        ## sum((y-mean(sim(nsim,theta,sdsim)))^2)
        (mean(y)-mean(sim(nsim,theta,sdsim)))^2
    }
    c(ymean=mean(y),
      standard=optimize(objective,lower=-10,upper=10, common=FALSE)$minimum,
      common=optimize(objective,lower=-10,upper=10, common=TRUE)$minimum)
}
set.seed(12345)
par(mfrow=c(2,2))
sims <- t(sapply(1:100,simulate,ntarget=100,nsim=100))
matplot(sims[,1],sims[,2:3], type="p",pch=1:2, xlab="Target sample mean",ylab="Mean from simulated least squares",
        main="ntarget=100, nsim=100")
abline(0,1)
legend("topleft",bty="n",legend=c("No common random numbers","Common random numbers"),pch=1:2,col=1:2)
sims <- t(sapply(1:100,simulate,ntarget=100,nsim=1000))
matplot(sims[,1],sims[,2:3], type="p",pch=1:2, xlab="Target sample mean",ylab="Mean from simulated least squares",
        main="ntarget=100, nsim=1000")
abline(0,1)
sims <- t(sapply(1:100,simulate,ntarget=1000,nsim=100))
matplot(sims[,1],sims[,2:3], type="p",pch=1:2, xlab="Target sample mean",ylab="Mean from simulated least squares",
        main="ntarget=1000, nsim=100")
abline(0,1)
sims <- t(sapply(1:100,simulate,ntarget=1000,nsim=1000))
matplot(sims[,1],sims[,2:3], type="p",pch=1:2, xlab="Target sample mean",ylab="Mean from simulated least squares",
        main="ntarget=1000, nsim=1000")
abline(0,1)



## More unit tests required

system.time(callFhcrc(1e5))
system.time(callFhcrc(1e5,mc.cores=4))
system.time(callFhcrc(1e6,mc.cores=4))

## Reading in the data from FHCRC
fhcrcData <- lapply(dir("~/src/fhcrc/data")[-10],
               function(name) structure(read.table(paste("~/src/fhcrc/data/",
                                                         name,sep=""),
                                                   head=TRUE,sep=","),
                                        filename=name))
lookup <- data.frame(filename=c("all_cause_mortality.csv",
   "biopsy_frequency.csv",
   "biopsy_sensitivity_smoothed.csv",
   "seer_incidence_imputed.csv",
   "hormone_frequency.csv",
   "primary_treatment_frequency.csv",
   "dre_sensitivity.csv",
   "gleason_7_frequency.csv",
   "stage_T2a_frequency.csv",
   "prostate_cancer_survival_local-regional.csv",
   "prostate_cancer_survival_distant.csv"),
  enum=c("all_cause_mortality",
    "biopsy_frequency",
    "biopsy_sensitivity",
    "obs_incidence",
    "pradt",
    "prtx",
    "dre",
    "prob_grade7",
    "prob_earlystage",
    "survival_local",
    "survival_dist"), stringsAsFactors = FALSE)
lookup <- subset(lookup, enum!="obs_incidence")
names(fhcrcData) <- with(lookup, enum[match(lapply(fhcrcData,attr,"filename"),
                                            filename)])
save("fhcrcData",file="~/src/R/microsimulation/data/fhcrcData.rda")
## lapply(fhcrcData,head)
##
## biopsy frequency
## with(fhcrcData[[2]],data.frame(psa=rep(PSA.beg,5),
##                           age=rep(c(55,60,65,70,75),each=3),
##                           biopsy_frequency=unlist(temp[[2]][,-(1:2)])))
## fhcrcData[[2]]

## testing using parallel
require(parallel)
require(microsimulation)
n <- 1e4
system.time(test <- mclapply(1:10,
                             function(i) callFhcrc(n,screen="noScreening"),
                             mc.cores=1))
system.time(test <- mclapply(1:10,
                             function(i) callFhcrc(n,screen="noScreening"),
                             mc.cores=4))
##
test <- lapply(1:10, function(i) callFhcrc(10,screen="noScreening"))
test2 <- list(lifeHistories=do.call("rbind", lapply(test,function(obj) obj$lifeHistories)),
              enum=test[[1]]$enum,
              n=sum(sapply(test,function(obj) obj$n)),
              parameters=do.call("rbind", lapply(test,"[[", "parameters")),
              summary=list(pt=do.call("rbind", lapply(test,function(obj) obj$summary$pt)),
                events=do.call("rbind", lapply(test,function(obj) obj$summary$events)),
                prev=do.call("rbind", lapply(test,function(obj) obj$summary$prev))))

## baseline analysis
options(width=110)
require(microsimulation)
n <- 1e7
n.cores <- 4
compliance <- 0.75
participation <- 1.0
noScreening <- callFhcrc(n,screen="noScreening",mc.cores=n.cores)
## "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified"
uptake <- callFhcrc(n,screen="screenUptake",mc.cores=n.cores,
                    studyParticipation=participation,
                    screeningCompliance=compliance)
goteborg <- callFhcrc(n,screen="stockholm3_goteborg",mc.cores=n.cores,
                      studyParticipation=participation,
                      screeningCompliance=compliance)
riskStrat <- callFhcrc(n,screen="stockholm3_risk_stratified",mc.cores=n.cores,
                       studyParticipation=participation,
                       screeningCompliance=compliance)

## Lexis diagrams
plotLexis <- function(obj) {
    stopifnot(require(Epi))
    stopifnot(require(sqldf))
    history <- obj$lifeHistories
    param <- obj$parameters
    tab <- sqldf("select t1.*, ageAtCancerDiagnosis, cohort, t0 from (select id, end as ageAtDeath, (event='toCancerDeath') as cancerDeath from history where event in ('toOtherDeath','toCancerDeath')) as t1  inner join param as p on p.id=t1.id left join (select id, end as ageAtCancerDiagnosis from history where event in ('toClinicalDiagnosis','toScreenDiagnosis')) as t2 on t1.id=t2.id")
    lexis1 <- Lexis(entry=list(coh=cohort,age=0),exit=list(coh=cohort+ageAtDeath,age=ageAtDeath),
                    data=tab)
    plot(lexis1, xlab="Calendar period", ylab="Age (year)", ylim=c(0,100), asp=1)
    with(subset(tab,!is.na(ageAtCancerDiagnosis)),
         points(cohort+ageAtCancerDiagnosis,ageAtCancerDiagnosis,pch=19,cex=0.4,col="red"))
    with(subset(tab,t0+35<ageAtDeath),
         points(cohort+t0+35,t0+35,pch=19,cex=0.2,col="blue"))
    legend("topleft",legend=c("Latent cancer onset","Cancer diagnosis"),
           pch=19,col=c("blue","red"),bty="n")
}
pdf(file="~/work/lexis-20131128.pdf",width=5,height=4)
par(mar=c(5.1, 4.1, 4.1-2, 2.1))
plotLexis(noScreening)
dev.off()

## rate calculations
pop <- data.frame(age=0:100,pop=c(12589, 14785, 15373, 14899, 14667,
14437, 14076, 13386, 13425, 12971, 12366, 11659, 11383, 10913, 11059,
11040, 11429, 12303, 13368, 13388, 13670, 13539, 13886, 13913, 14269,
14508, 15073, 15419, 15767, 15721, 16328, 16489, 17126, 16345, 15573,
15702, 16017, 16251, 17069, 16853, 16898, 16506, 15738, 15151, 15224,
15960, 16248, 16272, 16325, 14963, 14091, 13514, 13000, 12758, 12521,
12534, 12333, 11699, 11320, 11167, 11106, 10427, 10889, 10732, 11042,
11367, 11269, 11210, 10982, 10115, 9000, 7652, 6995, 6680, 6144, 5473,
5108, 4721, 4130, 3911, 3756, 3507, 3249, 2803, 2708, 2355, 2188,
2020, 1734, 1558, 1183, 1064, 847, 539, 381, 277, 185, 90, 79, 48,
61))
w <- with(subset(pop,age>=50 & age<80),data.frame(age=age,wt=pop/sum(pop)))

require(sqldf)
eventRates <- function(obj,pattern="Diagnosis") {
    stopifnot(require(sqldf))
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events natural join ev group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year")
}
plotEvents <- function(pattern,ylab="Rate",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  with(eventRates(noScreening,pattern),
       plot(year, rate, type="l",ylim=c(0,max(eventRates(goteborg,pattern)$rate)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(eventRates(uptake,pattern), lines(year, rate, col="red"))
  with(eventRates(goteborg,pattern), lines(year, rate, col="green"))
  with(eventRates(riskStrat,pattern), lines(year, rate, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol (2+2)",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
prevRatios <- function(obj,predicate) {
  ## ev <- data.frame(event=grep(pattern,levels(obj$summary$prev$event),value=TRUE))
    stopifnot(require(sqldf))
  prev <- obj$summary$prev
  sqldf(sprintf("select year, sum(n) as n, sum(y) as y, sum(p*wt) as prev from (select cohort+age as year, age, t1.n as n, coalesce(t2.y,0.0) as y, 1.0*coalesce(t2.y,0.0)/t1.n*1.0 as p from (select cohort, age, sum(count) as n from prev group by cohort, age) as t1 natural left outer join (select cohort, age, sum(count) as y from prev where %s group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year", predicate))
}
plotPrev <- function(pattern,ylab="Prevalence",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  with(prevRatios(noScreening,pattern),
       plot(year, prev, type="l",ylim=c(0,max(prevRatios(goteborg,pattern)$prev)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(prevRatios(uptake,pattern), lines(year, prev, col="red"))
  with(prevRatios(goteborg,pattern), lines(year, prev, col="green"))
  with(prevRatios(riskStrat,pattern), lines(year, prev, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol (2+2)",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
table(goteborg$summary$events$event)
table(goteborg$summary$prev$dx)

##path <- function(filename) sprintf("/media/sf_C_DRIVE/usr/tmp/tmp/%s",filename)
##pdf(path("screening_20130425.pdf"),width=7,height=6)
##par(mfrow=c(2,2))
pdf(file="~/work/screening-comparison.pdf",width=6,height=5)
plotEvents("^toScreen$",main="PSA screen",legend.x="topleft")
dev.off()
pdf(file="~/work/biopsy-comparison.pdf",width=6,height=5)
plotEvents("Biopsy",main="Biopsies",legend.x="topleft")
dev.off()
pdf(file="~/work/diagnosis-comparison.pdf",width=6,height=5)
plotEvents("Diagnosis",main="Prostate cancer incidence",legend.x="topleft")
dev.off()
pdf(file="~/work/clinicaldx-comparison.pdf",width=6,height=5)
plotEvents("^toClinicalDiagnosis$",legend.x="bottomleft",
           main="PC incidence (Clinical Dx)")
dev.off()
pdf(file="~/work/prevalence-comparison.pdf",width=6,height=5)
plotPrev("dx!='NotDiagnosed'",main="PC diagnosis",legend.x=2010,legend.y=0.04)
dev.off()
pdf(file="~/work/mortality-comparison.pdf",width=6,height=5)
plotEvents("^toCancerDeath$",legend.x="bottomleft",main="PC mortality")
dev.off()
##dev.off()

## extend the plots to include general conditions
eventRatesCondition <- function(obj,condition,substitute.condition=FALSE) {
  if (substitute.condition)
    condition <- substitute(condition)
  events <- eval(substitute(subset(obj$summary$events, condition), list(condition=condition)))
  pt <- obj$summary$pt
  sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year")
}
plotEventsCondition <- function(condition,ylab="Rate",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  condition <- substitute(condition)
  with(eventRatesCondition(noScreening,condition),
       plot(year, rate, type="l",ylim=c(0,max(eventRatesCondition(goteborg,condition)$rate)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(eventRatesCondition(uptake,condition), lines(year, rate, col="red"))
  with(eventRatesCondition(goteborg,condition), lines(year, rate, col="green"))
  with(eventRatesCondition(riskStrat,condition), lines(year, rate, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
plotEventsCondition(grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"),
                    legend.x="bottomleft",main="PC incidence Gleason 7+")

## How many PSA tests, biopsies etc in the eight years of follow-up?
ratio <- 35000/sum(subset(goteborg$summary$prev,year==2013)$count)
lastYear <- 2014+8
describe <- function(a,b)
  sprintf("pchange=%.1f%%,change=%.1f",100*(1-b/a),(a-b)*ratio)
## PSA tests
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("toScreen",event))$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("toScreen",event))$n))
## biopsies
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("toBiopsy",event))$n),
         sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("toBiopsy",event))$n))
## Cancer Gleason 6
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade == "Gleason_le_6")$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade == "Gleason_le_6")$n))
## Cancer Gleason 7+
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"))$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"))$n))
## metastatic cancer
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & state=="Metastatic")$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & state == "Metastatic")$n))



## In summary, compared with the modified Göteborg protocol over eight years of follow-up, the risk-stratified protocol is expected to have 30% fewer PSA tests (approximately 15,000 fewer), 10% fewer biopsies (~2000 fewer), 3% fewer prostate cancer diagnoses for Gleason 6 cancers (~40 fewer) and 2% fewer cancer diagnoses for Gleason 7+ cancers (~15 fewer cases).
##
## These results are very similar to those obtained by Gulati et al (2013?).
  

  
plotPrev("dx='NotDiagnosed' and state!='Healthy'",main="Latent disease",
         legend.x=2010,legend.y=0.04)
plotPrev("dx='NotDiagnosed' and state!='Healthy' and psa='PSA>=3'",
         main="Latent screen-detectable disease",
         legend.x=2010,legend.y=0.04)

temp0 <- subset(uptake$summary$prev,year==2010 & dx=='NotDiagnosed')
temp <- subset(temp0,state!='Healthy' & psa=='PSA>=3')
temp2 <- sqldf("select pop.age,pop,coalesce(n,0) as n,coalesce(n,0)*1.0/pop as prev  from
(select age, sum(n) as pop from temp0 group by age) as pop natural left join
(select age, sum(n) as n from temp group by age) as cases")
##
temp <- subset(temp0,psa=='PSA>=3')
temp3 <- sqldf("select pop.age,pop,coalesce(n,0) as n,coalesce(n,0)*1.0/pop as prev  from
(select age, sum(n) as pop from temp0 group by age) as pop natural left join
(select age, sum(n) as n from temp group by age) as cases")
with(temp3, plot(age,prev,type="l",ylim=c(0,0.55)))
with(temp2, lines(age,prev,lty=2))

w <- with(subset(pop,age>=50 & age<70),data.frame(age=age,wt=pop/sum(pop)))

sqldf("select sum(wt*prev)/sum(wt) from temp3 natural join w") # prev of PSA 3+ | Not diagnosed
sqldf("select sum(wt*prev)/sum(wt) from temp2 natural join w") # prev of PSA 3+ & cancer | Not diagnosed
5e4*sqldf("select sum(wt*prev)/sum(wt) from temp3 natural join w")
5e4*sqldf("select sum(wt*prev)/sum(wt) from temp2 natural join w")



## Plot of the cohorts over the Lexis diagram
plot(c(1900,2030),c(0,100),type="n",xlab="Calendar period",ylab="Age (years)")
polygon(c(1900,1980,1980+50,1980+50,1980+50-30,1900),
        c(0,0,50,100,100,0))
polygon(c(1990,2030,2030,1990),
        c(50,50,80,80),
        lty=2, border="blue")



### FHCRC model ###
options(width=110)
require(microsimulation)
n <- 1e6
n.cores <- 4
temp <- callFhcrc(n,screen="noScreening",mc.cores=n.cores)
temp4 <- callFhcrc(n,screen="fourYearlyScreen50to70",mc.cores=n.cores)
temp2 <- callFhcrc(n,screen="twoYearlyScreen50to70",mc.cores=n.cores)
temp50 <- callFhcrc(n,screen="screen50",mc.cores=n.cores)
temp60 <- callFhcrc(n,screen="screen60",mc.cores=n.cores)
temp70 <- callFhcrc(n,screen="screen70",mc.cores=n.cores)
uptake <- callFhcrc(n,screen="screenUptake",mc.cores=n.cores)
## "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified"

## incremental life-expectancy calculations
LE <- function(obj) sum(obj$summary$pt$pt)/obj$n
IE <- function(obj,objref=temp) LE(obj)-LE(objref)
LE(temp)
IE(temp2)
IE(temp4)
IE(temp50)
IE(temp60)
IE(temp70)

require(data.table)
prev <- data.table(temp2$summary$prev,key="age")
totals <- prev[,sum(count),by="age"]
strat <- prev[,sum(count),by="age,state"]
m <- transform(merge(totals,strat,all=TRUE),prev=V1.y/V1.x)
plot(prev~age+state,m)


require(sqldf)
eventRatesOld <- function(obj,pattern="Diagnosis") {
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select year, age, sum(pt) as pt from pt group by year, age) as t1 natural left outer join (select year, age, sum(n) as n from events natural join ev group by year, age) as t2")
}
dx <- eventRatesOld(uptake)
require(mgcv)
dx$pred <- 1000*predict(gam(n~s(age,year,k=50),data=dx,offset=log(pt),family=poisson),type="response")
library("RColorBrewer"); library("lattice");
brewer.div <- colorRampPalette(brewer.pal(9, "Spectral"), interpolate = "spline")
pdf("~/work/levelplot-pc-incidence.pdf")
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = rev(brewer.div(100)), aspect = "iso",
  xlab="Age (years)", ylab="Calendar period", ylim=c(1980,2050))
dev.off()
with(list(res=600),
       jpeg(file="~/work/levelplot-pc-incidence.jpg",height=5*res,width=5*res,res=res,
            quality=100))
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = rev(brewer.div(100)), aspect = "iso",
  xlab="Age (years)", ylab="Calendar period", ylim=c(1980,2050))
dev.off()



## State occupancy: prevalence in different states
prevRatios <- function(obj,predicate) {
  ## ev <- data.frame(event=grep(pattern,levels(obj$summary$prev$event),value=TRUE))
    stopifnot(require(sqldf))
  prev <- obj$summary$prev
  sqldf(sprintf("select year, sum(n) as n, sum(y) as y, sum(p*wt) as prev from (select cohort+age as year, age, t1.n as n, coalesce(t2.y,0.0) as y, 1.0*coalesce(t2.y,0.0)/t1.n*1.0 as p from (select cohort, age, sum(count) as n from prev group by cohort, age) as t1 natural left outer join (select cohort, age, sum(count) as y from prev where %s group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year", predicate))
}
##plotPrev("dx!='NotDiagnosed'",main="PC diagnosis",legend.x=2010,legend.y=0.04)

## rate calculations
## do this using data.table
require(data.table)
eventRates <- function(obj,pattern="Diagnosis") {
  events <- data.table(obj$summary$events,key="event")
  ev <- grep(pattern,levels(events$event),value=TRUE)
  pt <- data.table(obj$summary$pt,key="age")
  with(merge(pt[,sum(pt),by=age],events[J(ev),sum(n),by=age], all=TRUE),
       transform(data.table(age=age,pt=V1.x,n=ifelse(is.na(V1.y),0.0,V1.y))[-1,],
                 rate=n/pt))
}
## the old way: using SQL
require(sqldf)
eventRatesOld <- function(obj,pattern="Diagnosis") {
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select age, sum(pt) as pt from pt group by age) as t1 natural left outer join (select age, sum(n) as n from events natural join ev group by age) as t2")
}

require(microsimulation)
require(dplyr)
require(splines)
base <- callFhcrc(1e6,screen="noScreening",mc.cores=3,cohort=1970)
new <- callFhcrc(1e6,screen="twoYearlyScreen50to70",mc.cores=3,cohort=1970)
baserate <- predict(base,"cancerdeath") %>% filter(age>=50)
newrate <- predict(new,"cancerdeath") %>% filter(age>=50)
merged <- rbind(transform(baserate,group=0), transform(newrate,group=1)) %>% transform(age50=age-50)
fit1 <- glm(n ~ ns(age,5)+group:ns(age,5)+offset(log(pt)), data=merged, family=poisson)
RR <- predict(fit1,newdata=transform(baserate,group=1,pt=1,age50=age-50),type="response") /
predict(fit1,newdata=transform(baserate,group=0,pt=1,age50=age-50),type="response")

pdf("~/work/mortality_rate_reduction_twoYearly50to69.pdf",width=7,height=4)
par(mfrow=1:2)
with(predict(new,"incidence") %>% filter(age>=40),
     plot(age,rate*1e5,type="l",xlab="Age (years)",ylab="Prostate cancer incidence rate per 100,000",main="(a)"))
legend("topleft",legend=c("No screening","Screening"),lty=2:1, bty="n")
with(predict(base,"incidence") %>% filter(age>=40),
     lines(age,rate*1e5,lty=2))
plot(baserate$age, RR, type="l",ylab="Prostate cancer mortality rate ratio",xlab="Age (years)",
ylim=c(0.5,1),main="(b)")
dev.off()

with(predict(new,"cancerdeath") %>% filter(age>=40),
     plot(age,rate*1e5,type="l",xlab="Age (years)",ylab="Rate per 100,000"))
with(predict(base,"cancerdeath") %>% filter(age>=40),
     lines(age,rate*1e5,lty=2))



## all(abs(eventRates(temp)-data.table(eventRatesOld(temp)))<1e-8)
## system.time(eventRatesOld(temp))
## system.time(eventRates(temp))

png("~/work/screening-comparison-20130222.png",height=2,width=4,res=1200,units="in",pointsize=3)
##x11(width=8,height=5)
##layout(matrix(1:2,nrow=1,byrow=TRUE))
par(mfrow=c(1,2),
  mar      = c(5+2, 4+2, 4+2, 1+2)+0.1,
  ##xaxs     = "i",
  ##yaxs     = "i",
  cex.main = 2,
  cex.axis = 2,
  cex.lab  = 2
)
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Prostate cancer incidence rate",
                            main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly\nscreening"), lty=1, col=c("black","red"), bty="n",
       cex=2)
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Prostate cancer incidence rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly\nscreening"), lty=1, col=c("black","blue"), bty="n",
       cex=2)
dev.off()


pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly screening"), lty=1, col=c("black","red"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly screening"), lty=1, col=c("black","blue"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 50"))
with(eventRates(temp50), lines(age, rate, col="green"))
legend("topleft", legend=c("No screening","Screening at age 50"), lty=1, col=c("black","green"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 60"))
with(eventRates(temp60), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 60"), lty=1, col=c("black","orange"), bty="n")
dev.off()

pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly screening"), lty=1, col=c("black","red"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly screening"), lty=1, col=c("black","blue"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 50"))
with(eventRates(temp50), lines(age, rate, col="green"))
legend("topleft", legend=c("No screening","Screening at age 50"), lty=1, col=c("black","green"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 60"))
with(eventRates(temp60), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 60"), lty=1, col=c("black","orange"), bty="n")
dev.off()

with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 70"))
with(eventRates(temp70), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 70"), lty=1, col=c("black","orange"), bty="n")

personTime <- function(obj) {
  pt <- obj$summary$pt
  n <- obj$n
  sum(pt$pt)/n
}
sapply(list(temp=temp,temp2=temp2,temp4=temp4,temp50=temp50,temp60=temp60,temp70=temp70),personTime) - personTime(temp)
sapply(list(temp=temp,temp2=temp2,temp4=temp4,temp50=temp50,temp60=temp60,temp70=temp70),personTime)




1-exp(-sum(subset(eventRates(temp),age<=75)$rate)) ## 8.4% cumulative risk to age 75 years
1-exp(-sum(subset(eventRates(temp2),age<=75)$rate)) ## 9.4% cumulative risk to age 75 with one random screen between 50 and 70
1-exp(-sum(subset(eventRates(temp3),age<=75)$rate)) ## 9.4% cumulative risk to age 75 with one random screen between 50 and 70

## Stockholm lan, males, 2012

pop <- data.frame(age=0:100,pop=c(12589, 14785, 15373, 14899, 14667,
14437, 14076, 13386, 13425, 12971, 12366, 11659, 11383, 10913, 11059,
11040, 11429, 12303, 13368, 13388, 13670, 13539, 13886, 13913, 14269,
14508, 15073, 15419, 15767, 15721, 16328, 16489, 17126, 16345, 15573,
15702, 16017, 16251, 17069, 16853, 16898, 16506, 15738, 15151, 15224,
15960, 16248, 16272, 16325, 14963, 14091, 13514, 13000, 12758, 12521,
12534, 12333, 11699, 11320, 11167, 11106, 10427, 10889, 10732, 11042,
11367, 11269, 11210, 10982, 10115, 9000, 7652, 6995, 6680, 6144, 5473,
5108, 4721, 4130, 3911, 3756, 3507, 3249, 2803, 2708, 2355, 2188,
2020, 1734, 1558, 1183, 1064, 847, 539, 381, 277, 185, 90, 79, 48,
61))
expected <- function(obj,pattern="Diagnosis",start,end)
sum(transform(subset(merge(obj,eventRates(temp)),age>=50 & age<70),e=rate*pop)$e)



sum(transform(subset(merge(pop,eventRates(temp2)),age>=50 & age<70),e=rate*pop)$e) ## n=748 cases aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp3)),age>=50 & age<70),e=rate*pop)$e) ## n=748 cases aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## at present, biopsies == cancers
sum(transform(subset(merge(pop,eventRates(temp2,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## n=1080 biopsies aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp3,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## n=6920 biopsies aged 50-69 years


plot.EventReport <- function(eventReport, ...) {
  data <- transform(merge(eventReport$pt,eventReport$events), rate=n/pt)
  data <- data[with(data,order(state,event,age)),]
  with(data, plot(age,rate,type="n",ylim=c(0,0.2), ...))
  set <- unique(subset(data,select=c(state,event)))
  invisible(lapply(1:nrow(set),
                   function(i)
                   with(subset(data, state==set$state[i] & event==set$event[i]),
                        lines(age,rate,lty=i))))
} ## Rates only. TODO: add prevalence.
plot.EventReport(temp$summary, xlab="Age (years)")

sapply(temp$parameters,summary)

summary.EventReport <- function(eventReport) {
  data <- transform(merge(eventReport$pt,eventReport$events), rate=n/pt)
  data[with(data,order(state,event,age)),]
}
head(summary.EventReport(temp$summary))


## totals <- data.frame(xtabs(n~age,temp$prev))
## names(totals)[2] <- "total"
## bystate <- data.frame(xtabs(n~age+state,temp$prev))
## names(bystate)[3] <- "n"
## merge(bystate,totals)

require(sqldf)
prev <- temp$prev
temp2 <- sqldf("select *, 1.0*n/total as prev from (select age, sum(n) as total from prev group by age) as t1 left outer join prev on t1.age=prev.age order by state, age")
with(temp2, plot(age,prev,type="n"))
set <- unique(subset(temp2,select=c(state)))
invisible(lapply(1:nrow(set),
                 function(i)
                 with(subset(temp2, state==set$state[i]),
                      lines(age,prev,lty=i))))
legend("bottomleft", legend=set$state, lty=1:7, bty="n")


set.seed(123)
system.time(df <- callSimplePerson(100000))


oldRNGkind <- RNGkind()
if (exists(".Random.seed")) old.Random.seed <- .Random.seed

RNGkind("user")
df <- callSimplePerson()

RNGkind("user")
df <- callPersonSimulation(n=100)

if (exists("old.Random.seed")) .Random.seed <- old.Random.seed
do.call("RNGkind",as.list(oldRNGkind))


require(microsimulation)
Simulation <-
  setRefClass("Simulation",
              contains = "BaseDiscreteEventSimulation",
              fields = list(id = "numeric", state = "character", report = "data.frame"),
              methods= list(initialize = function(id = 0) callSuper(id = id)))
Simulation$methods(init = function() {
  clear()
  id <<- id + 1
  state <<- "Healthy"
  scheduleAt(rweibull(1,8,85), "Death due to other causes")
  scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
})
Simulation$methods(handleMessage = function(event) {
  report <<- rbind(report, data.frame(id = id,
                                      state = state,
                                      begin = previousEventTime,
                                      end = currentTime,
                                      event=event,
                                      stringsAsFactors = FALSE))
  if (event %in% c("Death due to other causes", "Cancer death")) {
    clear()
  }
  else if (event == "Cancer diagnosis") {
    state <<- "Cancer"
    if (runif(1) < 0.5)
      scheduleAt(now() + rweibull(1,2,10), "Cancer death")
  }
})
RNGkind("Mersenne-Twister")
if (exists(".Random.seed")) rm(.Random.seed)
set.seed(123)
sim <- Simulation$new()
system.time(for (i in 1:1000) sim$run())
subset(sim$report,id<=4)

RNGkind("Mersenne-Twister") # cf. "L'Ecuyer-CMRG"!
set.seed(123)
head(.Random.seed)
rng1 <- RNGStream(nextStream = FALSE)
rng2 <- RNGStream()
with(rng1,rexp(1))
with(rng2,rexp(1))
rng1$nextSubStream()
with(rng1,rexp(1))
##
rng1$resetStream()
rng2$resetStream()
with(rng1,rexp(2))
with(rng2,rexp(2))
rng1$nextSubStream()
with(rng1,rexp(2))
rng1$resetRNGkind() # be a good citizen
head(.Random.seed)



temp <- callSimplePerson2(100)
temp2 <- transform(merge(temp$pt,temp$events), rate=n/pt)
temp2 <- temp2[with(temp2,order(state,event,age)),]
with(temp2, plot(age,rate,type="n"))
set <- unique(subset(temp2,select=c(state,event)))
invisible(lapply(1:nrow(set),
                 function(i)
                 with(subset(temp2, state==set$state[i] & event==set$event[i]),
                      lines(age,rate,lty=i))))


muWeibull <- function(a,b) b*gamma(1+1/a)
varWeibull <- function(a,b)  b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
bWeibull <- function(a,mu) mu/gamma(1+1/a)
plotWeibull <- function(a,b,max=60) { x <- seq(0,max,length=301); plot(x,dweibull(x,a,b),type="l") }
##
muWeibull(2,10)
sqrt(varWeibull(2,10))
plotWeibull(2,10)
##
muWeibull(2,3)
sqrt(varWeibull(2,3))
plotWeibull(2,3)
