/**
 * @file
 * @author  Mark Clements <mark.clements@ki.se>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * Prostate cancer screening simulation model.
 */


#include <microsimulation.h>

#include <boost/algorithm/cxx11/iota.hpp>

namespace fhcrc_example {

  using namespace std;
  using namespace Rcpp;
  using namespace ssim;

  // declarations

  namespace base {
    enum grade_t {Gleason_le_7,Gleason_ge_8,Healthy};
  }
  namespace ext {
    enum grade_t {Gleason_le_6,Gleason_7,Gleason_ge_8,Healthy};
    enum state_t {Healthy_state,T1_T2,T3plus,Metastatic};
  }

  enum state_t {Healthy,Localised,Metastatic}; // stage?

  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};

  enum event_t {toLocalised, toMetastatic, toClinicalDiagnosis, toCancerDeath,
                toOtherDeath, toScreen, toBiopsyFollowUpScreen,
                toScreenInitiatedBiopsy, toClinicalDiagnosticBiopsy,
                toScreenDiagnosis, toOverDiagnosis, toOrganised, toTreatment,
                toCM, toRP, toRT, toADT, toUtilityChange, toUtilityRemove,
                toSTHLM3, toOpportunistic, toT3plus, toCancelScreens,
                toYearlyActiveSurveillance, toYearlyPostTxFollowUp, toMRI, toPalliative, toTerminal};

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70,
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified,
		 goteborg, risk_stratified, mixed_screening, regular_screen, single_screen,
                 introduced_screening_only, introduced_screening_preference, introduced_screening,
		 stopped_screening, cap_control, cap_study, sthlm3_mri_arm,
		 grs_stratified, grs_stratified_age};

  enum treatment_t {no_treatment, CM, RP, RT};

  enum survival_t {StageShiftBased, LeadTimeBased};

  enum biomarker_model_t {random_correction, psa_informed_correction};

  enum cost_t {Direct,Indirect};

  enum utility_scale_t {UtilityAdditive, UtilityMultiplicative, UtilityMinimum};

  namespace FullState {
    typedef boost::tuple<short,short,short,bool,double> Type;
    enum Fields {ext_state, ext_grade, dx, psa_ge_3, cohort};
    // string names[5] = {"ext_state","ext_grade","dx","psa_ge_3","cohort"};
  }
  namespace LifeHistory {
    typedef boost::tuple<int, short, short, int, short, double, double, double, double, double> Type;
    enum Fields {id, ext_state, ext_grade, dx, event, begin, end, year, psa, utility};
  }

  RcppExport SEXP rllogis_(SEXP shape, SEXP scale) {
    RNGScope scope;
    return wrap(R::rllogis(as<double>(shape),as<double>(scale)));
  }
  RcppExport SEXP rllogis_trunc_(SEXP shape, SEXP scale, SEXP left) {
    RNGScope scope;
    return wrap(R::rllogis_trunc(as<double>(shape),as<double>(scale),as<double>(left)));
  }

  typedef std::pair<int,string> CostKey; // (cost_type,cost_name)
  // typedef boost::tuple<int,string,double> CostKey; // (cost_type,cost_name,cohort)
  class SimOutput {
  public:
    EventReport<FullState::Type,short,double> report;
    EventReport<int,short,double> shortReport;
    CostReport<CostKey> costs;
    vector<LifeHistory::Type> lifeHistories;
    SimpleReport<double> outParameters;
    SimpleReport<double> psarecord, bxrecord, falsePositives;
    SimpleReport<double> diagnoses;
    Means tmc_minus_t0;
  };
  // SimOutput * out; // in callFhcrc
  // &out in the Person object
  // out->report etc


  typedef std::pair<double,double> Double;
  typedef Table<double,double,int,double> TablePrtx; // Age, DxY, G
  typedef Table<double,int,int,double> TableLocoHR; // Age, G, PSA10+
  typedef Table<double,double> TableDD; // Age
  typedef Table<double,double> TableMetastaticHR; // Age
  typedef Table<int,double,double,int,double> TablePradt;
  typedef Table<double,double,double> TableBiopsyCompliance;
  typedef Table<double,double,double> TableDDD; // as per TableBiopsyCompliance
  typedef map<int,NumericInterpolate> H_dist_t;
  typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  class SimInput {
  public:
    TableLocoHR hr_locoregional;
    TableMetastaticHR hr_metastatic;
    TableDD tableBiopsySensitivity, tableSecularTrendTreatment2008OR,
      tableNegBiopsyToPSAmeanlog, tableNegBiopsyToPSAsdlog, tableNegBiopsyToBiopsymeanlog,
      tableNegBiopsyToBiopsysdlog, tableCMtoRPpnever, tableCMtoRPmeanlog, tableCMtoRPsdlog,
      tableCMtoRTpnever, tableCMtoRTmeanlog, tableCMtoRTsdlog;
    TablePrtx prtxCM, prtxRP;
    TablePradt pradt;
    TableBiopsyCompliance tableOpportunisticBiopsyCompliance, tableFormalBiopsyCompliance;
    TableDDD rescreen_shape, rescreen_scale, rescreen_cure;
    NumericInterpolate interp_prob_grade7;
    H_dist_t H_dist;
    H_local_t H_local;
    set<double,greater<double> > H_local_age_set;

    Rng * rngNh, * rngOther, * rngScreen, * rngTreatment, * rngSurv, * rngBx, * rngPrelude;
    Rpexp rmu0;

    NumericVector parameter;
    LogicalVector bparameter;

    // read in the parameters
    NumericVector cost_parameters, utility_estimates, utility_duration, lost_production_years;
    NumericVector mubeta2, sebeta2; // otherParameters["mubeta2"] rather than as<NumericVector>(otherParameters["mubeta2"])
    int screen, nLifeHistories;
    bool panel, debug;
    Table<double,double> production;
    // DataFrame background_utilities;
    NumericVector bg_lower, bg_upper, bg_utility;
    // cumulative hazard for screening uptake (currently only cap_control and cap_study)
    NumericInterpolate H_screen_uptake;
    NumericVector cap_pScreened;

    ~SimInput() {
      if (rngNh != NULL) delete rngNh;
      if (rngOther != NULL) delete rngOther;
      if (rngScreen != NULL) delete rngScreen;
      if (rngTreatment != NULL) delete rngTreatment;
      if (rngSurv != NULL) delete rngSurv;
      if (rngBx != NULL) delete rngBx;
      if (rngPrelude != NULL) delete rngPrelude;
    }
  };
  // SimInput in; // callFhcrc
  // SimInput* in; // Person
  // in->hr_locoregional(...);

  class cMessageUtility : public cMessage {
  public:
    cMessageUtility(short kind, int id, double utility = 1.0) :
      cMessage(kind), id(id), utility(utility) { }
    int id;
    double utility;
  };
  class Utilities {
  public:
    typedef boost::unordered_map<int,double> UMap;
    UMap umap;
    int counter;
    utility_scale_t scale;
    bool truncate;
    Utilities(utility_scale_t scale = UtilityMultiplicative, bool truncate = true) :counter(0), scale(scale), truncate(truncate) {}
    double utility() {
      // case: no utilities?
      // case: value>1.0?
      double value = 1.0;
      for (UMap::iterator it = umap.begin(); it!=umap.end(); it++) {
	if (scale == UtilityAdditive)
	  value -= (1.0 - it->second);
	else if (scale == UtilityMultiplicative)
	  value *= it->second;
	else if (scale == UtilityMinimum) {
	  if (it->second < value) value = it->second;
	}
      }
      if (truncate && value < 0.0) value = 0.0;
      return value;
    }
    void clear() { umap.clear(); counter = 0; }
    void handleMessage(const cMessage* msg) {
      const cMessageUtility * umsg;
      if ((umsg = dynamic_cast<const cMessageUtility *>(msg)) != 0) {
	if (umsg->kind == toUtilityChange)
	  umap[umsg->id] = umsg->utility;
	else if (umsg->kind == toUtilityRemove)
	  umap.erase(umsg->id);
      }
    }
  };
  void Rprint(Utilities utilities) {
    Rprintf("utilities: [");
    for (Utilities::UMap::iterator it = utilities.umap.begin();
	 it!=utilities.umap.end();
	 it++) {
      if (it != utilities.umap.begin()) Rprintf(",");
      Rprintf("%f",it->second);
    }
    Rprintf("]\n");
  }

  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
  }

  class FhcrcPerson : public cProcess
  {
  public:
    SimInput* in;
    SimOutput* out;
    Utilities* utilities;
    double beta0, beta1, beta2;
    double t0, y0, t3p, tm, tc, tmc, aoc;
    state_t state;
    ext::state_t ext_state;
    diagnosis_t dx;
    base::grade_t future_grade, grade;
    ext::grade_t future_ext_grade, ext_grade;
    treatment_t tx;
    bool adt;
    double txhaz, psa_last_screen;
    int id, index;
    double cohort, rescreening_frailty, ageEntry, grs_frailty, other_frailty, ageFirstScreen;
    bool everPSA, previousNegativeBiopsy, organised, previousFollowup, MRIpos, everGRS;
    FhcrcPerson(SimInput* in, SimOutput* out, Utilities* utilities, const int id = 0, const double cohort = 1950, const int index = 0) :
      in(in), out(out), utilities(utilities), id(id), index(index), cohort(cohort) { };
    double utility() { return utilities->utility(); }
    double psamean(double age);
    double psameasured(double age);
    treatment_t calculate_treatment(double u, double age, double year);
    double calculate_mortality_hr(double age_diag);
    double calculate_survival(double u, double age_diag, double age_c, treatment_t tx);
    double calculate_transition_time(double u, double t_enter, double gamma);
    void opportunistic_rescreening(double psa);
    void opportunistic_uptake_if_ever();
    bool screening_preference();
    double callenderStartAge(double frailty, double threshold=0.04);
    void cancel_events_after_diagnosis();
    void rescreening_schedules(double psa, bool organised, bool mixed_programs);
    bool detectable(double now, double year);
    void init();
    void add_costs(string item, cost_t cost_type = Direct, double weight=1.0);
    void lost_productivity(string item, double weight=1.0);
    virtual void handleMessage(const cMessage* msg);
    void scheduleUtilityChange(double at, std::string category);
    void scheduleUtilityChange(double at, double utility);
    void scheduleUtilityChange(double from, double to, double utility);
    bool onset_p();
    double onset();
  };

  /**
      Calculate the (geometric) mean PSA value at a given age (** NB: this used to be t=age-35.0 **)
  */
  double FhcrcPerson::psamean(double age) {
    double t = age<35.0 ? 0.0 : age - 35.0;
    double yt = t<t0 ? exp(beta0+beta1*t) : exp(beta0+beta1*t+beta2*(t-t0));
    return yt;
  }

  /**
      Calculate the *measured* PSA value at a given age (** NB: this used to be t=age-35 **)
  */
  double FhcrcPerson::psameasured(double age) {
    return FhcrcPerson::psamean(age)*exp(R::rnorm(0.0, sqrt(double(in->parameter["tau2"]))));
    }

  /**
      Report on costs for a given item
  */
  void FhcrcPerson::add_costs(string item, cost_t cost_type, double weight) {
    out->costs.add(CostKey((int) cost_type,item),now(),in->cost_parameters(item) * weight,
		   index);
  }

  /**
      Report on lost productivity
  */
  void FhcrcPerson::lost_productivity(string item, double weight) {
    double loss = in->lost_production_years[item] * in->production(now()) * weight;
    out->costs.add(CostKey((int) Indirect,item),now(),loss,index);
  }

  /**
     Schedule a utility change.
   **/
  void FhcrcPerson::scheduleUtilityChange(double at, std::string category) {
    utilities->counter++; // increment
    scheduleAt(at, new cMessageUtility(toUtilityChange,
				       utilities->counter,
				       in->utility_estimates[category]));
    scheduleAt(at + in->utility_duration[category],
	       new cMessageUtility(toUtilityRemove, utilities->counter));
  }
  void FhcrcPerson::scheduleUtilityChange(double from, double to, double utility) {
    utilities->counter++; // increment
    scheduleAt(from, new cMessageUtility(toUtilityChange, utilities->counter, utility));
    scheduleAt(to, new cMessageUtility(toUtilityRemove, utilities->counter));
  }

  treatment_t FhcrcPerson::calculate_treatment(double u, double age, double year) {
    double pCM, pRP, pRT;
    // treatment probabilities in 2008
    if (in->bparameter["stockholmTreatment"]) {
       pCM = in->prtxCM(bounds<double>(age,50.0,85.0),
		    bounds<double>(year,2008.0,2012.0),
		    int(ext_grade));
       pRP = in->prtxRP(bounds<double>(age,50.0,85.0),
		    bounds<double>(year,2008.0,2012.0),
		    int(ext_grade));
    } else { // original FHCRC table prtx
       pCM = in->prtxCM(bounds<double>(age,50.0,79.0),
		    bounds<double>(year,1973.0,2004.0),
		    int(grade));
       pRP = in->prtxRP(bounds<double>(age,50.0,79.0),
		    bounds<double>(year,1973.0,2004.0),
		    int(grade));
    }
    pRT = 1.0 - pCM - pRP;
    // adjust for secular trend in odds of RP or RT
    {
      double OR = in->tableSecularTrendTreatment2008OR(bounds(year,1988.0,2008.0));
      double gamma = log(OR);
      double betaRT = log(pRT/pCM);
      double betaRP = log(pRP/pCM);
      pCM = 1.0/(1.0+exp(betaRT+gamma)+exp(betaRP+gamma));
      pRP = exp(betaRP+gamma)*pCM;
      pRT = exp(betaRT+gamma)*pCM;
    }
    treatment_t tx = (u<pCM) ? CM :
      (u<pCM+pRP) ? RP : RT;
    if (in->debug) Rprintf("id=%i, Age=%3.0f, DxY=%4.0f, stage=%i, grade=%i, tx=%i, u=%8.6f, pCM=%8.6f, pRP=%8.6f\n",id,age,year,state,grade,(int)tx,u,pCM,pRP);
    return tx;
  }

  /** @brief calculate survival taking account of screening and treatment.

      For the Nordic natural history model, we have calibrated
      survival to observed survival from PCBase. These calibrations
      are represented by the hr_locoregional and hr_metastatic tables.

      Note that we do not use now(), as the time points change
      depending on how we represent screening. This explains why we
      pass all of the ages as parameters.

      The ustar is an approach to adjust survival for different
      factors. Rather than solve a revised survival function, we solve
      the baseline survival function for a revised ustar = u^(1/HR).
   **/

  double FhcrcPerson::calculate_mortality_hr(double age_diag) { // also: tm, ext_grade
    double age_m = tm + 35.0;   // age at onset of metastatic cancer
    bool localised = (age_diag < age_m);
    double mort_hr;
    if (localised) {
      mort_hr = in->hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, psamean(age_diag)>10 ? 1 : 0);
      if (ext_state == ext::T3plus) mort_hr*=double(in->parameter["RR_T3plus"]);
    }
    else
      mort_hr = in->hr_metastatic(age_diag);
    return mort_hr;
  }

  double FhcrcPerson::calculate_survival(double u, double age_diag, double age_c, treatment_t tx) { // also: tc, tm, tmc, grade, now()
    double age_d = -1.0;        // age at death (output)
    double age_m = tm + 35.0;   // age at onset of metastatic cancer
    bool localised = (age_diag < age_m);
    double txhaz = (localised && (tx == RP || tx == RT)) ? in->parameter["RP_mortHR"] : 1.0; // assume same HR for RP & RT
    // calibration HR(age_diag,PSA,ext_grade) for loco-regional or HR(age_diag) for metastatic cancer
    double lead_time = age_c - age_diag;
    double txbenefit = exp(log(txhaz)+log(double(in->parameter["c_txlt_interaction"]))*lead_time); // treatment lead-time interaction
    double mort_hr = calculate_mortality_hr(age_diag);
    double ustar = pow(u,1/(in->parameter["c_baseline_specific"]*mort_hr*txbenefit*in->parameter["sxbenefit"]));
    if (localised)
      age_d = age_c + in->H_local[H_local_t::key_type(* in->H_local_age_set.lower_bound(bounds<double>(age_diag,50.0,80.0)),grade)].invert(-log(ustar));
    else
      age_d = age_c + in->H_dist[grade].invert(-log(ustar));
    if (in->debug) Rprintf("id=%i, lead_time=%f, ext_grade=%i, psamean=%f, tx=%i, txbenefit=%f, u=%f, ustar=%f, age_diag=%f, age_m=%f, age_c=%f, age_d=%f, mort_hr=%f\n",
		       id, lead_time, (int)ext_grade, psamean(age_diag), (int)tx, txbenefit, u, ustar, age_diag, age_m, age_c, age_d, mort_hr);
    return age_d;
  }

  double FhcrcPerson::onset() { return this->t0+35.0; }
  bool FhcrcPerson::onset_p() { return onset() <= now(); }

  /** @brief Calculate transition times for h(t) = y(t)*gamma = exp(beta0+beta1*t+beta2*(t-t0))*gamma
      This is equivalent to solving H(t) = gamma/(beta1+beta2)*(y(t)-y(s)) = -log(U) for entry time s.
  **/
  double FhcrcPerson::calculate_transition_time(double u, double t_enter, double gamma) {
    double y_enter = psamean(35.0 + t_enter);
    return (log(-log(u)*(beta1+beta2)/gamma + y_enter) - beta0 + beta2*t0) / (beta1+beta2);
  }

  void FhcrcPerson::opportunistic_rescreening(double psa) {
    double prescreened = 1.0 - in->rescreen_cure(bounds<double>(now(),30.0,90.0),psa);
    double shape = in->rescreen_shape(bounds<double>(now(),30.0,90.0),psa);
    double scale = in->rescreen_scale(bounds<double>(now(),30.0,90.0),psa);
    double u = R::runif(0.0,1.0);
    double t = now() + R::rweibull(shape,scale);
    if (u<prescreened) {
      scheduleAt(t, toScreen);
    }
  }

  void FhcrcPerson::cancel_events_after_diagnosis() {
    RemoveKind(toLocalised); // ?
    RemoveKind(toMetastatic);
    RemoveKind(toT3plus);
    RemoveKind(toScreen);
    RemoveKind(toOrganised);
    RemoveKind(toBiopsyFollowUpScreen);
    RemoveKind(toScreenInitiatedBiopsy);
    RemoveKind(toSTHLM3);
    RemoveKind(toOpportunistic);
    RemoveKind(toCancelScreens);
    RemoveKind(toScreenDiagnosis);
    RemoveKind(toOverDiagnosis);
    RemoveKind(toClinicalDiagnosis);
    RemoveKind(toClinicalDiagnosticBiopsy);
    RemoveKind(toMRI);
  }

  void FhcrcPerson::opportunistic_uptake_if_ever() {
    // Defaults values assume:
    // (i)   cohorts aged <35 in 1995 have a llogis(3.8,15) from age 35 (cohort > 1960)
    // (ii)  cohorts aged 50+ in 1995 have a llogis(2,10) distribution from 1995 (cohort < 1945)
    // (iii) intermediate cohorts are a weighted mixture of (i) and (ii)
    double first_screen;
    if (cohort > double(in->parameter["endUptakeMixture"])) {
      first_screen = in->parameter["uptakeStartAge"] + R::rllogis(in->parameter["shapeA"],
				       in->parameter["scaleA"]); // (i) age
    } else if (cohort < double(in->parameter["startUptakeMixture"])) {
      first_screen = (double(in->parameter["screeningIntroduced"]) - cohort) +
	R::rllogis(in->parameter["shapeT"],in->parameter["scaleT"]); // (ii) period
    } else {
      double age0 = double(in->parameter["screeningIntroduced"]) - cohort;
      double u = R::runif(0.0,1.0);
      if ((age0 - in->parameter["uptakeStartAge"]) / (double(in->parameter["endUptakeMixture"]) -
			   double(in->parameter["startUptakeMixture"])) < u) // (iii) mixture
	first_screen = age0 + R::rllogis_trunc(in->parameter["shapeA"],
					       in->parameter["scaleA"],
					       age0-in->parameter["uptakeStartAge"]);
      else first_screen = age0 + R::rllogis(in->parameter["shapeT"],
					    in->parameter["scaleT"]);
    }
    scheduleAt(first_screen, toScreen);
  }

  bool FhcrcPerson::screening_preference() {
    double pscreening = double(cohort>=in->parameter["startFullUptake"]) ?
      double(in->parameter["fullUptakePortion"]) : (double(in->parameter["fullUptakePortion"])
      - (double(in->parameter["startUptakeMixture"]) - cohort) * double(in->parameter["yearlyUptakeIncrease"]));
    // decrease for previous year instead of increase for next year
    double uscreening = R::runif(0.0,1.0);
    return (uscreening<pscreening);
  }

  void FhcrcPerson::rescreening_schedules(double psa, bool organised, bool mixed_programs) {
    // Check for organised screens - opportunistic screens are described later
    if (R::runif(0.0,1.0) < in->parameter["rescreeningParticipation"]) {
      switch (in->screen) {
      case mixed_screening:
      case stockholm3_goteborg:
      case goteborg:
        {
          if (organised && now() >= in->parameter["start_screening"] && now() < in->parameter["stop_screening"]) { // age groups
            if (psa<1.0 && now()+4.0 <= in->parameter["stop_screening"]) // re-screen late for low psa
              scheduleAt(now() + 4.0, toScreen);
            else if (psa>=1.0 && now()+2.0 <= in->parameter["stop_screening"]) // re-screen soon for moderate psa
              scheduleAt(now() + 2.0, toScreen);
            // else do nothing
          }
        }
        break;
      case introduced_screening: //rescreen
      case introduced_screening_preference:
      case introduced_screening_only:
        if (organised && now() + in->parameter["screening_interval"] <= in->parameter["stop_screening"]) {// within age?
          scheduleAt(now() + in->parameter["screening_interval"], toOrganised); //if there are planned opportunistic screens
        }
        break;
      case stockholm3_risk_stratified:
      case risk_stratified:
        if (now() >= in->parameter["start_screening"]) {
          if (psa < in->parameter["risk_psa_threshold"] &&
	      now()+in->parameter["risk_lower_interval"] <= in->parameter["stop_screening"])
            scheduleAt(now() + in->parameter["risk_lower_interval"], toScreen);
          if (psa >= in->parameter["risk_psa_threshold"] &&
	      now()+in->parameter["risk_higher_interval"] <= in->parameter["stop_screening"])
            scheduleAt(now() + in->parameter["risk_higher_interval"], toScreen);
        }
        break;
      case regular_screen:
        if (in->parameter["start_screening"] <= now() &&
            now() + in->parameter["screening_interval"] <= in->parameter["stop_screening"])
          scheduleAt(now() + in->parameter["screening_interval"], toScreen);
        break;
      case grs_stratified:
        if (now() + in->parameter["screening_interval"] <= in->parameter["stop_screening"])
          scheduleAt(now() + in->parameter["screening_interval"], toScreen);
        break;
      case grs_stratified_age:
	if (ageFirstScreen < in->parameter("screening_interval_split")) {
	  if (now() + in->parameter["screening_interval1"] <= in->parameter["stop_screening"])
	    scheduleAt(now() + in->parameter["screening_interval1"], toScreen);
	} else {
	  if (now() + in->parameter["screening_interval2"] <= in->parameter["stop_screening"])
	    scheduleAt(now() + in->parameter["screening_interval2"], toScreen);
	}
        break;
      case twoYearlyScreen50to70:
        if (50.0 <= now() && now() < 70.0)
          scheduleAt(now() + 2.0, toScreen);
        break;
      case sthlm3_mri_arm:
      case fourYearlyScreen50to70:
        if (50.0 <= now() && now() < 70.0)
          scheduleAt(now() + 4.0, toScreen);
        break;
      case screenUptake:
      case randomScreen50to70:
      case single_screen:
      case screen50:
      case screen60:
      case screen70:
      case stopped_screening:
      case cap_control:
      case cap_study:
	// case sthlm3_mri_arm:
        break;
      default:
        REprintf("Screening not matched: %s\n",in->screen);
        break;
      }
    } // rescreening participation
    if (in->screen == screenUptake ||
	in->screen == cap_control ||
	in->screen == cap_study ||
	// in->screen == sthlm3_mri_arm ||
	(mixed_programs && !organised))
      opportunistic_rescreening(psa); // includes rescreening participation
  } // rescreening

  bool FhcrcPerson::detectable(double now, double year) {
    return (state == Metastatic ||
        (state == Localised && ext_state == ext::T3plus) ||
        (state == Localised && ext_state == ext::T1_T2 &&
         (now > t3p + 35.0 - (t3p - t0) *
          in->parameter["biopsySensitivityTimeProportionT1T2"] *
          in->tableBiopsySensitivity(bounds(year,1987.0,2000.0)) /
          in->tableBiopsySensitivity(2000.0))));
  }

  double FhcrcPerson::callenderStartAge(double frailty, double threshold) {
    double tenYearRisk[] = {0.013, 0.015, 0.018, 0.02, 0.023, 0.026, 0.03, 0.033, 0.037,
			    0.041, 0.045, 0.049, 0.052, 0.056, 0.059, 0.062, 0.065,
			    0.067, 0.069, 0.071};
    for (int i=0; i<20; i++) {
      double rr = log(1.0-threshold)/log(1.0-tenYearRisk[i]);
      if (frailty>rr) return i+50;
    }
    return R_PosInf;
  }

  Double rbinorm(Double mean, Double sd, double rho) {
    double z1 = R::rnorm(0.0,1.0);
    double z2 = R::rnorm(0.0,1.0);
    z2 = rho*z1+sqrt(1-rho*rho)*z2;
    return Double(z1*sd.first+mean.first, mean.second+sd.second*z2);
  }
  Double rbinormPos(Double mean, Double sd, double rho, Double lbound = Double(0.0,0.0)) {
    Double val;
    do {
      val = rbinorm(mean,sd,rho);
    } while (val.first<lbound.first || val.second<lbound.second);
    return val;
  }
  RcppExport SEXP rbinorm_test() {
    RNGScope rng;
    Double x = rbinorm(Double(1.0,1.0), Double(2.0,2.0), 0.62);
    vector<double> v; v.push_back(x.first); v.push_back(x.second);
    return wrap(v);
  }
  RcppExport SEXP rbinormPos_test() {
    RNGScope rng;
    Double x = rbinormPos(Double(1.0,1.0), Double(2.0,2.0), 0.62);
    vector<double> v; v.push_back(x.first); v.push_back(x.second);
    return wrap(v);
  }
  double rweibull_frailty(double shape, double scale, double frailty=1.0) {
    return R::rweibull(shape, scale/std::pow(frailty,1.0/shape));
  }


/**
    Initialise a simulation run for an individual
 */
void FhcrcPerson::init() {

  // declarations
  double ym;

  // utilities
  utilities->clear();

  // background utilities
  // NumericVector bg_lower = in->background_utilities["lower"];
  // NumericVector bg_upper = in->background_utilities["upper"];
  // NumericVector bg_utility = in->background_utilities["utility"];
  for (int i=0; i < in->bg_utility.size(); i++)
	 scheduleUtilityChange(in->bg_lower[i], in->bg_upper[i], in->bg_utility[i]);

  // change state variables
  state = Healthy;
  ext_state = ext::Healthy_state;
  grade = base::Healthy;
  ext_grade = ext::Healthy;
  dx = NotDiagnosed;
  ageFirstScreen = -1.0;
  everPSA = previousNegativeBiopsy = organised = adt = previousFollowup = MRIpos = everGRS = false;
  in->rngNh->set();
  // https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.1002998&type=printable
  // genetic risk score with T ~ LogNormal(mu,sigma^2) with mean(T) =1 and sigma^2=0.68
  // E(T)=exp(mu+sigma^2/2)=1 (assumed mean)
  // => mu=-sigma^2/2 , sd=sqrt(sigma^2)
  // R: local({y=rlnorm(1e5,-1.68/2,sqrt(1.68)); c(mean(y), var(y))})
  // R: local({y=rlnorm(1e5,-1.68/2,sqrt(1.68)); plot(density(y,from=0),xlim=c(0,5))})
  if (in->bparameter["frailty"]) {
    double grs_log_mean = -in->parameter["grs_variance"]/2;
    double grs_log_sd = std::sqrt(in->parameter["grs_variance"]);
    double other_log_mean = -in->parameter["other_variance"]/2;
    double other_log_sd = std::sqrt(in->parameter["other_variance"]);
    grs_frailty = R::rlnorm(grs_log_mean, grs_log_sd);
    other_frailty = R::rlnorm(other_log_mean, other_log_sd);
  } else {
    grs_frailty = 1.0;
    other_frailty = 1.0;
  }
  if (R::runif(0.0, 1.0) <= in->parameter["susceptible"]) { // portion susceptible
    if (in->bparameter["weibull_onset"]) {
      t0 = rweibull_frailty(in->parameter["weibull_onset_shape"],
			    in->parameter["weibull_onset_scale"],
			    grs_frailty*other_frailty);
    } else {
      t0 = sqrt(2*R::rexp(1.0)/(in->parameter["g0"]*grs_frailty*other_frailty)); // is susceptible
    }
  }
  else
    t0 = 200.0; // not susceptible
  if (!in->bparameter["revised_natural_history"]){
    future_grade = (R::runif(0.0, 1.0)>=1+in->parameter["c_low_grade_slope"]*t0) ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(in->mubeta2[future_grade],in->sebeta2[future_grade]);
  }
  else {
    // multinomial logistic regression
    double u = R::runif(0.0,1.0);
    double denom = 1.0 +
      exp(in->parameter["alpha7"] + in->parameter["beta7"] * t0) +
      exp(in->parameter["alpha8"] + in->parameter["beta8"] * t0);
    double p6 = 1.0/denom;
    double p7 = exp(in->parameter["alpha7"] + in->parameter["beta7"] * t0) / denom;
    // double p8 = exp(in->parameter["alpha8"] + in->parameter["beta8"] * t0) / denom;
    if (u < p6) future_ext_grade = ext::Gleason_le_6;
    else if (u < p6+p7) future_ext_grade = ext::Gleason_7;
    else future_ext_grade = ext::Gleason_ge_8;
    future_grade = future_ext_grade == ext::Gleason_ge_8 ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(in->mubeta2[future_ext_grade],in->sebeta2[future_ext_grade]);
  }
  beta0 = R::rnorm(in->parameter["mubeta0"],in->parameter["sebeta0"]);
  beta1 = R::rnormPos(in->parameter["mubeta1"],in->parameter["sebeta1"]);

  y0 = psamean(t0+35); // depends on: t0, beta0, beta1, beta2
  t3p = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["g3p"]);
  tm = calculate_transition_time(R::runif(0.0,1.0), t3p, in->parameter["gm"]);
  ym = psamean(tm+35);
  if (future_grade==base::Gleason_le_7) { // Annals
    tc = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["gc"]);
    tmc = calculate_transition_time(R::runif(0.0,1.0), tm, in->parameter["gc"]*in->parameter["thetac"]);
  } else {
    tc = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["gc"]*in->parameter["grade.clinical.rate.high"]);
    tmc = calculate_transition_time(R::runif(0.0,1.0), tm, in->parameter["gc"]*in->parameter["thetac"]*in->parameter["grade.clinical.rate.high"]);
  }
  out->tmc_minus_t0 += (tmc - t0);
  aoc = in->rmu0.rand(R::runif(0.0,1.0));
  if (!in->bparameter["revised_natural_history"]){
    future_ext_grade= (future_grade==base::Gleason_le_7) ?
      (R::runif(0.0,1.0) <= in->interp_prob_grade7.approx(beta2) ? ext::Gleason_7 : ext::Gleason_le_6) :
      ext::Gleason_ge_8;
  }
  ageEntry = 0.0; // used by cap_control, cap_study and sthlm3_mri_arm

  if (in->debug) {
    Rprintf("id=%i, future_grade=%i, future_ext_grade=%i, beta0=%f, beta1=%f, beta2=%f, mubeta0=%f, sebeta0=%f, mubeta1=%f, sebeta1=%f, mubeta2=%f, sebeta2=%f\n", id, future_grade, future_ext_grade, beta0, beta1, beta2, double(in->parameter["mubeta0"]), double(in->parameter["sebeta0"]), double(in->parameter["mubeta1"]), double(in->parameter["sebeta1"]), in->mubeta2[future_grade], in->sebeta2[future_grade]);
  }

  tx = no_treatment;
  txhaz = -1.0;
  // schedule natural history events
  scheduleAt(t0+35.0,toLocalised);
  scheduleAt(aoc,toOtherDeath);

  // schedule screening events that depend on screeningParticipation
  in->rngScreen->set();
  rescreening_frailty = R::rgamma(1.25, 1.25);
  double u1 = R::runif(0.0,1.0);
  double u2 = R::runif(0.0,1.0);
  if (R::runif(0.0,1.0)<in->parameter["screeningParticipation"]) {
    switch(in->screen) {
    case noScreening:
      break; // no screening
    case randomScreen50to70:
      scheduleAt(50.0 + u1 * 20.0, toScreen);
      break;
    case single_screen:
    case regular_screen:
    case goteborg:
    case risk_stratified:
      scheduleAt(in->parameter["start_screening"],toScreen);
      break;
    case fourYearlyScreen50to70: // 50,54,58,62,66,70
    case twoYearlyScreen50to70:  // 50,52, ..., 68,70
    case screen50:
      scheduleAt(50.0,toScreen);
      break;
    case screen60:
      scheduleAt(60.0,toScreen);
      break;
    case screen70:
      scheduleAt(70.0,toScreen);
      break;
    case grs_stratified_age:
    case grs_stratified: {
      double startAge = callenderStartAge(grs_frailty, in->parameter("grs_risk_threshold"));
      if (startAge != R_PosInf)
	scheduleAt(startAge,toScreen);
    } break;
    case mixed_screening:
    case introduced_screening:
    case introduced_screening_preference:
    case introduced_screening_only:
    case stopped_screening:
    case stockholm3_goteborg:
    case stockholm3_risk_stratified:
    case screenUptake:
    case cap_control:
    case cap_study:
    case sthlm3_mri_arm:
      // see below (models screening participation)
      break;
    default:
      REprintf("Screening not matched: %i\n",in->screen);
      break;
    }
  }

  // schedule screening events that already incorporate screening participation
  switch(in->screen) {
  case mixed_screening:
    if (screening_preference())
      opportunistic_uptake_if_ever();
    scheduleAt(in->parameter["start_screening"], toOrganised);
    break;
  case stopped_screening:
    if (screening_preference())
      opportunistic_uptake_if_ever();
    scheduleAt(in->parameter["introduction_year"] - cohort, toCancelScreens);
    break;
  case introduced_screening: //first screen
    // One participation during opportunistic and another during the regular screening
    if (screening_preference())
      opportunistic_uptake_if_ever(); // 'toOrganised' will remove opportunistic screens
    if ( in->parameter["introduction_year"] - cohort <= in->parameter["start_screening"]) { // under screen age at 2015
      scheduleAt(in->parameter["start_screening"], toOrganised); //screen all in age interval
    } else if( in->parameter["introduction_year"] - cohort >= in->parameter["start_screening"] && //between screen ages
               in->parameter["introduction_year"] - cohort <= in->parameter["stop_screening"]) {
      scheduleAt(u1 + in->parameter["introduction_year"] - cohort, toOrganised); //in 1 year screen all in age interval
    }
    break;
  case introduced_screening_preference: //first screen
    // Only those who would have had a opportunistic screen will have a regular screen
    if (screening_preference()) {
      opportunistic_uptake_if_ever(); // 'toOrganised' will remove opportunistic screens
      if ( in->parameter["introduction_year"] - cohort <= in->parameter["start_screening"]) { // under screen age at 2015
	scheduleAt(in->parameter["start_screening"], toOrganised); //screen all in age interval
      } else if( in->parameter["introduction_year"] - cohort >= in->parameter["start_screening"] && //between screen ages
		 in->parameter["introduction_year"] - cohort <= in->parameter["stop_screening"]) {
	scheduleAt(u1 + in->parameter["introduction_year"] - cohort, toOrganised); //in 1 year screen all in age interval
      }
    }
    break;
  case introduced_screening_only: //first screen
    if ( in->parameter["introduction_year"] - cohort <= in->parameter["start_screening"]) {
      scheduleAt(in->parameter["start_screening"], toOrganised); //screen all in age interval
    } else if( in->parameter["introduction_year"] - cohort >= in->parameter["start_screening"] && //between screen ages
	       in->parameter["introduction_year"] - cohort <= in->parameter["stop_screening"]) {
      scheduleAt(u1 + in->parameter["introduction_year"] - cohort, toOrganised); //in 1 year screen all in age interval
    }
    break;
  case stockholm3_goteborg:
  case stockholm3_risk_stratified:
    if (screening_preference())
      opportunistic_uptake_if_ever();
    if (u1 < in->parameter["studyParticipation"] &&
	(2013.0-cohort >= in->parameter["start_screening"] &&
	 2013.0-cohort < in->parameter["stop_screening"]))
      scheduleAt((u2 * 2.0 + 2013.0) - cohort, toSTHLM3);
    break;
  case screenUptake:
    if (screening_preference())
      opportunistic_uptake_if_ever();
    break;
  case cap_control:
  case cap_study: {
    scheduleAt(in->H_screen_uptake.invert(-log(R::runif(0.0,1.0))), toScreen);
    ageEntry = 2005.0 - cohort + R::runif(0.0,5.0); // 50-54, 55-59, 60-64, 65-69
    if (in->screen == cap_study) {
      double pScreened = in->cap_pScreened[cohort==1955.0 ? 0 :
					    cohort==1950.0 ? 1 :
					    cohort==1945.0 ? 2 :
					    3];
      if (R::runif(0.0,1.0)<pScreened)
	scheduleAt(ageEntry, toOrganised);
    } else {
      // keep the RNG in sync
      R::runif(0.0,1.0);
    }
  } break;
  case sthlm3_mri_arm:
    if (screening_preference())
      opportunistic_uptake_if_ever();
    ageEntry = R::runif(2019.0,2020.0) - cohort;
    scheduleAt(ageEntry, toOrganised);
    break;
  default:
    break;
  }

  in->rngNh->set();

  // record some parameters using SimpleReport - too many for a tuple
  if (id < in->nLifeHistories) {
    out->outParameters.record("id",double(id));
    out->outParameters.record("beta0",beta0);
    out->outParameters.record("beta1",beta1);
    out->outParameters.record("beta2",beta2);
    out->outParameters.record("t0",t0);
    out->outParameters.record("t3p",t0);
    out->outParameters.record("tm",tm);
    out->outParameters.record("tc",tc);
    out->outParameters.record("tmc",tmc);
    out->outParameters.record("y0",y0);
    out->outParameters.record("ym",ym);
    out->outParameters.record("aoc",aoc);
    out->outParameters.record("cohort",cohort);
    out->outParameters.record("future_ext_grade",future_ext_grade);
    out->outParameters.record("ext_grade",ext_grade);
    out->outParameters.record("age_psa",-1.0);
    out->outParameters.record("age_pca",-1.0);
    out->outParameters.record("pca_death",0.0);
    out->outParameters.record("psa55",psameasured(55.0));
    out->outParameters.record("psa65",psameasured(65.0));
    out->outParameters.record("psa75",psameasured(75.0));
    out->outParameters.record("psa85",psameasured(85.0));
    out->outParameters.record("rescreening_frailty",rescreening_frailty);
    out->outParameters.record("ageEntry",ageEntry);
  }

  if (in->debug) Rprint_actions();

}

/**
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // by default, use the natural history RNG
  in->rngNh->set();

  // declarations
  in->rngOther->set();
  double psa = psameasured(now()); // includes measurement error
  in->rngNh->set();
  // double test = panel ? biomarker : psa;
  double Z = psamean(now());
  double age = now();
  double year = age + cohort;
  double compliance;
  bool mixed_programs = (in->screen == mixed_screening) ||
    (in->screen == introduced_screening) ||
    (in->screen == introduced_screening_preference) ||
    (in->screen == stopped_screening) ||
    (in->screen == cap_study);
  bool formal_costs = in->parameter["formal_costs"]==1.0 && (!mixed_programs || organised);
  bool formal_compliance = in->parameter["formal_compliance"]==1.0 && (!mixed_programs || organised);
  double utility = FhcrcPerson::utility();
  in->rngPrelude->set();
  bool detectable = FhcrcPerson::detectable(now(), year);
  if (in->parameter["rand_biopsy_sensitivityG6"]<1.0) {
    detectable = detectable && R::runif(0.0,1.0) < in->parameter["rand_biopsy_sensitivityG6"];
  }
  in->rngNh->set();

  // record information
  switch(int(in->parameter["full_report"])) {
  case 1:
    out->report.add(FullState::Type(ext_state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, age, utility, index);
    break;
  case 0:
    out->shortReport.add(1, msg->kind, previousEventTime, age, utility);
    break;
  default:
    out->shortReport.addBrief(previousEventTime, age, utility);
    break;
  }

  if (in->bparameter["includeEventHistories"] && id < in->nLifeHistories) { // only record up to the first n individuals
    out->lifeHistories.push_back(LifeHistory::Type(id, ext_state, ext_grade, dx, msg->kind, previousEventTime, age, year, psa, utility));
  }

  if (in->debug)
    Rprint(*utilities);

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath:
    if (in->bparameter["Andreas"]) {
      lost_productivity("Terminal illness");
      add_costs("Cancer death");
    }
    if (id < in->nLifeHistories) {
      out->outParameters.record("age_d",now());
      out->outParameters.revise("pca_death",1.0);
    }
    out->report.individualReset();
    out->shortReport.individualReset();
    out->costs.individualReset();
    Sim::stop_simulation();
    break;

  case toOtherDeath:
    // add_costs("Death"); // cost for death, should this be zero???

    if (id < in->nLifeHistories) {
      out->outParameters.record("age_d",now());
    }
    out->report.individualReset();
    out->shortReport.individualReset();
    out->costs.individualReset();
    Sim::stop_simulation();
    break;

  case toLocalised:
    state = Localised; ext_state = ext::T1_T2;
    ext_grade = future_ext_grade;
    grade = future_grade;
    if (now()<tc+35.0-6.0/52.0)
      scheduleAt(tc+35.0-6.0/52.0,toClinicalDiagnosticBiopsy);
    if (now()<tc+35.0-3.0/52.0)
      scheduleAt(tc+35.0-3.0/52.0,toClinicalDiagnosticBiopsy);
    scheduleAt(tc+35.0,toClinicalDiagnosis);
    scheduleAt(t3p+35.0,toT3plus);
    scheduleAt(tm+35.0,toMetastatic);
    break;

  case toT3plus:
    ext_state = ext::T3plus;
    break;

  case toMetastatic:
    state = Metastatic; ext_state = ext::Metastatic;
    RemoveKind(toClinicalDiagnosis);
    RemoveKind(toClinicalDiagnosticBiopsy);
    if (in->bparameter["Andreas"]) {
      if (now()<tc+35.0-6.0/52.0) // should this be tmc?
	scheduleAt(tmc+35.0-6.0/52.0,toClinicalDiagnosticBiopsy);
      if (now()<tc+35.0-3.0/52.0) // should this be tmc?
	scheduleAt(tmc+35.0-3.0/52.0,toClinicalDiagnosticBiopsy);
    } else {
      if (now()<tmc+35.0-6.0/52.0)
	scheduleAt(tmc+35.0-6.0/52.0,toClinicalDiagnosticBiopsy);
      if (now()<tmc+35.0-3.0/52.0)
	scheduleAt(tmc+35.0-3.0/52.0,toClinicalDiagnosticBiopsy);
    }
    scheduleAt(tmc+35.0,toClinicalDiagnosis);
    // Remove possible secondary Tx
    RemoveKind(toRP);
    RemoveKind(toRT);
    break;

  case toOrganised:
  case toSTHLM3:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    scheduleAt(now(), toScreen); // now start organised screening
    break;

  case toCancelScreens:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    break;

  case toScreen:
  case toBiopsyFollowUpScreen: {
    if (ageFirstScreen < 0.0) ageFirstScreen = now();
    in->rngBx->set();
    this->psa_last_screen = psa;
    if (in->bparameter["includePSArecords"]) {
      out->psarecord.record("id",id);
      out->psarecord.record("state",state);
      out->psarecord.record("ext_grade",ext_grade);
      out->psarecord.record("organised",organised); // only meaningful for mixed_programs
      out->psarecord.record("dx",dx);
      out->psarecord.record("age",age);
      out->psarecord.record("cohort",cohort);
      out->psarecord.record("psa",psa_last_screen);
      out->psarecord.record("t0",t0);
      out->psarecord.record("beta0",beta0);
      out->psarecord.record("beta1",beta1);
      out->psarecord.record("beta2",beta2);
      out->psarecord.record("Z",Z);
      out->psarecord.record("onset",double(onset_p()));
      out->psarecord.record("detectable",double(detectable));
    }
    if ((in->screen==grs_stratified || in->screen==grs_stratified_age) && !everGRS) {
      add_costs("Polygenic risk stratification"); // once-only cost
      everGRS = true;
    }
    if (!everPSA) {
      if (id < in->nLifeHistories) {
	out->outParameters.revise("age_psa",now());
	// outParameters.revise("first_psa",psa);
      }
      everPSA = true;
    }
    if (formal_costs) {
      add_costs("Invitation");
      lost_productivity(in->panel && psa>=in->parameter["panelReflexThreshold"] ? "Formal panel" : "Formal PSA");
      add_costs(in->panel && psa>=in->parameter["panelReflexThreshold"] ? "Formal panel" : "Formal PSA");
      scheduleUtilityChange(now(), "Formal PSA");
    } else { // opportunistic costs
      add_costs(in->panel && psa>=in->parameter["panelReflexThreshold"] ? "Opportunistic panel" : "Opportunistic PSA");
      lost_productivity(in->panel && psa>=in->parameter["panelReflexThreshold"] ? "Opportunistic panel" : "Opportunistic PSA");
      scheduleUtilityChange(now(), "Opportunistic PSA");
    }
    compliance = formal_compliance ?
      in->tableFormalBiopsyCompliance(bounds<double>(psa,3.0,10.0),
				  bounds<double>(age,40,80)) :
      in->tableOpportunisticBiopsyCompliance(bounds<double>(psa,3.0,10.0),
					 bounds<double>(age,40,80));
    bool positive_test =
      (msg->kind == toScreen && psa >= in->parameter["psaThreshold"]) ? true :
      (msg->kind == toBiopsyFollowUpScreen && psa >= in->parameter["psaThresholdBiopsyFollowUp"]) ? true :
      false;
    // Important case: PSA<1 (to check)
    // Reduce false positives wrt Gleason 7+ by 1-rFPF: which BPThreshold?
    if (in->panel && positive_test && (!in->bparameter["Andreas"] || psa < 10.)) {
      if (int(in->parameter("biomarker_model"))==random_correction) { // simplistic model for the biomarker
	if (R::runif(0.0,1.0) < 1.0 - in->parameter["rFPF"])
	  positive_test = false;
      }
      else if (int(in->parameter("biomarker_model"))==psa_informed_correction) { // PSA based model for the biomarker
	if ((ext_grade == ext::Gleason_le_6 &&
	     detectable && psa < in->parameter["PSA_FP_threshold_GG6"]) // FP GG 6 PSA threshold
	    ||  (!detectable && psa < in->parameter["PSA_FP_threshold_nCa"]) // FP no cancer PSA threshold
	    || ((ext_grade == ext::Gleason_7 || ext_grade == ext::Gleason_ge_8) &&
		detectable && psa < in->parameter["PSA_FP_threshold_GG7plus"])) { // FP GG >= 7 PSA threshold
	  positive_test = false; // assumption relying on PSA being a strong panel component
          if (in->debug) Rprintf("Panel adjusted tests id=%i, psa=%8.6f, ext_grade=%i, future_ext_grade=%i, onset=%d, detectable=%d\n", id, psa, ext_grade, future_ext_grade, onset_p(), detectable);
	}
      }
      else {
	REprintf("Parameter biomarker_model not matched: %i\n", int(in->parameter("biomarker_model")));
      }
    }
    // Case: PSA>=10. The man has a positive_test.
    if (in->bparameter["includePSArecords"] && !onset_p() && positive_test) {
      out->falsePositives.record("id",id);
      out->falsePositives.record("psa",psa);
      out->falsePositives.record("age",now());
      out->falsePositives.record("age0",t0+35.0);
      out->falsePositives.record("ext_grade",ext_grade);
    }
    // if (panel && !positive_test && t0<now()-35.0 && ext_grade > ext::Gleason_le_6) {
    //   if (R::runif(0.0,1.0) < 1.0-parameter["rTPF"]) positive_test = true;
    // }
    in->rngBx->set();
    if (positive_test && R::runif(0.0,1.0) < compliance) {
      if (in->bparameter["MRI_screen"]) {
	scheduleAt(now()+1.0/52.0, toMRI); // MRI in one month
      } else {
	scheduleAt(now()+1.0/52.0, toScreenInitiatedBiopsy); // biopsy in one month
      }
    } // assumes similar biopsy compliance, reasonable? An option to different psa-thresholds would be to use different biopsyCompliance. /AK
    else {
          in->rngScreen->set();
	  if ((in->screen == cap_study || in->screen == sthlm3_mri_arm) && organised)
	    organised = false;
	  rescreening_schedules(psa, organised, mixed_programs);
    }
    in->rngNh->set();
  } break;

  case toClinicalDiagnosis:
    scheduleUtilityChange(now(), "Cancer diagnosis");
    dx = ClinicalDiagnosis;
    cancel_events_after_diagnosis();
    scheduleAt(now()+1.0/12.0, toTreatment);
    if (id < in->nLifeHistories) {
      out->outParameters.revise("age_pca",now());
    }
    break;

  case toScreenDiagnosis:
    scheduleUtilityChange(now(), "Cancer diagnosis");
    // add cost for half a subsequent consultation
    add_costs("Assessment", Direct, 0.5);
    dx = ScreenDiagnosis;
    cancel_events_after_diagnosis();
    scheduleAt(now()+1.0/12.0, toTreatment); // treatment one month after the diagnosis
    if (aoc < (35.0 + min(tc,tmc))) {
      scheduleAt(now(), toOverDiagnosis);
    }
    if (id < in->nLifeHistories) {
      out->outParameters.revise("age_pca",now());
    }
    break;

  case toOverDiagnosis:
    // only for recording
    break;

  case toMRI: {
    in->rngBx->set();
    add_costs("MRI"); // does this include costs for the consultation?
    lost_productivity("MRI");
    // scheduleUtilityChange(now(), "MRI");
    double pMRIpos = (this->ext_grade == ext::Healthy || !detectable) ? in->parameter["pMRIposG0"] :
      (this->ext_grade == ext::Gleason_le_6) ? in->parameter["pMRIposG1"] :
      in->parameter["pMRIposG2"];
    this->MRIpos = (R::runif(0.0,1.0) < pMRIpos); // we need to know if they are MRI+ at toScreenInitiatedBiopsy
    if (this->MRIpos || in->bparameter("MRInegSBx")) {
      scheduleAt(now(), toScreenInitiatedBiopsy);
    } else {
      if ((in->screen == cap_study || in->screen == sthlm3_mri_arm) && organised)
	organised = false;
      rescreening_schedules(psa, organised, mixed_programs);
    }
    in->rngNh->set();
  } break;
    
  // record additional biopsies for clinical diagnoses
  case toClinicalDiagnosticBiopsy:
    if (in->bparameter["MRI_clinical"]) {
      add_costs("Combined biopsy");
      lost_productivity("Combined biopsy");
    } else {
      add_costs("Biopsy");
      lost_productivity("Biopsy");
    }
    // common values
    add_costs("Assessment");
    lost_productivity("Assessment");
    scheduleUtilityChange(now(), "Biopsy");
    break;

  case toScreenInitiatedBiopsy: {
    in->rngBx->set();
    double u1 = R::runif(0.0,1.0);
    double u2 = R::runif(0.0,1.0);
    // the general case for the following block follows the same pattern as toClinicalDiagnosticBiopsy
    if (in->bparameter["MRI_screen"] && this->MRIpos) {
      // specific case: MRI+S3M+ or MRI-S3M>=25 (split by type of biopsy: combined and systematic, respectively)
      if (in->bparameter.containsElementNamed("SplitS3M25plus") &&
	  in->bparameter["SplitS3M25plus"]) {
	double pS3Mpos = (detectable && this->ext_grade == ext::Gleason_le_6) ? in->parameter["PrS3MposIfBx_GG6"] :
	  (detectable && this->ext_grade == ext::Gleason_7) ? in->parameter["PrS3MposIfBx_GG7plus"] :
	  (detectable && this->ext_grade == ext::Gleason_ge_8) ? in->parameter["PrS3MposIfBx_GG7plus"] :
	  in->parameter["PrS3MposIfBx_nCa"];
	if (u1 < pS3Mpos) {
	  add_costs("Combined biopsy");
	  lost_productivity("Combined biopsy");
	} else {
	  add_costs("Biopsy");
	  lost_productivity("Biopsy");
	}
      }
      else { // general case
	add_costs("Combined biopsy");
	lost_productivity("Combined biopsy");
      }
    } else {
      add_costs("Biopsy");
      lost_productivity("Biopsy");
    }
    // common values
    add_costs("Assessment");
    lost_productivity("Assessment");
    scheduleUtilityChange(now(), "Biopsy");

    // output biopsy record
    if (in->bparameter["includeBxrecords"]) {
      out->bxrecord.record("id",id);
      out->bxrecord.record("state",state);
      out->bxrecord.record("ext_state",ext_state);
      out->bxrecord.record("ext_grade",ext_grade);
      out->bxrecord.record("organised",organised); // only meaningful for mixed_programs
      out->bxrecord.record("dx",dx);
      out->bxrecord.record("age",age);
      out->bxrecord.record("cohort",cohort);
      out->bxrecord.record("psa",psa_last_screen);
      out->bxrecord.record("t0",t0);
      out->bxrecord.record("beta0",beta0);
      out->bxrecord.record("beta1",beta1);
      out->bxrecord.record("beta2",beta2);
      out->bxrecord.record("Z",Z);
      out->bxrecord.record("onset",double(onset_p()));
      out->bxrecord.record("detectable",double(detectable));
    }

    bool Bx_missed = false; // allow for randomly missing "detectable cancers" from Bx
    if (detectable) {
      if (in->bparameter["Andreas"]) {
	// pass
      } else if (in->bparameter["MRI_screen"] && this->MRIpos) {
	if (this->ext_grade == ext::Gleason_le_6) {
	  Bx_missed = (u2 < in->parameter["pTBxG0ifG1_MRIpos"]);
	}
	if (this->ext_grade == ext::Gleason_7) {
	  Bx_missed = (u2 < in->parameter["pTBxG0ifG2_MRIpos"]);
	}
      } else { // SBx compared with MRI
	if (this->ext_grade == ext::Gleason_le_6) {
	  Bx_missed = (u2 < in->parameter["pSBxG0ifG1"]);
	}
	if (this->ext_grade == ext::Gleason_7) {
	  Bx_missed = (u2 < in->parameter["pSBxG0ifG2"]);
	}
      }
      if (!Bx_missed)
	scheduleAt(now()+3.0/52.0, toScreenDiagnosis); // diagnosis three weeks after biopsy
      if (in->panel && state==Localised && ext_grade == ext::Gleason_le_6) {
	// fixed costs etc for men who were S3M+/PE-
	add_costs("Assessment", Direct, 766.0/722.0 - 1.0);
	lost_productivity("Assessment", 766.0/722.0 - 1.0);
      }
    }
    if (!detectable || Bx_missed) { // negative biopsy
      if (in->panel) { // fixed costs etc for men who were S3M+/PE-
        add_costs("Assessment", Direct, 1535.0/1381.0 - 1.0);
        lost_productivity("Assessment", 1535.0/1381.0 - 1.0);
      }
      if (previousNegativeBiopsy || (in->bparameter("MRI_screen") && !this->MRIpos && in->bparameter["rescreenDoubleNeg"])) { // first negative biopsy and not MRI-
	if ((in->screen == cap_study || in->screen == sthlm3_mri_arm) && organised)
	  organised = false;
        rescreening_schedules(psa_last_screen, organised, mixed_programs);
        previousNegativeBiopsy = false; // if going to rescreening, should their previous negative Bx be forgotten? (Currently only used here)
      } else {
        previousNegativeBiopsy=true;
        // Competing risk for event following a negative biopsy
        double timeToPSA = R::rlnorm(in->tableNegBiopsyToPSAmeanlog(age),
                                     in->tableNegBiopsyToPSAsdlog(age));
        double timeToBiopsy = R::rlnorm(in->tableNegBiopsyToBiopsymeanlog(age),
                                        in->tableNegBiopsyToBiopsysdlog(age));
        if (timeToPSA <= timeToBiopsy) { // PSA was the first event
          scheduleAt(now() + timeToPSA, toScreen);
        } else { // Biopsy was the first event
	  this->MRIpos = false; // HACK!!! This ensures that they do SBx only.
          scheduleAt(now() + timeToBiopsy, toScreenInitiatedBiopsy);
        }
      }
    }
    in->rngNh->set();
  } break;

  case toTreatment: { // To diagnoses, treatment & survival
    in->rngTreatment->set();
    double u_tx = R::runif(0.0,1.0);
    double u_adt = R::runif(0.0,1.0);
    if (state == Metastatic && in->bparameter["Andreas"]) {
      lost_productivity("Metastatic cancer");
      // utilities->clear(); // should this be age-specific??
    }
    else { // Loco-regional
      if (!in->bparameter["Andreas"] && now() < 65.0)
	lost_productivity("Long-term sick leave");
      tx = calculate_treatment(u_tx,now(),year);
      if (tx == CM) scheduleAt(now(), toCM);
      if (tx == RP) scheduleAt(now(), toRP);
      if (tx == RT) scheduleAt(now(), toRT);
      if (in->bparameter["Andreas"]) {
	// check for ADT
	double pADT =
	  in->pradt(tx,
		    bounds<double>(now(),50,79),
		    bounds<double>(year,1973,2004),
		    grade);
	if (u_adt < pADT)  {
	  adt = true;
	  scheduleAt(now(), toADT);
	}
	if (in->debug) Rprintf("id=%i, adt=%d, u=%8.6f, pADT=%8.6f\n",id,adt,u_adt,pADT);
      }
    }
    // reset the random number stream
    in->rngSurv->set();
    // check for cure
    bool cured = false;
    double age_c = (state == Localised) ? tc + 35.0 : tmc + 35.0;
    double lead_time = age_c - now();
    // calculate the age at cancer death by c_benefit_type
    double age_cancer_death=R_PosInf;
    double age_cd = R_PosInf, age_sd = R_PosInf, weight = R_PosInf;
    if (in->parameter["c_benefit_type"]==LeadTimeBased) { // [new paper ref]
      double pcure = pow(1 - exp(-lead_time * in->parameter["c_benefit_value1"]),
      			 calculate_mortality_hr(age_c));
      if (in->debug) Rprintf("hr for lead time=%f\n", calculate_mortality_hr(age_c));
      cured = (R::runif(0.0,1.0) < pcure);
      if (!cured) {
	double u_surv = R::runif(0.0,1.0);
        age_cancer_death = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      }
    }
    else if (in->parameter["c_benefit_type"]==StageShiftBased) { // [annals paper ref]
      // calculate survival
      double u_surv = R::runif(0.0,1.0);
      age_cd = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      age_sd = calculate_survival(u_surv,now(),age_c,tx);
      weight = exp(- in->parameter["c_benefit_value0"]*lead_time);
      age_cancer_death = weight*age_cd + (1.0-weight)*age_sd;
    }
    else REprintf("c_benefit_type not matched.");
    if (!cured) {
      scheduleAt(age_cancer_death, toCancerDeath);
      // Disutilities prior to a cancer death
      if (in->bparameter["Andreas"]) {
	double age_palliative = age_cancer_death - in->utility_duration["Palliative therapy"] - in->utility_duration["Terminal illness"];
	double age_terminal = age_cancer_death - in->utility_duration["Terminal illness"];
	if (age_palliative>now()) { // cancer death more than 36 months after diagnosis
	  scheduleUtilityChange(age_palliative, age_terminal,
				in->utility_estimates["Palliative therapy"]);
	  scheduleUtilityChange(age_terminal, "Terminal illness");
	}
	else if (age_terminal>now()) { // cancer death between 36 and 6 months of diagnosis
	  scheduleUtilityChange(now(), age_terminal, in->utility_estimates["Palliative therapy"]);
	  scheduleUtilityChange(age_terminal,"Terminal illness");
	}
	else // cancer death within 6 months of diagnosis/treatment
	  scheduleUtilityChange(now(), "Terminal illness");
      } else { // Shuang:
	double age_adt = age_cancer_death - in->utility_duration["Palliative therapy"] - in->utility_duration["Terminal illness"]
	   - in->utility_duration["ADT+chemo"];
	double age_palliative = age_cancer_death - in->utility_duration["Palliative therapy"] - in->utility_duration["Terminal illness"];
	double age_terminal = age_cancer_death - in->utility_duration["Terminal illness"];
	// check age_adt (two cases)
	if (now() <= age_adt) { 
	  scheduleUtilityChange(age_adt, age_palliative,
				in->utility_estimates["ADT+chemo"]);
	  scheduleAt(age_adt, toADT);
	} else if (now() < age_palliative) {
	  scheduleUtilityChange(now(), age_palliative,
				in->utility_estimates["ADT+chemo"]);
	  scheduleAt(now(), toADT);
	}
	// check age_palliative (two cases)
	if (now() <= age_palliative) { 
	  scheduleUtilityChange(age_palliative, age_terminal,
				in->utility_estimates["Palliative therapy"]);
	  scheduleAt(age_palliative, toPalliative);
	} else if (now() < age_palliative) {
	  scheduleUtilityChange(now(), age_palliative,
				in->utility_estimates["Palliative therapy"]);
	  scheduleAt(now(), toPalliative);
	}
	// check age_terminal (two cases)
	if (now() <= age_terminal) { 
	  scheduleUtilityChange(age_terminal, "Terminal illness");
	  scheduleAt(age_terminal, toTerminal);
	} else if (now() < age_terminal) {
	  scheduleUtilityChange(now(), "Terminal illness");
	  scheduleAt(now(), toTerminal);
	}
      }
    }
    if (in->bparameter["includeDiagnoses"]) {
      out->diagnoses.record("id",id);
      out->diagnoses.record("age",age);
      out->diagnoses.record("year",year);
      out->diagnoses.record("psa",psa_last_screen);
      out->diagnoses.record("ext_grade",ext_grade);
      out->diagnoses.record("ext_state",ext_state);
      out->diagnoses.record("organised",organised); // only meaningful for mixed_program, keep this?
      out->diagnoses.record("dx",dx);
      out->diagnoses.record("tx",tx);
      out->diagnoses.record("cancer_death",(aoc>age_cancer_death) ? 1.0 : 0.0);
      out->diagnoses.record("age_at_death", (aoc>age_cancer_death) ? age_cancer_death : aoc);
      out->diagnoses.record("age_cancer_death", age_cancer_death);
      out->diagnoses.record("aoc", aoc);
      out->diagnoses.record("age_cd", age_cd);
      out->diagnoses.record("age_sd", age_sd);
      out->diagnoses.record("weight", weight);
      out->diagnoses.record("lead_time", lead_time);
    }
    in->rngNh->set();
  } break;

  case toRP:
    add_costs("Prostatectomy");
    this->previousFollowup = false;
    scheduleAt(now() + 1.0, toYearlyPostTxFollowUp);
    lost_productivity("Prostatectomy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "Prostatectomy part 1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + in->utility_duration["Prostatectomy part 1"],
			  "Prostatectomy part 2");
    scheduleUtilityChange(now() + in->utility_duration["Prostatectomy part 1"] +
                          in->utility_duration["Prostatectomy part 2"], "Postrecovery period");
    // Remove yearly active surveillance if the RP is the secondary Tx
    RemoveKind(toYearlyActiveSurveillance); // breaks recursive call
    RemoveKind(toRT);
    break;

  case toRT:
    add_costs("Radiation therapy");
    this->previousFollowup = false;
    scheduleAt(now() + 1.0, toYearlyPostTxFollowUp);
    lost_productivity("Radiation therapy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "Radiation therapy part 1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + in->utility_duration["Radiation therapy part 1"],
			  "Radiation therapy part 2");
    scheduleUtilityChange(now() + in->utility_duration["Radiation therapy part 1"] +
                          in->utility_duration["Radiation therapy part 2"], "Postrecovery period");
    RemoveKind(toYearlyActiveSurveillance); // breaks recursive call
    break;

  case toCM:
    in->rngTreatment->set();
    if (in->bparameter["Andreas"])
      add_costs("Active surveillance - single MR"); // expand here
    scheduleAt(now(), toYearlyActiveSurveillance);
    scheduleUtilityChange(now(), "Active surveillance");
    // Modelling for possible subsequent RP and RT. P(RP|RT) ~ P(RP)
    // whereas P(RT|RP) << P(RT). As a simplification, we simulate
    // separately for RP and RT and remove an RT following an RP.
    if (R::runif(0.0,1.0) > in->tableCMtoRPpnever(age)) {// pnever -> pever
      scheduleAt(now() + R::rlnorm(in->tableCMtoRPmeanlog(age), in->tableCMtoRPsdlog(age)), toRP);
    }
    if (R::runif(0.0,1.0) > in->tableCMtoRTpnever(age)) {// pnever -> pever
      scheduleAt(now() + R::rlnorm(in->tableCMtoRTmeanlog(age), in->tableCMtoRTsdlog(age)), toRT);
    }
    in->rngNh->set();
    break;

  case toYearlyActiveSurveillance:
    if (in->bparameter["Andreas"]) {
      add_costs("Active surveillance - yearly");
      lost_productivity("Active surveillance - yearly");
    } else {
      if (in->bparameter["MRI_screen"] || in->bparameter["MRI_clinical"]) {
	add_costs("Active surveillance - yearly - with MRI");
	lost_productivity("Active surveillance - yearly - with MRI");
      } else {
	add_costs("Active surveillance - yearly - w/o MRI");
	lost_productivity("Active surveillance - yearly - w/o MRI");
      }
    }
    scheduleAt(now() + 1.0, toYearlyActiveSurveillance);
    break;

  case toYearlyPostTxFollowUp: // not active surveillance
    if (in->bparameter["Andreas"])
      add_costs("Post-Tx follow-up - yearly");
    else {
      lost_productivity("Post-Tx follow-up - yearly");
      if (this->previousFollowup)
	add_costs("Post-Tx follow-up - yearly after");
      else 
	add_costs("Post-Tx follow-up - yearly first");
    }
    previousFollowup = true;
    scheduleAt(now() + 1.0, toYearlyPostTxFollowUp);
    break;

  case toADT:
    if (!in->bparameter("Andreas")) {
      add_costs("ADT+chemo");
    }
    break;

  case toPalliative:
    if (!in->bparameter("Andreas")) {
      add_costs("Palliative therapy - yearly");
    }
    break;

  case toTerminal:
    if (!in->bparameter("Andreas")) {
      add_costs("Terminal illness");
      lost_productivity("Terminal illness");
    }
    break;

  case toUtilityChange:
  case toUtilityRemove:
    {
      utilities->handleMessage(msg);
    } break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;

  } // switch

} // handleMessage()


RcppExport SEXP callFhcrc(SEXP parmsIn) {

  // declarations
  SimInput in;
  SimOutput out;

  in.rngNh = new Rng();
  in.rngOther = new Rng();
  in.rngScreen = new Rng();
  in.rngTreatment = new Rng();
  in.rngSurv = new Rng();
  in.rngBx = new Rng();
  in.rngPrelude = new Rng();
  in.rngNh->set();
  Utilities utilities;

  // read in the parameters
  List parms(parmsIn);
  List tables = parms["tables"];
  in.parameter = parms["parameter"];
  in.bparameter = parms["bparameter"]; // scalar bools
  List otherParameters = parms["otherParameters"];
  in.debug = as<bool>(parms["debug"]);
  if (! in.bparameter["revised_natural_history"]) {
    in.mubeta2 = as<NumericVector>(otherParameters["mubeta2"]);
    in.sebeta2 = as<NumericVector>(otherParameters["sebeta2"]);
  } else {
    in.mubeta2 = as<NumericVector>(otherParameters["rev_mubeta2"]);
    in.sebeta2 = as<NumericVector>(otherParameters["rev_sebeta2"]);
  }
  NumericVector mu0 = as<NumericVector>(otherParameters["mu0"]);
  in.cost_parameters = as<NumericVector>(otherParameters["cost_parameters"]);
  in.utility_estimates = as<NumericVector>(otherParameters["utility_estimates"]);
  in.utility_duration = as<NumericVector>(otherParameters["utility_duration"]);
  utilities.truncate = as<bool>(in.bparameter["utility_truncate"]);
  utilities.scale = utility_scale_t(as<int>(in.parameter["utility_scale"]));

  in.production = Table<double,double>(as<DataFrame>(otherParameters["production"]), "ages", "values");
  in.lost_production_years = as<NumericVector>(otherParameters["lost_production_years"]);

  int n = as<int>(parms["n"]);
  int firstId = as<int>(parms["firstId"]);
  DataFrame background_utilities =
    as<DataFrame>(tables["background_utilities"]);
  in.bg_lower = background_utilities["lower"];
  in.bg_upper = background_utilities["upper"];
  in.bg_utility = background_utilities["utility"];
  
  in.interp_prob_grade7 =
    NumericInterpolate(as<DataFrame>(tables["prob_grade7"]));
  in.prtxCM = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","CM"); // NB: Grade is now {0,1[,2]} coded cf {1,2[,3]}
  in.prtxRP = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","RP");
  in.pradt = TablePradt(as<DataFrame>(tables["pradt"]),"Tx","Age","DxY","Grade","ADT");
  in.hr_locoregional = TableLocoHR(as<DataFrame>(otherParameters["hr_locoregional"]),"age","ext_grade","psa10","hr");
  in.hr_metastatic = TableMetastaticHR(as<DataFrame>(otherParameters["hr_metastatic"]),"age","hr");
  in.tableBiopsySensitivity = TableDD(as<DataFrame>(otherParameters["biopsy_sensitivity"]),"Year","Sensitivity");
  in.tableNegBiopsyToPSAmeanlog = TableDD(as<DataFrame>(otherParameters["neg_biopsy_to_psa"]), "age", "meanlog");
  in.tableNegBiopsyToPSAsdlog = TableDD(as<DataFrame>(otherParameters["neg_biopsy_to_psa"]), "age", "sdlog");
  in.tableNegBiopsyToBiopsymeanlog = TableDD(as<DataFrame>(otherParameters["neg_biopsy_to_biopsy"]), "age", "meanlog");
  in.tableNegBiopsyToBiopsysdlog = TableDD(as<DataFrame>(otherParameters["neg_biopsy_to_biopsy"]), "age", "sdlog");
  in.tableCMtoRPpnever = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RP"]), "age", "pnever");
  in.tableCMtoRPmeanlog = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RP"]), "age", "meanlog");
  in.tableCMtoRPsdlog = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RP"]), "age", "sdlog");
  in.tableCMtoRTpnever = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RT"]), "age", "pnever");
  in.tableCMtoRTmeanlog = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RT"]), "age", "meanlog");
  in.tableCMtoRTsdlog = TableDD(as<DataFrame>(otherParameters["cure_m_CM_to_RT"]), "age", "sdlog");
  in.tableSecularTrendTreatment2008OR = TableDD(as<DataFrame>(tables["secularTrendTreatment2008OR"]),"year","OR");
  in.tableOpportunisticBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyOpportunisticComplianceTable"]),
						"psa","age","compliance");
  in.tableFormalBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyFormalComplianceTable"]),
						"psa","age","compliance");
  in.rescreen_shape = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "shape");
  in.rescreen_scale = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "scale");
  in.rescreen_cure  = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "cure");

  in.H_dist.clear();
  DataFrame df_survival_dist = as<DataFrame>(tables["survival_dist"]); // Grade,Time,Survival
  DataFrame df_survival_local = as<DataFrame>(tables["survival_local"]); // Age,Grade,Time,Survival
  // extract the columns from the survival_dist data-frame
  IntegerVector sd_grades = df_survival_dist["Grade"];
  NumericVector
    sd_times = df_survival_dist["Time"],
    sd_survivals = df_survival_dist["Survival"];
  for (int i=0; i<sd_grades.size(); ++i)
    in.H_dist[sd_grades[i]].push_back(Double(sd_times[i],-log(sd_survivals[i])));
  for (H_dist_t::iterator it_sd = in.H_dist.begin(); it_sd != in.H_dist.end(); it_sd++)
    it_sd->second.prepare();
  // now we can use: H_dist[grade].invert(-log(u))
  in.H_local.clear();
  // extract the columns from the data-frame
  IntegerVector sl_grades = df_survival_local["Grade"];
  NumericVector
    sl_ages = df_survival_local["Age"],
    sl_times = df_survival_local["Time"],
    sl_survivals = df_survival_local["Survival"];
  // push to the map values and set of ages
  for (int i=0; i<sl_grades.size(); ++i) {
    in.H_local_age_set.insert(sl_ages[i]);
    in.H_local[H_local_t::key_type(sl_ages[i],sl_grades[i])].push_back
      (Double(sl_times[i],-log(sl_survivals[i])));
  }
  // prepare the map values for lookup
  for (H_local_t::iterator it_sl = in.H_local.begin();
       it_sl != in.H_local.end();
       it_sl++)
    it_sl->second.prepare();
  // now we can use: H_local[H_local_t::key_type(*H_local_age_set.lower_bound(age),grade)].invert(-log(u))

  if (in.debug) {
    Rprintf("SurvTime: %f\n",exp(-in.H_local[H_local_t::key_type(65.0,0)].approx(63.934032)));
    Rprintf("SurvTime: %f\n",in.H_local[H_local_t::key_type(* in.H_local_age_set.lower_bound(65.0),0)].invert(-log(0.5)));
    Rprintf("SurvTime: %f\n",exp(-in.H_dist[0].approx(5.140980)));
    Rprintf("SurvTime: %f\n",in.H_dist[0].invert(-log(0.5)));
    // Rprintf("Biopsy compliance: %f\n",tableBiopsyCompliance(bounds<double>(1.0,4.0,10.0), bounds<double>(100.0,55,75)));
    Rprintf("Interp for grade 6/7 (expecting approx 0.3): %f\n",in.interp_prob_grade7.approx(0.143));
    Rprintf("prtxCM(80,2008,1) [expecting 0.970711]: %f\n",in.prtxCM(80.0,2008.0,1));
    {
      double age_diag=51.0;
      ext::grade_t ext_grade = ext::Gleason_ge_8;
      // FhrcPerson person = FhcrcPerson(0,1960);
      // Rprintf("hr_localregional(50,0,)=%g\n",hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, person.psamean(age_diag)>10 ? 1 : 0));
      Rprintf("hr_localregional(50,8+,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 0));
      Rprintf("hr_localregional(50,8+,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 1));
      Rprintf("hr_localregional(50,7,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 0));
      Rprintf("hr_localregional(50,7,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 1));
      Rprintf("hr_localregional(50,<=6,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 0));
      Rprintf("hr_localregional(50,<=6,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 1));
      Rprintf("screeningParticipation=%g\n",as<double>(in.parameter["screeningParticipation"]));
    }
  }

  in.nLifeHistories = as<int>(otherParameters["nLifeHistories"]);
  in.screen = as<int>(otherParameters["screen"]);
  if (in.debug) Rprintf("screen=%i\n",in.screen);
  in.panel = as<bool>(parms["panel"]);
  NumericVector cohort = as<NumericVector>(parms["cohort"]); // at present, this is the only chuck-specific data
  bool indiv_reports = as<bool>(in.bparameter["indiv_reports"]);

  // set up the parameters
  double ages0[mu0.size()];
  boost::algorithm::iota(ages0, ages0+mu0.size(), 0.0);
  in.rmu0 = Rpexp(&mu0[0], ages0, mu0.size());
  vector<double> ages(101);
  boost::algorithm::iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // setup for cap_control and cap_study
  // if (in.screen == cap_control || in.screen == cap_study) {
  //   DataFrame uk_screen_uptake = as<DataFrame>(tables["uk_screen_uptake"]); // age,H
  //   in.H_screen_uptake = NumericInterpolate(uk_screen_uptake);
  // }
  if (in.screen == cap_study) {
    in.cap_pScreened = as<NumericVector>(otherParameters["cap_pScreened"]);
  }

  // re-set the output objects
  out.report.clear();
  out.shortReport.clear();
  out.costs.clear();
  out.outParameters.clear();
  out.lifeHistories.clear();
  out.psarecord.clear();
  out.bxrecord.clear();
  out.falsePositives.clear();
  out.diagnoses.clear();

  out.report.discountRate = in.parameter["discountRate.effectiveness"];
  out.report.setPartition(ages);
  out.report.setStartReportAge(in.parameter["startReportAge"]);
  out.shortReport.discountRate = in.parameter["discountRate.effectiveness"];
  out.shortReport.setPartition(ages);
  out.shortReport.setStartReportAge(in.parameter["startReportAge"]);
  out.costs.discountRate = in.parameter["discountRate.costs"];
  out.costs.setPartition(ages);
  out.costs.setStartReportAge(in.parameter["startReportAge"]);
  if (indiv_reports) {
    out.costs.setIndivN(n);
    out.report.setIndivN(n);
    out.shortReport.setIndivN(n);
  }

  // main loop
  FhcrcPerson person(&in, &out, &utilities, 1, 2000, 0);
  for (int i = 0; i < n; ++i) {
    person = FhcrcPerson(&in, &out, &utilities, i+firstId, cohort[i], indiv_reports ? i : 0);
    Sim::create_process(&person);
    Sim::run_simulation();
    Sim::clear();
    in.rngNh->nextSubstream();
    in.rngOther->nextSubstream();
    in.rngScreen->nextSubstream();
    in.rngTreatment->nextSubstream();
    in.rngSurv->nextSubstream();
    in.rngBx->nextSubstream();
    in.rngPrelude->nextSubstream();
    if (i % 10000 == 0) Rcpp::checkUserInterrupt(); /* be polite -- did the user hit ctrl-C? */
  }

  // output
  // TODO: clean up these objects in C++ (cf. R)
  return List::create(_("costs") = out.costs.wrap(),                // CostReport
		      _("summary") = out.report.wrap(),             // EventReport
		      _("shortSummary") = out.shortReport.wrap(),   // EventReport
		      _("lifeHistories") = wrap(out.lifeHistories), // vector<LifeHistory::Type>
		      _("parameters") = out.outParameters.wrap(),   // SimpleReport<double>
		      _("psarecord")=out.psarecord.wrap(),          // SimpleReport<double>
		      _("bxrecord")=out.bxrecord.wrap(),            // SimpleReport<double>
		      _("falsePositives")=out.falsePositives.wrap(),// SimpleReport<double>
		      _("diagnoses")=out.diagnoses.wrap(),          // SimpleReport<double>
		      _("tmc_minus_t0")=out.tmc_minus_t0.wrap(),    // Means
		      _("indiv_costs")=out.costs.wrap_indiv(),      // vector<double>
		      _("indiv_utilities")=out.report.wrap_indiv(), // vector<double>
		      _("mean_utilities")=(in.parameter["full_report"] == 1.0) ? out.report.wrap_means() : out.shortReport.wrap_means(),  // Rcpp::DataFrame
		      _("mean_costs")=out.costs.wrap_means()        // Rcpp::DataFrame
		      );
}

} // anonymous namespace
