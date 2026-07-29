#include <gtest/gtest.h>

#include <Rembedded.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

// Compile the production unit into this test target so we can access
// SimInput/SimOutput and internal helpers directly.
#include "../../src/prostata.cpp"

namespace {

Rcpp::List evalList(const std::string& expr) {
  Rcpp::Function parse("parse", R_BaseEnv);
  Rcpp::Function eval("eval", R_BaseEnv);
  Rcpp::ExpressionVector parsed = parse(Rcpp::Named("text", expr));
  return eval(parsed);
}

Rcpp::List getCachedSimParms() {
  static const std::string parms_expr = R"RSCRIPT(
local({
  sim <- suppressWarnings(prostata::callFhcrc(
    n = 10,
    screen = "noScreening",
    mc.cores = 1,
    print.timing = FALSE
  ))
  parameter <- sim$simulation.parameters
  pind <- sapply(parameter, class) == "numeric" & sapply(parameter, length) == 1
  bInd <- sapply(parameter, class) == "logical" & sapply(parameter, length) == 1
  list(
    panel = FALSE,
    debug = FALSE,
    cohort = as.double(rep(1950.0, 10)),
    parameter = unlist(parameter[pind]),
    bparameter = unlist(parameter[bInd]),
    otherParameters = parameter[!pind & !bInd]
  )
})
)RSCRIPT";

  static Rcpp::List cached = evalList(parms_expr);
  return Rcpp::clone(cached);
}

Rcpp::IntegerVector eventCountsFromShortReport(fhcrc_example::SimOutput& output) {
  enum EventIndex {
    kToScreen,
    kToSTHLM3,
    kToMRI,
    kToClinicalDiagnosis,
    kToScreenDiagnosis,
    kToCancerDeath,
    kToOverDiagnosis,
    kToOpportunistic
  };

  Rcpp::IntegerVector aggregated(8);
  aggregated.names() = Rcpp::CharacterVector::create(
    "toScreen",
    "toSTHLM3",
    "toMRI",
    "toClinicalDiagnosis",
    "toScreenDiagnosis",
    "toCancerDeath",
    "toOverDiagnosis",
    "toOpportunistic"
  );

  for (const auto& life_history : output.lifeHistories) {
    const short event_code = std::get<fhcrc_example::LifeHistory::event>(life_history);
    switch (event_code) {
      case fhcrc_example::toScreen:
        aggregated[kToScreen] += 1;
        break;
      case fhcrc_example::toSTHLM3:
        aggregated[kToSTHLM3] += 1;
        break;
      case fhcrc_example::toMRI:
        aggregated[kToMRI] += 1;
        break;
      case fhcrc_example::toClinicalDiagnosis:
        aggregated[kToClinicalDiagnosis] += 1;
        break;
      case fhcrc_example::toScreenDiagnosis:
        aggregated[kToScreenDiagnosis] += 1;
        break;
      case fhcrc_example::toCancerDeath:
        aggregated[kToCancerDeath] += 1;
        break;
      case fhcrc_example::toOverDiagnosis:
        aggregated[kToOverDiagnosis] += 1;
        break;
      case fhcrc_example::toOpportunistic:
        aggregated[kToOpportunistic] += 1;
        break;
    }
  }

  return aggregated;
}

class EmbeddedREnvironment : public ::testing::Environment {
public:
  static void loadPackage(const char* package_name) {
    SEXP call = PROTECT(Rf_lang2(Rf_install("library"), Rf_mkString(package_name)));
    Rf_eval(call, R_GlobalEnv);
    UNPROTECT(1);
  }

  void SetUp() override {
    const char* r_home = std::getenv("R_HOME");
    if (r_home == nullptr || std::strlen(r_home) == 0) {
      FILE* fp = popen("R RHOME", "r");
      if (fp != nullptr) {
        char buffer[1024];
        if (fgets(buffer, sizeof(buffer), fp) != nullptr) {
          std::size_t len = std::strlen(buffer);
          while (len > 0 && (buffer[len - 1] == '\n' || buffer[len - 1] == '\r')) {
            buffer[len - 1] = '\0';
            --len;
          }
          if (len > 0) {
            setenv("R_HOME", buffer, 1);
          }
        }
        pclose(fp);
      }
    }

    if (std::getenv("R_HOME") == nullptr || std::strlen(std::getenv("R_HOME")) == 0) {
      FAIL() << "R_HOME is not set and could not be inferred for embedded R.";
    }

    int argc = 4;
    char arg0[] = "R";
    char arg1[] = "--silent";
    char arg2[] = "--no-save";
    char arg3[] = "--no-restore";
    char* argv[] = {arg0, arg1, arg2, arg3};
    Rf_initEmbeddedR(argc, argv);

    // Ensure Rcpp C-callables are registered before any Rcpp wrapper objects are created.
    loadPackage("Rcpp");
    loadPackage("microsimulation");
    loadPackage("prostata");
  }

  void TearDown() override {
    Rf_endEmbeddedR(0);
  }
};

::testing::Environment* const embedded_r_env =
  ::testing::AddGlobalTestEnvironment(new EmbeddedREnvironment());

} // namespace

class CallFhcrc : public ::testing::TestWithParam<std::tuple<int, int>> {
protected:
  void SetUp() override {
    // Ensure that the R environment is initialized before any tests run.
    ASSERT_NE(embedded_r_env, nullptr);
  }
};

INSTANTIATE_TEST_SUITE_P(
  CallFhcrcTests,
  CallFhcrc,
  ::testing::Values(
    std::make_tuple(10, 1),
    std::make_tuple(10, 2),
    std::make_tuple(10, 3),
    std::make_tuple(10, 4),
    std::make_tuple(10, 5),
    std::make_tuple(10, 6),
    std::make_tuple(10, 7),
    std::make_tuple(10, 8)
  )
);

void RunPerThreadLoopCaseWithAssertions(const int n, const int numThreads) {
  using namespace fhcrc_example;

  Rcpp::List parms = getCachedSimParms();
  SimInput in = initialize(parms);

  // enable RNG
  std::array<double, 6> currentSeed {12345,12345,12345,12345,12345,12345};
  set_user_random_seed(currentSeed);
  r_create_current_stream();
  //Advance RNG with U01 to mimic "sample" call in R code
  // cohort <- sample(pop$cohort,n,prob=pop$pop/sum(pop$pop),replace=TRUE)
  for(int i=0; i<n; ++i) {
    user_unif_rand();
  }

  Rcpp::NumericVector cohort = Rcpp::as<Rcpp::NumericVector>(parms["cohort"]);
  ASSERT_EQ(cohort.size(), n);

  auto initialSeeds = getInitialSeeds(numThreads, n);
  
  std::vector<SimOutput> outputs(numThreads);
  for(int i=0; i<numThreads; i++) {
      set_user_random_seed(initialSeeds[i]);
      in.resetRngs();
      int firstId = static_cast<int>(std::floor(static_cast<double>(i)/static_cast<double>(numThreads)*n));
      int nextId = static_cast<int>(std::floor(static_cast<double>(i+1)/static_cast<double>(numThreads)*n));
      outputs[i] = callFhcrc_inner(in, nextId - firstId, firstId, cohort);
  }

  // aggregate outputs from each thread
  SimOutput out = outputs[0];
  for(int i=1; i<numThreads; ++i) {
      // out.costs.append(outputs[i].costs);
      out.report.append(outputs[i].report);
      out.shortReport.append(outputs[i].shortReport);
      out.lifeHistories.insert(out.lifeHistories.end(), outputs[i].lifeHistories.begin(), outputs[i].lifeHistories.end());
      out.outParameters.append(outputs[i].outParameters);
      out.psarecord.append(outputs[i].psarecord);
      out.bxrecord.append(outputs[i].bxrecord);
      out.falsePositives.append(outputs[i].falsePositives);
      out.diagnoses.append(outputs[i].diagnoses);
      out.tmc_minus_t0.combine(outputs[i].tmc_minus_t0);
      (in.parameter("full_report") == 1.0) ? out.report.mean_utilities.combine(outputs[i].report.mean_utilities) : out.shortReport.mean_utilities.combine(outputs[i].shortReport.mean_utilities);
      out.costs.mean_costs.combine(outputs[i].costs.mean_costs);
  }

  ASSERT_EQ(in.screen, static_cast<int>(noScreening));

  ASSERT_DOUBLE_EQ(in.parameter("full_report"), 1.0);

  Rcpp::IntegerVector actual_event_counts = eventCountsFromShortReport(out);
  Rcpp::DataFrame actual_mean_utilities = out.report.wrap_means();
  Rcpp::DataFrame actual_mean_costs = out.costs.wrap_means();

  Rcpp::NumericVector utility_n = actual_mean_utilities["n"];
  Rcpp::NumericVector utility_mean = actual_mean_utilities["mean"];
  Rcpp::NumericVector utility_var = actual_mean_utilities["var"];
  Rcpp::NumericVector utility_sd = actual_mean_utilities["sd"];
  Rcpp::NumericVector utility_se = actual_mean_utilities["se"];
  Rcpp::NumericVector utility_sum = actual_mean_utilities["sum"];
  Rcpp::NumericVector utility_sumsq = actual_mean_utilities["sumsq"];

  Rcpp::NumericVector cost_n = actual_mean_costs["n"];
  Rcpp::NumericVector cost_mean = actual_mean_costs["mean"];
  Rcpp::NumericVector cost_var = actual_mean_costs["var"];
  Rcpp::NumericVector cost_sd = actual_mean_costs["sd"];
  Rcpp::NumericVector cost_se = actual_mean_costs["se"];
  Rcpp::NumericVector cost_sum = actual_mean_costs["sum"];
  Rcpp::NumericVector cost_sumsq = actual_mean_costs["sumsq"];

  EXPECT_EQ(actual_event_counts[0], 0);
  EXPECT_EQ(actual_event_counts[1], 0);
  EXPECT_EQ(actual_event_counts[2], 0);
  EXPECT_EQ(actual_event_counts[3], 1);
  EXPECT_EQ(actual_event_counts[4], 0);
  EXPECT_EQ(actual_event_counts[5], 0);
  EXPECT_EQ(actual_event_counts[6], 0);
  EXPECT_EQ(actual_event_counts[7], 0);

  EXPECT_DOUBLE_EQ(utility_n[0], 10.0);
  EXPECT_DOUBLE_EQ(utility_mean[0], 27.267603597388469);
  EXPECT_DOUBLE_EQ(utility_var[0], 7.7947179660890438);
  EXPECT_DOUBLE_EQ(utility_sd[0], 2.7919022128450424);
  EXPECT_DOUBLE_EQ(utility_se[0], 0.88287699970545408);
  EXPECT_DOUBLE_EQ(utility_sum[0], 272.67603597388467);
  EXPECT_DOUBLE_EQ(utility_sumsq[0], 7505.3745211379264);

  EXPECT_DOUBLE_EQ(cost_n[0], 10.0);
  EXPECT_DOUBLE_EQ(cost_mean[0], 8.4759924384488858);
  EXPECT_DOUBLE_EQ(cost_var[0], 718.42447816642687);
  EXPECT_DOUBLE_EQ(cost_sd[0], 26.803441535863019);
  EXPECT_DOUBLE_EQ(cost_se[0], 8.4759924384488858);
  EXPECT_DOUBLE_EQ(cost_sum[0], 84.759924384488855);
  EXPECT_DOUBLE_EQ(cost_sumsq[0], 7184.244781664268);
}

TEST_P(CallFhcrc, ReproducesPerThreadLoopWithInitializeAndInnerCall) {
  const int n = std::get<0>(GetParam());
  const int numThreads = std::get<1>(GetParam());

  RunPerThreadLoopCaseWithAssertions(n, numThreads);
}
