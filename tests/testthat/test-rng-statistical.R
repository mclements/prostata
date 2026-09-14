library(testthat)
library(prostata)

context("RNG statistical equivalence")

rng_reference <- data.frame(
    endpoint = c(
        "clinical_diagnosis_risk",
        "cancer_death_risk",
        "discounted_utility",
        "discounted_cost"
    ),
    mean = c(
        0.1549166666666667,
        0.05258333333333334,
        27.9481397000961884,
        310.07395801365288
    ),
    standard_error = c(
        0.0024963137073701,
        0.00196888538565244,
        0.0124429210801198,
        9.43291089637793
    ),
    equivalence_margin = c(0.02, 0.01, 0.15, 75.0),
    row.names = NULL
)

run_rng_statistical_sample <- function(seeds = 11001:11024, cohort_size = 500L) {
    run_one <- function(seed) {
        sim <- suppressWarnings(callFhcrc(
            n = cohort_size,
            nLifeHistories = 0,
            screen = "noScreening",
            seed = seed,
            mc.cores = 1,
            print.timing = FALSE,
            parms = list(full_report = 0)
        ))

        events <- sim$summary$events
        event_risk <- function(event) {
            sum(events$n[events$event == event]) / cohort_size
        }

        c(
            clinical_diagnosis_risk = event_risk("toClinicalDiagnosis"),
            cancer_death_risk = event_risk("toCancerDeath"),
            discounted_utility = sim$mean_utilities$mean[[1]],
            discounted_cost = sim$mean_costs$mean[[1]]
        )
    }

    t(vapply(seeds, run_one, numeric(nrow(rng_reference))))
}

test_that("RNG replacement preserves important outcome distributions", {
    skip_on_cran()

    samples <- run_rng_statistical_sample()
    candidate_mean <- colMeans(samples)
    candidate_se <- apply(samples, 2, sd) / sqrt(nrow(samples))
    combined_se <- sqrt(candidate_se^2 + rng_reference$standard_error^2)
    observed_difference <- abs(candidate_mean - rng_reference$mean)
    difference_upper_bound <- observed_difference +
        qnorm(0.995) * combined_se

    for (index in seq_len(nrow(rng_reference))) {
        expect_true(
            unname(difference_upper_bound[index]) <=
                unname(rng_reference$equivalence_margin[index]),
            info = paste(
                rng_reference$endpoint[index],
                "mean:", format(candidate_mean[index], digits = 8),
                "reference:", format(rng_reference$mean[index], digits = 8),
                "99% difference bound:",
                format(difference_upper_bound[index], digits = 8),
                "equivalence margin:",
                format(rng_reference$equivalence_margin[index], digits = 8)
            )
        )
    }
})