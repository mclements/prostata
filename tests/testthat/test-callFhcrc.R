library(testthat)
library(prostata)

context("callFhcrc regression")

update_fixtures <- tolower(Sys.getenv("PROSTATA_UPDATE_FIXTURES", "false")) %in% c("1", "true", "yes")
full_fixture_dir <- file.path("tests", "testthat", "fixtures", "callFhcrc-full-v1")

scenario_params <- function(screen) {
    if (screen %in% c("cap_control", "cap_study")) {
        return(list(cap_pScreened = c(0.3235438, 0.3531057, 0.3582405, 0.3236445)))
    }
    if (screen == "germany_2018") {
        return(prostata:::TrustParameters())
    }
    NULL
}

scenario_signatures <- list(
    noScreening = list(screen_index = 0L, event_counts = c(toScreen = 0L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    randomScreen50to70 = list(screen_index = 1L, event_counts = c(toScreen = 62L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 1L, toCancerDeath = 5L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    twoYearlyScreen50to70 = list(screen_index = 2L, event_counts = c(toScreen = 480L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 6L, toScreenDiagnosis = 7L, toCancerDeath = 4L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    fourYearlyScreen50to70 = list(screen_index = 3L, event_counts = c(toScreen = 299L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 8L, toScreenDiagnosis = 5L, toCancerDeath = 5L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    screen50 = list(screen_index = 4L, event_counts = c(toScreen = 59L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    screen60 = list(screen_index = 5L, event_counts = c(toScreen = 60L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    screen70 = list(screen_index = 6L, event_counts = c(toScreen = 58L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 9L, toScreenDiagnosis = 4L, toCancerDeath = 5L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    screenUptake = list(screen_index = 7L, event_counts = c(toScreen = 433L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 7L, toScreenDiagnosis = 6L, toCancerDeath = 4L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    stockholm3_goteborg = list(screen_index = 8L, event_counts = c(toScreen = 67L, toSTHLM3 = 2L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    stockholm3_risk_stratified = list(screen_index = 9L, event_counts = c(toScreen = 118L, toSTHLM3 = 2L, toMRI = 0L, toClinicalDiagnosis = 10L, toScreenDiagnosis = 2L, toCancerDeath = 5L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    goteborg = list(screen_index = 10L, event_counts = c(toScreen = 59L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    risk_stratified = list(screen_index = 11L, event_counts = c(toScreen = 206L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 8L, toScreenDiagnosis = 5L, toCancerDeath = 5L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    mixed_screening = list(screen_index = 12L, event_counts = c(toScreen = 526L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 7L, toScreenDiagnosis = 5L, toCancerDeath = 5L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    regular_screen = list(screen_index = 13L, event_counts = c(toScreen = 472L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 6L, toScreenDiagnosis = 7L, toCancerDeath = 4L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    single_screen = list(screen_index = 14L, event_counts = c(toScreen = 59L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 0L, toCancerDeath = 5L, toOverDiagnosis = 0L, toOpportunistic = 0L)),
    introduced_screening_only = list(screen_index = 15L, event_counts = c(toScreen = 474L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 6L, toScreenDiagnosis = 6L, toCancerDeath = 4L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    introduced_screening_preference = list(screen_index = 16L, event_counts = c(toScreen = 503L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 7L, toScreenDiagnosis = 5L, toCancerDeath = 4L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    introduced_screening = list(screen_index = 17L, event_counts = c(toScreen = 556L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 7L, toScreenDiagnosis = 5L, toCancerDeath = 4L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    stopped_screening = list(screen_index = 18L, event_counts = c(toScreen = 173L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 11L, toScreenDiagnosis = 1L, toCancerDeath = 5L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    cap_control = list(screen_index = 19L, event_counts = c(toScreen = 449L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 6L, toScreenDiagnosis = 8L, toCancerDeath = 4L, toOverDiagnosis = 3L, toOpportunistic = 0L)),
    cap_study = list(screen_index = 20L, event_counts = c(toScreen = 375L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 7L, toScreenDiagnosis = 6L, toCancerDeath = 5L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    sthlm3_mri_arm = list(screen_index = 21L, event_counts = c(toScreen = 148L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 10L, toScreenDiagnosis = 2L, toCancerDeath = 5L, toOverDiagnosis = 1L, toOpportunistic = 0L)),
    grs_stratified = list(screen_index = 22L, event_counts = c(toScreen = 294L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 8L, toScreenDiagnosis = 5L, toCancerDeath = 4L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    grs_stratified_age = list(screen_index = 23L, event_counts = c(toScreen = 294L, toSTHLM3 = 0L, toMRI = 0L, toClinicalDiagnosis = 8L, toScreenDiagnosis = 5L, toCancerDeath = 4L, toOverDiagnosis = 2L, toOpportunistic = 0L)),
    germany_2018 = list(screen_index = 24L, event_counts = c(toScreen = 402L, toSTHLM3 = 0L, toMRI = 78L, toClinicalDiagnosis = 9L, toScreenDiagnosis = 5L, toCancerDeath = 4L, toOverDiagnosis = 1L, toOpportunistic = 0L))
)

run_signature <- function(screen) {
    sim <- suppressWarnings(callFhcrc(
        n = 80,
        screen = screen,
        seed = 12345,
        mc.cores = 1,
        print.timing = FALSE,
        parms = scenario_params(screen)
    ))

    ev <- sim$summary$events
    event_counts <- c(
        toScreen = as.integer(sum(ev$n[ev$event == "toScreen"])),
        toSTHLM3 = as.integer(sum(ev$n[ev$event == "toSTHLM3"])),
        toMRI = as.integer(sum(ev$n[ev$event == "toMRI"])),
        toClinicalDiagnosis = as.integer(sum(ev$n[ev$event == "toClinicalDiagnosis"])),
        toScreenDiagnosis = as.integer(sum(ev$n[ev$event == "toScreenDiagnosis"])),
        toCancerDeath = as.integer(sum(ev$n[ev$event == "toCancerDeath"])),
        toOverDiagnosis = as.integer(sum(ev$n[ev$event == "toOverDiagnosis"])),
        toOpportunistic = as.integer(sum(ev$n[ev$event == "toOpportunistic"]))
    )

    list(
        class = class(sim),
        screen_index = as.integer(sim$simulation.parameters$screen),
        event_counts = event_counts
    )
}

run_full_sim <- function(screen) {
    suppressWarnings(callFhcrc(
        n = 80,
        screen = screen,
        seed = 12345,
        mc.cores = 1,
        print.timing = FALSE,
        parms = scenario_params(screen)
    ))
}

strip_start_attr <- function(x) {
    if (is.list(x)) {
        for (i in seq_along(x)) {
            x[[i]] <- strip_start_attr(x[[i]])
        }
        return(x)
    }

    if (!is.null(attributes(x)) && "start" %in% names(attributes(x))) {
        attr(x, "start") <- NULL
    }

    x
}

expect_same_type_tree <- function(got, expected, path = "root") {
    expect_identical(typeof(got), typeof(expected), info = paste(path, "typeof"))
    expect_identical(class(got), class(expected), info = paste(path, "class"))

    if (is.list(got) && is.list(expected)) {
        got_names <- names(got)
        expected_names <- names(expected)
        expect_identical(got_names, expected_names, info = paste(path, "names"))

        if (!is.null(got_names)) {
            for (nm in got_names) {
                expect_same_type_tree(got[[nm]], expected[[nm]], paste0(path, "$", nm))
            }
        } else {
            for (i in seq_along(got)) {
                expect_same_type_tree(got[[i]], expected[[i]], paste0(path, "[[", i, "]]"))
            }
        }
    }
}

fixture_path <- function(screen) {
    file.path(full_fixture_dir, paste0(screen, ".rds"))
}

test_that("All screening scenarios match regression signatures", {
    skip_on_cran()

    for (screen in names(scenario_signatures)) {
        got <- run_signature(screen)
        expected <- scenario_signatures[[screen]]

        expect_equal(got$class, "fhcrc", info = screen)
        expect_equal(got$screen_index, expected$screen_index, info = c("screen", screen, "got index", got$screen_index, "expected index", expected$screen_index))
        expect_equal(got$event_counts, expected$event_counts, info = c("screen", screen, "got counts", got$event_counts, "expected counts", expected$event_counts))
    }
})

test_that("All screening scenarios match full sim list fixtures", {
    skip_on_cran()

    if (update_fixtures) {
        dir.create(full_fixture_dir, recursive = TRUE, showWarnings = FALSE)
    }

    for (screen in names(scenario_signatures)) {
        path <- fixture_path(screen)
        got <- run_full_sim(screen)

        if (update_fixtures || !file.exists(path)) {
            if (!update_fixtures) {
                fail(paste0(
                    "Missing regression fixture for screen '", screen, "': ", path,
                    ". Re-run tests with PROSTATA_UPDATE_FIXTURES=true to create fixtures."
                ))
            }
            saveRDS(got, path, version = 2)
        }

        expected <- readRDS(path)
        expect_same_type_tree(got, expected, path = paste0("screen=", screen))
        expect_equal(strip_start_attr(got), strip_start_attr(expected), info = screen)
    }
})
