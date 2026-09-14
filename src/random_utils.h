#pragma once

#include "microsimulation_patched.h"

namespace fhcrc_example {

double rlnorm(double meanlog, double sdlog, const std::unique_ptr<ssim::Rng>& rng);
double rnorm(double mu, double sigma, const std::unique_ptr<ssim::Rng>& rng);
double runif(double a, double b, const std::unique_ptr<ssim::Rng>& rng);

} // namespace fhcrc_example
