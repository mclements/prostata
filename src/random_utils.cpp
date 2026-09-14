#include "random_utils.h"

#include <cmath>
#include <random>
#include <stdexcept>

#include <mutex>
#include "Rmath.h"

namespace fhcrc_example {
  
  std::mutex mtx;

  double rlnorm(double meanlog, double sdlog, const std::unique_ptr<ssim::Rng>& rng) {
    std::lock_guard<std::mutex> lock(mtx);
    rng->set();
    return R::rlnorm(meanlog, sdlog);
  }
  double rnorm(double mu, double sigma, const std::unique_ptr<ssim::Rng>& rng) {
    std::lock_guard<std::mutex> lock(mtx);
    rng->set();
    return R::rnorm(mu, sigma);
  }
  double runif(double a, double b, const std::unique_ptr<ssim::Rng>& rng) {
    std::lock_guard<std::mutex> lock(mtx);
    rng->set();
    return R::runif(a, b);
  }

}  // namespace fhcrc_example
