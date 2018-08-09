#ifndef WRAP_OPTIM_HPP
#define WRAP_OPTIM_HPP

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <iomanip>
#include "optim.hpp"
#include "common.hpp"

namespace optim {

typedef arma::vec native;

struct Settings : optim::algo_settings {
  common::parameters start;
  std::string algorithm{"de"};
  std::string file{"./output/optim.dat"};
  bool write_to_disk{true};
  bool overwrite{true};
  bool verbose{true};
  int n_dims{10};
  bool default_limits{true};
  bool default_start{true};
  bool write_unphysical{false};
};

void run(common::loglike loglike,
         common::prior prior,
         Settings settings) {
  /*
    Wrapper to optim library.
  */

  std::ofstream file;
  if (settings.write_to_disk) {
    file << std::setprecision(20);
    if (settings.overwrite) {
      file.open(settings.file.c_str());
    } else {
      file.open(settings.file.c_str(), std::ios_base::app);
    }
  }
  const char sep = ' ';
  double threshold = -std::numeric_limits<double>::max();

  auto combined = [loglike, prior, settings, &file, sep, threshold](native x, native* grad, void* context) mutable {
    const auto physical = prior(arma::conv_to<common::parameters>::from(x));
    common::pair pair = loglike(physical, threshold);
    double loglike_ = pair.second;
    threshold = std::max(threshold, loglike_);
    auto extended = pair.first;

    if (settings.write_to_disk && (settings.write_unphysical || (loglike_ > common::unphysical))) {
      file << loglike_ << sep;
      for (common::parameters::const_iterator p = extended.begin(); p != extended.end(); ++p) {
        file << *p << sep;
      }
      file << std::endl;
    }
    return -loglike_;
  };

  if (settings.default_limits) {
    settings.lower_bounds = arma::zeros(settings.n_dims, 1);
    settings.upper_bounds = arma::ones(settings.n_dims, 1);
    settings.vals_bound = true;
  }

  auto x = settings.default_start ?
    arma::zeros(settings.n_dims, 1) + 0.5 : native(settings.start);

  if (settings.algorithm == "de") {
    optim::de(x, combined, nullptr, settings);
  } else if (settings.algorithm == "pso") {
    optim::pso(x, combined, nullptr, settings);
  } else if (settings.algorithm == "nm") {
    optim::nm(x, combined, nullptr, settings);
  } else {
    throw std::invalid_argument("unknown algorithm");
  }

  file.close();

  if (settings.verbose) {
    std::cout << "log-likelihood = " << settings.opt_value << std::endl;
    std::cout << "x = ";
    const auto physical = prior(arma::conv_to<common::parameters>::from(x));
    for (common::parameters::const_iterator p = physical.begin(); p != physical.end(); ++p) {
      std::cout << *p << std::endl << "    ";
    }
    std::cout << std::endl;
  }
}

}  // namespace optim

#endif  // WRAP_OPTIM_HPP
