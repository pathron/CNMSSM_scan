/*
  Example of using general algorithm.
*/

#include "sample.hpp"

common::pair loglike(common::parameters x, double threshold) {
  const double a = x[0];
  const double b = x[1];
  double loglike_ = -(pow(a - 0.42, 2) + pow(b - 0.57, 2));
  return {x, loglike_};
}

common::parameters prior(common::parameters x) {
  x[0] = priors::flat(x[0], 10., 100.);
  return x;
}

int main() {
  auto settings = Settings();
  optim::run(loglike, prior, settings);
  nested::run(loglike, prior, settings);
  return 0;
}
