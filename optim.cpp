/*
  Example of using optim algorithm.
*/

#include "wrap_optim.hpp"

common::pair loglike(common::parameters x, double threshold) {
  const double a = x[0];
  const double b = x[1];
  double loglike_ = -(pow(a - 0.42, 2) + pow(b - 0.57, 2));
  x.push_back(100.); // Add something to be printed.
  return {x, loglike_};
}

common::parameters prior(common::parameters x) {
  return x;
}

int main() {
  auto settings = optim::Settings();
  optim::run(loglike, prior, settings);
  return 0;
}
