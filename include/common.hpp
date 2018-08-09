#ifndef COMMON_HPP
#define COMMON_HPP

#include <functional>
#include <vector>
#include <utility>
#include <limits>


namespace common {
  typedef std::vector<double> parameters;
  typedef std::pair<common::parameters, double> pair;
  typedef std::function<pair(parameters, double)> loglike;
  typedef std::function<parameters(parameters)> prior;
  double unphysical = -std::numeric_limits<double>::max();
}

#endif  // COMMON_HPP
