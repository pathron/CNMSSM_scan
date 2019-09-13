

#ifndef PRIORS_HPP
#define PRIORS_HPP

#include <boost/math/distributions/normal.hpp>


namespace priors {

  double flat(double unit, double min, double max) {
    return (max - min) * unit + min;
  }

  double log(double unit, double min, double max) {
    return min * pow(max / min, unit);
  }

  double gaussian(double unit, double mu, double sigma) {
    boost::math::normal dist(mu, sigma);
    return quantile(dist, unit);
  }

  double signed_log(double unit, double min, double max) {
    double mapped = (unit < 0.5) ? 2. * unit : 2. * unit - 1.;
    double sign = (unit < 0.5) ? -1. : 1.;
    return sign * min * pow(max / min, mapped);
  }

}

#endif  // PRIORS_HPP
