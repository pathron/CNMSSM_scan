#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "multinest.hpp"
#include "priors.hpp"
#include "common.hpp"

struct Settings : nested::Settings {
  using nested::Settings::n_dims;
  using nested::Settings::write_to_disk;
  using nested::Settings::overwrite;
  using nested::Settings::verbose;
};

#endif  // SAMPLE_HPP
