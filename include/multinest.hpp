#ifndef MULTINEST_HPP
#define MULTINEST_HPP

#ifdef __INTEL_COMPILER
  #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__
  #define NESTRUN __nested_MOD_nestrun
#endif

#include <limits>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include "common.hpp"


namespace nested {

struct Settings {
  bool IS_mode{true};
  bool multimodal_mode{false};
  bool ceff_mode{false};
  int n_live{1000};
  double efficiency{0.8};
  double tol{0.5};
  int n_dims{2};
  int n_pars{0};
  int n_separate{2};
  int n_update_every{1000};
  double mode_evidence_tol{-std::numeric_limits<double>::max()};
  int n_max_modes{100};
  std::string root{"./output/multinest-"};
  int seed{-1};
  bool verbose{true};
  bool overwrite{false};
  bool write_to_disk{true};
  bool initialize_MPI{true};
  bool write_calls_to_disk{true};
  bool write_unphysical{false};
  std::string file{"./output/multinest.dat"};
  double log_zero{-std::numeric_limits<double>::max()};
  int n_max_iter{0};
};

extern "C" {
  void NESTRUN(int &IS_mode,
               int &multimodal_mode,
               int &ceff_mode,
               int &n_live,
               double &tol,
               double &efficiency,
               int &n_dims,
               int &nPar,
               int &n_separate,
               int &n_max_modes,
               int &n_update_every,
               double &mode_evidence_tol,
               char *root,
               int &seed,
               int *pWrap,
               int &verbose,
               int &resume,
               int &write_to_disk,
               int &initialize_MPI,
               double &log_zero,
               int &n_max_iter,
               void (*loglike)(double *, int *, int *, double *, void *),
               void (*dumper)(int *, int *, int *, double **, double **,
                 double **, double *, double *, double *, double *, void *),
               void *context,
               int &root_len);
}

void dumper(int *nSamples,
            int *n_live,
            int *nPar,
            double **physLive,
            double **posterior,
            double **paramConstr,
            double *maxLogLike,
            double *logZ,
            double *INSlogZ,
            double *logZerr,
            void *context) {
}

static void to_c_api(double* cube, int *n_dims, int *n_pars, double *loglike_, void* context) {
  auto threshold = *loglike_;
  auto combined = static_cast<common::loglike*>(context);
  common::parameters unphysical(cube, cube + *n_dims);
  common::pair pair = (*combined)(unphysical, threshold);
  *loglike_ = pair.second;

  if (static_cast<unsigned int>(*n_pars) < pair.first.size()) {
    throw std::invalid_argument("n_pars < size of cube");
  }

  for (unsigned int i = 0; i < pair.first.size(); i += 1) {
    cube[i] = pair.first[i];
  }
}

void run(common::loglike loglike,
         common::prior prior,
         Settings settings) {
  /*
    Wrapper to MultiNest with C++ style arguments.
  */

  std::ofstream file;
  if (settings.write_calls_to_disk) {
    file << std::setprecision(20);
    if (settings.overwrite) {
      file.open(settings.file.c_str());
    } else {
      file.open(settings.file.c_str(), std::ios_base::app);
    }
  }

  const char sep = ' ';

  common::loglike combined = [prior, loglike, settings, &file, sep](common::parameters x, double threshold) {
    auto physical = prior(x);
    common::pair pair = loglike(physical, threshold);
    double loglike_ = pair.second;
    auto extended = pair.first;
    if (settings.write_calls_to_disk && (settings.write_unphysical || (loglike_ > common::unphysical))) {
      file << loglike_ << sep;
      for (common::parameters::const_iterator p = extended.begin(); p != extended.end(); ++p) {
        file << *p << sep;
      }
      file << std::endl;
    }
    return pair;
  };

  // Guess number of parameters

  if (settings.n_pars == 0) {
    common::parameters x(settings.n_dims);
    auto trial = combined(x, -std::numeric_limits<double>::max());
    settings.n_pars = trial.first.size();
  }

  int IS_mode = static_cast<int>(settings.IS_mode);
  int multimodal_mode = static_cast<int>(settings.multimodal_mode);
  int ceff_mode = static_cast<int>(settings.ceff_mode);
  int verbose = static_cast<int>(settings.verbose);
  int resume = static_cast<int>(!settings.overwrite);
  int write_to_disk = static_cast<int>(settings.write_to_disk);
  int initialize_MPI = static_cast<int>(settings.initialize_MPI);
  int root_len = settings.root.length();

  constexpr int kLen = 1000;
  auto pad = settings.root;
  pad.append(kLen, ' ');

  int wrap[settings.n_dims];
  for (int i = 0; i < settings.n_dims; i += 1) {
    wrap[i] = 0.;
  }

  NESTRUN(IS_mode,
          multimodal_mode,
          ceff_mode,
          settings.n_live,
          settings.tol,
          settings.efficiency,
          settings.n_dims,
          settings.n_pars,
          settings.n_separate,
          settings.n_max_modes,
          settings.n_update_every,
          settings.mode_evidence_tol,
          const_cast<char*>(pad.c_str()),
          settings.seed,
          wrap,
          verbose,
          resume,
          write_to_disk,
          initialize_MPI,
          settings.log_zero,
          settings.n_max_iter,
          to_c_api,
          dumper,
          &combined,
          root_len);
}

}  // namespace nested

#endif  // MULTINEST_HPP
