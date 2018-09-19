/*
  Example multinest scan of FS-> Lilith.
*/

#include "sample.hpp"
#include "higgs.hpp"
#include "common.hpp"
#include <cmath>

#include "config.h"
#include "CNMSSM_input_parameters.hpp"
#include "CNMSSM_effective_couplings.hpp"
#include "CNMSSM_slha_io.hpp"
#include "CNMSSM_spectrum_generator.hpp"
#include "CNMSSM_semi_analytic_spectrum_generator.hpp"

#include "standard_model.hpp"
#include "standard_model_effective_couplings.hpp"

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"
#include "wrappers.hpp"
#include <iostream>
#include <cstring>
#include "error.h"

using namespace flexiblesusy;
using namespace standard_model;

typedef CNMSSM_spectrum_generator<Semi_analytic> NMSSM;
typedef CNMSSM_effective_couplings bsm_eff_couplings;
typedef Standard_model_effective_couplings sm_eff_couplings;

constexpr double MA_unphysical = 5.;


common::pair error(common::parameters x, std::string what, std::string where) {
  /*
    Treat errors in a consistent way.
  */
  std::cerr << where << std::endl;
  std::cerr << what << std::endl;

  for (common::parameters::const_iterator p = x.begin(); p != x.end(); ++p) {
    std::cerr << *p << std::endl;
  }

  // Return back to sampler
  return {x, common::unphysical};
}

common::pair loglike(common::parameters x, double threshold) {
  // We will sum the log-likelihood piece by piece
  double loglike_ = 0.;

  CNMSSM_input_parameters input;

  input.LambdaInput = x[0];
  input.m12 = x[1];
  input.TanBeta = x[2];
  input.Azero = x[3];
  input.SignvS = x[4];
  Physical_input physical_input;  // extra non-SLHA physical input
  // constructed with default values may want to fill
  softsusy::QedQcd qedqcd;
  // arrange settings
  Spectrum_generator_settings settings;
  settings.set(Spectrum_generator_settings::precision, 1.0e-4);
  settings.set(Spectrum_generator_settings::calculate_sm_masses, 1.0e-4);
  NMSSM spectrum_generator;
  spectrum_generator.set_settings(settings);

  // Run FS spectrum generator
  try {
    spectrum_generator.run(qedqcd, input);
  } catch (flexiblesusy::Error& e) {
    return error(x, e.what(), "run");
  } catch (std::exception& e) {
    return error(x, e.what(), "run");
  }

  // Fetch model from FlexibleSUSY
  auto model = std::get<0>(spectrum_generator.get_models_slha());

  // Check that model is physical
  bool unphysical = spectrum_generator.get_problems().have_problem();
  if (unphysical) {
    return {x, common::unphysical};
  }

  // Higgs masses
  auto Mh = model.get_Mhh_pole_slha();
  higgs::real m_h{Mh(0), Mh(1), Mh(2)};

  // Higgs couplings from mixing angles
  auto ZH = model.get_ZH_pole_slha();
  const double vu = model.get_vu();
  const double vd = model.get_vd();
  const double tb = vu / vd;
  const double cb = std::cos(std::atan(tb));
  const double sb = std::sin(std::atan(tb));

  higgs::complex Cu{-1., -1., -1.};
  higgs::complex Cd{-1., -1., -1.};
  higgs::real CV{-1., -1., -1.};
  higgs::real Cgammagamma{-1., -1., -1.};
  higgs::real Cgg{-1., -1., -1.};

  // Pseudoscalar masses
  auto MA = model.get_MAh_pole_slha();
  auto MA_min = std::min(std::min(MA(0), MA(1)), MA(2));

  if (MA_min < MA_unphysical) {
    return {x, common::unphysical};
  }

  // Higgs effective couplings from FlexibleSUSY
  bsm_eff_couplings effective_couplings(model, qedqcd, physical_input);

  try {
    effective_couplings.calculate_effective_couplings();
  } catch (std::exception& e) {
    return error(x, e.what(), "calculate_effective_couplings");
  } catch (flexiblesusy::Error& e) {
    return error(x, e.what(), "calculate_effective_couplings");
  }

  for (int i = 0; i < higgs::knHiggs; i++) {
    // h_i uu / h_{SM}uu = ZH(i,1) / sin (beta)
    // c.f. MSSM/THDM huu /h_{SM}uu = ca/sb = ZH(0,1)/sb
    //                Huu /h_{SM}uu = sa /sb = ZH(1,1) / sb
    Cu[i] = ZH(i, 1) / sb;
    // h_idd / h_{SM}dd = ZH(i,0) / cos (beta)
    // c.f. MSSM/THDM hdd /h_{SM}dd = -sa/cb = ZH(0,0)/cb
    //                Hdd /h_{SM}dd = ca /sb = ZH(1,0) / cb
    Cd[i] = ZH(i, 0) / cb;

    // h_iVV / h_{SM}VV = sb ZH(i,1) + cb ZH(i,0)
    CV[i] = sb * ZH(i, 1) + cb * ZH(i, 0);

    // Don't calculate sm effective couplings if outside lilith range
    if (Mh(i) < higgs::m_h_min || Mh(i) > higgs::m_h_max) {
      continue;
    }

    // set up sm TODO(PA) can the next 4 lines be outside the loop?
    Standard_model sm;
    sm.set_thresholds(1);
    sm.set_pole_mass_loop_order(1);
    sm.set_precision(1e-5);
    physical_input.set(Physical_input::mh_pole, Mh(i));
    sm.set_physical_input(physical_input);

    // calculate drbar parameters and pole masses
    try {
      sm.initialise_from_input(qedqcd);
      sm.calculate_spectrum();
    } catch (std::exception& e) {
      return error(x, e.what(), "calculate_spectrum");
    } catch (flexiblesusy::Error& e) {
      return error(x, e.what(), "calculate_spectrum");
    }

    // should never be a problem when in range (min_mh, max_mh)
    bool unphysical = sm.get_problems().have_problem();
    if (unphysical) {
      return error(x, "unphysical point", "get_problems().have_problem()");
    }

    // calculate sm effective couplings with fs
    sm_eff_couplings sm_effective_couplings(sm, qedqcd, physical_input);
    try {
      sm_effective_couplings.calculate_effective_couplings();
    } catch (std::exception& e) {
      return error(x, e.what(), "sm_effective_couplings");
    } catch (flexiblesusy::Error& e) {
      return error(x, e.what(), "sm_effective_couplings");
    }

    Cgammagamma[i] = std::abs(effective_couplings.get_eff_CphhVPVP(i))
      / std::abs(sm_effective_couplings.get_eff_CphhVPVP());
    Cgg[i] = std::abs(effective_couplings.get_eff_CphhVGVG(i))
      / std::abs(sm_effective_couplings.get_eff_CphhVGVG());
  }

  // Higgs log-likelihood from Lilith
  double loglike_higgs = higgs::loglike(m_h, Cd, Cu, CV, Cgammagamma, Cgg);

  loglike_ += loglike_higgs;

  // TODO(AF) - use this pattern to abort the calculations if the loglike is no
  // good. The threshold argument contains the current lowest loglike in the
  // live points - points with worse loglike are anyway discarded by MultiNest.
  // if (loglike_ < threshold) {
  //   return {x, loglike_};
  // }

  // Higgs masses and couplings
  for (int i = 0; i < higgs::knHiggs; i += 1) {
    x.push_back(m_h[i]);
    x.push_back(Cu[i].real());
    x.push_back(Cu[i].imag());
    x.push_back(Cd[i].real());
    x.push_back(Cd[i].imag());
    x.push_back(CV[i]);
    x.push_back(Cgammagamma[i]);
    x.push_back(Cgg[i]);
  }

  // Phase-transition properties TODO

  return {x, loglike_};
}

common::parameters scan_prior(common::parameters x) {
  /*
    Priors for parameters in general scan.

    @param x - Unit hypercube
    @returns Physical parameters
  */
  x[0] = priors::flat(x[0], 0., 2. * M_PI);  // \lambda
  x[1] = priors::log(x[5], 1.e2, 3.e4);  // m12
  x[2] = priors::signed_log(x[2], 1.0, 1.e4);  // A_0
  x[3] = 1;
  return x;
}

common::parameters bm_prior(common::parameters x) {
  /*
    Priors for parameters in scan around benchmark point in 1407.4134.

    @param x - Unit hypercube
    @returns Physical parameters
  */
  x[0] = priors::flat(x[0], 0., 2. * M_PI);  // \lambda
  x[1] = 5000.;  //m12
  x[3] = 200.;  // A_0
  x[4] = 1;  // sign[vS]
  return x;
}

int main() {
  auto settings = Settings();
  settings.n_dims = 10;
  settings.n_pars = 100;  // TODO(AF) set properly
  Py_Initialize();
  nested::run(loglike, bm_prior, settings);
  Py_Finalize();
  return 0;
}
