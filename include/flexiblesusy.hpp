#ifndef FLEXIBLESUSY_HPP
#define FLEXIBLESUSY_HPP

#include "config.h"
#include "THDMIISNMSSMBCsimple_input_parameters.hpp"
#include "THDMIISNMSSMBCsimple_observables.hpp"
#include "THDMIISNMSSMBCsimple_slha_io.hpp"
#include "THDMIISNMSSMBCsimple_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "THDMIISNMSSMBCsimple_two_scale_spectrum_generator.hpp"
#endif

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"
#include "wrappers.hpp"
#include <iostream>
#include <cstring>

#include "common.hpp"

namespace nmssm {

using namespace flexiblesusy;

typedef THDMIISNMSSMBCsimple_spectrum_generator<Two_scale> NMSSM;

NMSSM spectrum(common::parameters x) {

  THDMIISNMSSMBCsimple_input_parameters input;

  input.lambdaNMSSM = x[0];
  input.kappaNMSSM = x[1];
  input.AlambdaNMSSM = x[2];
  input.AkappaNMSSM = x[3];
  input.AtopNMSSM = x[4];
  input.mstopL = x[5];
  input.mstopR = x[6];
  input.vSIN = x[7];
  input.TanBeta = x[8];
  input.MEWSB = x[9];
  
  Physical_input physical_input; // extra non-SLHA physical input
  // constructed with default values may want to fill
  softsusy::QedQcd qedqcd;
  // arrange settings
  Spectrum_generator_settings settings;
  settings.set(Spectrum_generator_settings::precision, 1.0e-4);
  NMSSM spectrum_generator;
  spectrum_generator.set_settings(settings);
  // run FS spectrum generator
  spectrum_generator.run(qedqcd, input);

	return spectrum_generator;
}

}  // namespace nmssm

#endif  // FLEXIBLESUSY_HPP
