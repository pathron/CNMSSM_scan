#include "sample.hpp"
#include "higgs.hpp"
#include "common.hpp"
#include <cmath>
#include "math.h"

#include "config.h"
#include "CNMSSM_input_parameters.hpp"
#include "CNMSSM_effective_couplings.hpp"
#include "CNMSSM_observables.hpp"
#include "CNMSSM_slha_io.hpp"
#include "CNMSSM_spectrum_generator.hpp"
#include "CNMSSM_semi_analytic_spectrum_generator.hpp"
#include "CNMSSM_utilities.hpp"

#include "problems.hpp"
#include "spectrum_generator_problems.hpp"

#include "standard_model.hpp"
#include "standard_model_effective_couplings.hpp"

#include "command_line_options.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"
#include "wrappers.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include "error.h"

#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include <dirent.h>
#include <pthread.h>
#include <sys/types.h>
#include <cstdlib>
#include <iomanip>
#include <random>
//#include <mpi.h>


#include "../MicrOMEGAs/include/micromegas.h"
#include "../MicrOMEGAs/include/micromegas_aux.h"
#include "../MicrOMEGAs/CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"

using namespace flexiblesusy;
using namespace standard_model;

typedef CNMSSM_spectrum_generator<Semi_analytic> NMSSM; //What does the <Two_scale> mean? Why do we call it NMSSM?
//typedef THDMIISNMSSMBCsimple_spectrum_generator<Two_scale> NMSSM;
typedef CNMSSM_effective_couplings bsm_eff_couplings;
//typedef THDMIISNMSSMBCsimple_effective_couplings bsm_eff_couplings;
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


std::string random_string(std::string::size_type length)
{
    static auto& chrs = "0123456789"
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    thread_local static std::mt19937 rg{std::random_device{}()};
    thread_local static std::uniform_int_distribution<std::string::size_type> pick(0, sizeof(chrs) - 2);

    std::string s;

    s.reserve(length);

    while(length--)
        s += chrs[pick(rg)];

    return s;
}


common::pair loglike(common::parameters x, double threshold) {
  // We will sum the log-likelihood piece by piece
  double loglike_ = 0.;


CNMSSM_input_parameters input;
//THDMIISNMSSMBCsimple_input_parameters input;

input.m12 = x[0];
input.TanBeta = x[1];
input.Azero = x[2];
input.LambdaInput = x[3];
input.SignvS = x[4]; //fix sign

//std::cout << "-------------- START SCAN------------------------" << std::endl;

std::cout << "M12 = "  << input.m12 << " "
          << "TanBeta = "  << input.TanBeta << " "
          << "Azero = "  << input.Azero << " "
          << "Lambda = " << input.LambdaInput << " "
          << "SignvS = "  << input.SignvS << std::endl;



Physical_input physical_input;  // extra non-SLHA physical input
// constructed with default values may want to fill
softsusy::QedQcd qedqcd;
// arrange settings
Spectrum_generator_settings settings;
settings.set(Spectrum_generator_settings::precision, 1.0e-4);
settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
NMSSM spectrum_generator;
spectrum_generator.set_settings(settings);

//////////////////////////////////

  // std::cout.precision(20);
  //  /// SMINPUTS
  // std::cout << "Pole Mt  = " << qedqcd.displayPoleMt() << std::endl;
  // std::cout << "Pole Mtau  = " << qedqcd.displayPoleMtau() << std::endl;
  // std::cout << "Pole Mb  = " << qedqcd.displayPoleMb() << std::endl;
  // std::cout << "Pole MW  = " << qedqcd.displayPoleMW() << std::endl;
  // std::cout << "Pole MZ  = " << qedqcd.displayPoleMZ() << std::endl;
  // std::cout << "Pole Mmuon  = " << qedqcd.displayPoleMmuon() << std::endl;
  // std::cout << "Pole Mel  = " << qedqcd.displayPoleMel() << std::endl;
  //
  // std::cout << "alpha_s input  = " << qedqcd.displayAlphaSInput() << std::endl;
  // std::cout << "alpha_e input  = " << qedqcd.displayAlphaEmInput() << std::endl;
  // std::cout << "G_F  = " << qedqcd.displayFermiConstant() << std::endl;
  //
  // std::cout << "MbMb = " << qedqcd.displayMbMb() << std::endl;
  // std::cout << "McMc = " << qedqcd.displayMcMc() << std::endl;
  // std::cout << "Mu at 2 GeV = " << qedqcd.displayMu2GeV() << std::endl;
  // std::cout << "Md at 2 GeV = " << qedqcd.displayMd2GeV() << std::endl;
  // std::cout << "Ms at 2 GeV = " << qedqcd.displayMs2GeV() << std::endl;
  //
  //   /// FlexibleSUSY settings
  // std::cout << "beta loop order = "  << settings.get(settings.beta_loop_order) << std::endl;
  // std::cout << "threshold_corrections_loop_order = "  <<  settings.get(settings.threshold_corrections_loop_order) << std::endl;
  // std::cout << "top_pole_qcd_corrections = "  << settings.get(settings.top_pole_qcd_corrections) << std::endl;
  // std::cout << "individual threshold_corrections = "  << settings.get(settings.threshold_corrections) << std::endl;
  //
  // std::cout << "eft_matching_loop_order_up = "  << settings.get(settings.eft_matching_loop_order_up) << std::endl;
  //
  // std::cout << "eft_matching_loop_order_down = "  << settings.get(settings.eft_matching_loop_order_down) << std::endl;

///////////////////////////////////////

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

 std::cout << "SM_Higgs = " << Mh(0) << " " // m_h(0)
          << "m_h(1) = " << Mh(1) << " "
          << "m_h(2) = " << Mh(2) << " "
          //<< "m_h[0] = "  << m_h[0] << " "
          //<< "m_h[1] = "  << m_h[1] << " "
          //<< "m_h[2] = "  << m_h[2]
          << std::endl;

  //get singlet vev (vS) and mu_effective
  const double lambda_out = model.get_Lambdax();
  const double vS = model.get_vS();
  const double mu_eff = 1.0/sqrt(2.0)*lambda_out*vS;

  //get universal scalar mass squared, m0sq, from flexiblesusy output
  const double m0sq = model.get_extra_parameters()(0);
  const double m0 = sqrt(m0sq);

  //Check Universal Scalar Mass
  if (m0sq < 0) {
    return {x, common::unphysical};
  }

  //get squark masses from flexiblesusy output
  auto msu = model.get_MSu_pole_slha();
  const double msu1 = msu(0);
  const double msu6 = msu(1);


    //Print Squark Masses
  std::cout << "m0sq = " << m0sq << " "
            << "m0 = " << m0 << " "
            << "msu1 = " << msu1 << " "
            << "msu6 = " << msu6 << std::endl;
  std::cout << "vS =  " << vS << " "
            << "mu_eff = " << mu_eff << std::endl;



  std::cout << "-------------- Write SLHA File------------------------" << std::endl;


  std::string random = random_string(40);
  //std::cout << random << '\n';
  std::string slha_output_file = "slhafiles/LesHouches.out.CNMSSM";
  //std::cout << slha_output_file << '\n';
  slha_output_file.append(random);
  //std::cout << slha_output_file << '\n';
  char *slha_output_file_tmp = &slha_output_file[0];
  std::cout << slha_output_file_tmp << '\n';


  Spectrum_generator_problems problems;

  CNMSSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   CNMSSM_observables observables;
   if (settings.get(Spectrum_generator_settings::calculate_observables))
      observables = calculate_observables(model, qedqcd, physical_input, scales.pole_mass_scale);



  CNMSSM_slha_io slha_io;

  slha_io.set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
  slha_io.set_input(input);
  slha_io.set_print_imaginary_parts_of_majorana_mixings(settings.get(Spectrum_generator_settings::force_positive_masses)); //Doesn't look like it is needed.
  slha_io.set_spectrum(model);
  slha_io.set_extra(model, scales, observables); // Needed for extra values in HMIX
  slha_io.write_to(slha_output_file);




  std::cout << "-------------- START MicrOMEGAs------------------------" << std::endl;


    //%rd             |slhaRead("/home/lukeb/Honours/CNMSSM_scan/slhafiles/LesHouches.out.CNMSSM",0)

  	int err, i;
	   	char lspname[10], nlspname[10];
		double Xf=-1;
    double Omega = -1;
    std::cout << "Initial OMEGA^2 = " << Omega << std::endl;
    //slhaRead("slhafiles/LesHouches.out.CNMSSM", 0);
    slhaRead(slha_output_file_tmp, 0);
    double w;
		double cut = 0.01;		// cut-off for channel output
		int fast = 1;			/* 0 = best accuracy, 1 = "fast option" accuracy ~1% 	     */
 		double Beps = 1.E-5;  		/* Criteqrium for including co-annihilations (1 = no coann.) */
 	  VZdecay=2; VWdecay=2; cleanDecayTable(); //Changed to 2 from 0 after the scan with 36000 points
    ForceUG=1;
    err = sortOddParticles(lspname);
    Omega = darkOmega(&Xf,fast,Beps,&err);
    std::cout << "Final OMEGA^2 = " << Omega << std::endl;

    // Observales to calculate
    double sum = printChannels(Xf,cut,Beps,1,stdout);
    printf("Result of printChannels = %.8E\n", sum);

    double pA0[2];
    double pA5[2];
    double nA0[2];
    double nA5[2];

    double Nmass=0.939; /*nucleon mass*/









    //   double SCcoeff;
    //
    //
    // printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");
    // printf("         TREE LEVEL\n");
    //
    //     nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    //     printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    //     printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    //     printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]);
    //
    //
    // printf("         BOX DIAGRAMS\n");
    //
    //
    //     nucleonAmplitudes(CDM1,  pA0,pA5,nA0,nA5);
    //     printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    //     printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    //     printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]);
    //
    //   SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    //     printf("CDM-nucleon cross sections[pb]:\n");
    //     printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    //     printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);
    //
    //
    //     fprintf(omega,"201 %6.15lf #\n",SCcoeff*pA0[0]*pA0[0]);
    //     fprintf(omega,"202 %6.15lf #\n",3*SCcoeff*pA5[0]*pA5[0]);
    //     fprintf(omega,"203 %6.15lf #\n",SCcoeff*nA0[0]*nA0[0]);
    //     fprintf(omega,"204 %6.15lf # \n",3*SCcoeff*nA5[0]*nA5[0]);
    // }
    //
    // {
    //   double dNdE[300];
    //   double nEvents;
    //   double nEventsCut;
    //
    // printf("\n======== Direct Detection ========\n");
    //
    //
    //
    //
    //   nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);
    //   printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
    //   nEventsCut=cutRecoilResult(dNdE,10,50);
    //   printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);                                   ;
    //   fprintf(omega,"301 %6.6lf #\n",nEvents);
    //
    //   nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);
    //   printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
    //   nEventsCut=cutRecoilResult(dNdE,10,50);
    //   printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);
    //   fprintf(omega,"302 %6.6lf #\n",nEvents);
    //
    //   nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);
    //   printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
    //   nEventsCut=cutRecoilResult(dNdE,10,50);
    //   printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);
    //   fprintf(omega,"303 %6.6lf #\n",nEvents);
    //
    //   nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);
    //   printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
    //   nEventsCut=cutRecoilResult(dNdE,10,50);
    //   printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);
    //   fprintf(omega,"304 %6.6lf #\n",nEvents);
    //
    //
    // }
    //
    //
    //        fclose(channels);
    //        fclose(omega);
    //

// std::remove("slhafiles/LesHouches.out.CNMSSM");





  // Higgs couplings from mixing angles
  auto ZH = model.get_ZH_pole_slha();
  const double vu = model.get_vu();
  //const double vu = model.get_v2();
  const double vd = model.get_vd();
  //const double vd = model.get_v1();
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

  //MO Loglike contribution

  const double Omega_planck = 0.1198; //0.1200
  const double Omega_planck_uncertainty = 0.0012; //0.0012

  double loglike_MO = - pow((((Omega)-(Omega_planck)) / (Omega_planck_uncertainty)), 2);

  loglike_ += loglike_MO;


  // TODO(AF) - use this pattern to abort the calculations if the loglike is no
  // good. The threshold argument contains the current lowest loglike in the
  // live points - points with worse loglike are anyway discarded by MultiNest.
  // if (loglike_ < threshold) {
  //   return {x, loglike_};
  // }



  // Higgs masses and couplings
  for (int i = 0; i < higgs::knHiggs; i += 1) {
    x.push_back(m_h[i]);
    //x.push_back(Cu[i].real());
    //x.push_back(Cu[i].imag());
    //x.push_back(Cd[i].real());
    //x.push_back(Cd[i].imag());
    //x.push_back(CV[i]);
    //x.push_back(Cgammagamma[i]);
    //x.push_back(Cgg[i]);
  }

  x.push_back(m0sq);
  x.push_back(m0);
  x.push_back(msu1);
  x.push_back(msu6);
  x.push_back(vS);
  x.push_back(mu_eff);
  x.push_back(Omega);



  std::cout << "-------------- SUCCESSFUL POINT------------------------" << std::endl;

  return {x, loglike_};
}


common::parameters test_fixed_prior(common::parameters x) {

  x[0] = 133.33;  //m12
  x[1] = 10.;   //TanBeta
  x[2] = -300.;   //Azero
  x[3] = -0.05;  // Lambda
  x[4] = 1.;   //SignvS
  return x;
}

common::parameters oned_test_prior(common::parameters x) {

  x[0] = 133.33; //m12
  x[1] = 10.;  //tanbeta
  x[2] = -300; //AZero
  x[3] = priors::flat(x[3], -0.10, 0.01); //priors::flat(x[0], 0., 2. * M_PI); //lambda
  x[4] = 1.; //signvS
  return x;
}



common::parameters oned_prior(common::parameters x) {

  x[0] = 500.; //m12
  x[1] = 10.;  //tanbeta
  x[2] = 0; //AZero
  x[3] = priors::flat(x[3], -0.03, 0.03); //priors::flat(x[0], 0., 2. * M_PI); //lambda
  x[4] = 1.; //signvS
  return x;
}

common::parameters oned_2_prior(common::parameters x) {

  x[0] = 500.; //m12
  x[1] = 10.;  //tanbeta
  x[2] = 0; //AZero
  x[3] = priors::flat(x[3], -0.03, 0.03); //priors::flat(x[0], 0., 2. * M_PI); //lambda
  x[4] = 1.; //signvS
  return x;
}

common::parameters fourd_test_prior(common::parameters x) {

  x[0] = priors::log(x[0], 1., 2.e4); //m12
  x[1] = priors::flat(x[1], 2., 25.); //tanbeta
  x[2] = priors::log(x[2], -1.e3, -1.); //Azero
  x[3] = priors::flat(x[3], -0.08, -0.02); //priors::flat(x[0], 0., 2. * M_PI); //lambda
  x[4] = 1.; //signvS
  return x;
}


common::parameters fourd_prior(common::parameters x) {

  x[0] = priors::log(x[0], 1., 5000); //m12
  x[1] = priors::flat(x[1], 2., 50.); //tanbeta
  x[2] = priors::log(x[2], -5.e3, -1.); //Azero
  x[3] = priors::flat(x[3], -0.20, 0.10); //lambda
  x[4] = 1.; //signvs
  return x;
}


int main() {
  auto settings = Settings();
  settings.n_dims = 5;
  settings.n_pars = 100;  // TODO(AF) set properly
  Py_Initialize();
  nested::run(loglike, fourd_prior, settings);
  Py_Finalize();
  return 0;
}
