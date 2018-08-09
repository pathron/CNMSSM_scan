// NB that callers should wrap programs that call this function with
// Py_Initialize();
// ...
// Py_Finalize();

#ifndef HIGGS_HPP
#define HIGGS_HPP

#include <Python.h>
#include <limits>
#include <array>
#include <complex>
#include <string>
#include "lilith.h"
#include <algorithm>

namespace higgs {

typedef std::array<std::complex<double>, 3> complex;
typedef std::array<double, 3> real;
constexpr int knHiggs = 3;
constexpr double m_h_min = 123.;
constexpr double m_h_max = 128.;

double loglike(real m_h,
               complex Ct,
               complex Cb,
               complex Cc,
               complex Ctau,
               real CW,
               real CZ,
               real Cgammagamma,
               real CZgamma,
               real Cgg_prod,
               real Cgg_decay,
               real CVBF,
               std::string precision = "BEST-QCD",
               bool ad_hoc_penalty = true
               ) {
  // Check there is a 125 GeV candidate

  int candidate = (m_h[0] > m_h_min && m_h[0] < m_h_max) ||
                  (m_h[1] > m_h_min && m_h[1] < m_h_max) ||
                  (m_h[2] > m_h_min && m_h[2] < m_h_max);

  if (!candidate) {
    if (ad_hoc_penalty) {
      constexpr double mu = 125.5;
      const double diff = std::min({std::abs(m_h[0] - mu), std::abs(m_h[1] - mu), std::abs(m_h[2] - mu)});
      const double penalty = -pow(diff, 2) * 1.e3;
      return penalty;
    } else {
      return -std::numeric_limits<double>::max();
    }
  }

  // Creating an object of the class Lilith with default ("") data

  char data[] = "";
  PyObject* lilithcalc = initialize_lilith(data);

  // Creating an XML input string for the complex coupling mode

  char XMLinputstring[6000] = "";
  char buffer[100];

  // Header

  snprintf(buffer, sizeof(buffer), "<?xml version=\"1.0\"?>\n");
  strcat(XMLinputstring, buffer);
  snprintf(buffer, sizeof(buffer), "<lilithinput>\n");
  strcat(XMLinputstring, buffer);

  // Setting mass and couplings of Higgs bosons

  for (int i = 0; i < knHiggs; i += 1) {
    // Higgs must be inside a window around 125 GeV

    int window_125 = m_h[i] > m_h_min && m_h[i] < m_h_max;

    if (!window_125) {
      continue;
    }

    int part = i + 1;

    snprintf(buffer, sizeof(buffer), "<reducedcouplings part=\"h%i\">\n", part);
    strcat(XMLinputstring, buffer);

    snprintf(buffer, sizeof(buffer), "<m_h>%f</m_h>\n", m_h[i]);
    strcat(XMLinputstring, buffer);

    snprintf(buffer, sizeof(buffer), "<C to=\"tt\" part=\"re\">%f</C>\n", Ct[i].real());
    strcat(XMLinputstring, buffer);
    snprintf(buffer, sizeof(buffer), "<C to=\"tt\" part=\"im\">%f</C>\n", Ct[i].imag());
    strcat(XMLinputstring, buffer);

    snprintf(buffer, sizeof(buffer), "<C to=\"bb\" part=\"re\">%f</C>\n", Cb[i].real());
    strcat(XMLinputstring, buffer);
    snprintf(buffer, sizeof(buffer), "<C to=\"bb\" part=\"im\">%f</C>\n", Cb[i].imag());
    strcat(XMLinputstring, buffer);

    snprintf(buffer, sizeof(buffer), "<C to=\"cc\" part=\"re\">%f</C>\n", Cc[i].real());
    strcat(XMLinputstring, buffer);
    snprintf(buffer, sizeof(buffer), "<C to=\"cc\" part=\"im\">%f</C>\n", Cc[i].imag());
    strcat(XMLinputstring, buffer);

    snprintf(buffer, sizeof(buffer), "<C to=\"tautau\" part=\"re\">%f</C>\n", Ctau[i].real());
    strcat(XMLinputstring, buffer);
    snprintf(buffer, sizeof(buffer), "<C to=\"tautau\" part=\"im\">%f</C>\n", Ctau[i].imag());
    strcat(XMLinputstring, buffer);

    // Couplings that would be calcualted if absent (indicated by negative coupling)

    if (CW[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"WW\">%f</C>\n", CW[i]);
      strcat(XMLinputstring, buffer);
    }
    if (CZ[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"ZZ\">%f</C>\n", CZ[i]);
      strcat(XMLinputstring, buffer);
    }
    if (Cgammagamma[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"gammagamma\">%f</C>\n", Cgammagamma[i]);
      strcat(XMLinputstring, buffer);
    }
    if (CZgamma[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"Zgamma\">%f</C>\n", CZgamma[i]);
      strcat(XMLinputstring, buffer);
    }
    if (Cgg_prod[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"gg\" for=\"prod\">%f</C>\n", Cgg_prod[i]);
      strcat(XMLinputstring, buffer);
    }
    if (Cgg_decay[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"gg\" for=\"decay\">%f</C>\n", Cgg_decay[i]);
      strcat(XMLinputstring, buffer);
    }
    if (CVBF[i] >= 0) {
      snprintf(buffer, sizeof(buffer), "<C to=\"VBF\">%f</C>\n", CVBF[i]);
      strcat(XMLinputstring, buffer);
    }

    snprintf(buffer, sizeof(buffer), "<precision>%s</precision>\n", precision.c_str());
    strcat(XMLinputstring, buffer);
    snprintf(buffer, sizeof(buffer), "</reducedcouplings>\n");
    strcat(XMLinputstring, buffer);
  }

  // Finish

  snprintf(buffer, sizeof(buffer), "</lilithinput>\n");
  strcat(XMLinputstring, buffer);

  // Reading user input XML string

  lilith_readuserinput(lilithcalc, XMLinputstring);

  // Calculating log-likelihood

  double loglike_ = -0.5 * lilith_computelikelihood(lilithcalc);
  return loglike_;
}

double loglike(real m_h,
               complex Cu,
               complex Cd,
               real CV,
               real Cgammagamma,
               real Cgg
               ) {
  real dummy{-1., -1., -1.};
  return loglike(m_h, Cu, Cd, Cu, Cd, CV, CV, Cgammagamma, dummy, Cgg, Cgg, dummy);
}

}  // namespace higgs

#endif  // HIGGS_HPP
