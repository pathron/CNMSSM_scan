/*
  Example of using lilith with multinest algorithm.
*/

#include "multinest.hpp"
#include "higgs.hpp"


common::pair loglike(common::parameters x, double threshold) {
  const double a = x[0];
  const double b = x[1];
  double loglike_ = -(pow(a - 0.42, 2) + pow(b - 0.57, 2));

  higgs::real m_h {129, 125., 125.};
  higgs::complex Ct{1., 1., 1.};
  higgs::complex Cb{1., 1., 1.};
  higgs::complex Cc{1., 1., 1.};
  higgs::complex Ctau{1., 1., 1.};
  higgs::real CW{1., 1., 1.};
  higgs::real CZ{1., 1., 1.};
  higgs::real Cgammagamma{1., 1., 1.};
  higgs::real CZgamma{1., 1., 1.};
  higgs::real Cgg_prod{1., 1., 1.};
  higgs::real Cgg_decay{1., 1., 1.};
  higgs::real CVBF{1., 1., 1.};

  double loglike_higgs = higgs::loglike(m_h, Ct, Cb, Cc, Ctau, CW, CZ,
    Cgammagamma, CZgamma, Cgg_prod, Cgg_decay, CVBF);
  return {x, loglike_ + loglike_higgs};
}

common::parameters prior(common::parameters x) {
  return x;
}

int main() {
  auto settings = nested::Settings();
  settings.n_dims = 2;
  Py_Initialize();
  nested::run(loglike, prior, settings);
  Py_Finalize();
  return 0;
}
