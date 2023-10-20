#include "TauAnalysis/PolarimetricVectorTau2a1/interface/PolarimetricVectorTau2a1.h"

#include "TMath.h"   // TMath::Pi()
#include "TString.h" // Form()

#include <assert.h>  // assert()
#include <cmath>     // pow(), std::cos(), std::fabs(), std::sin()
#include <iostream>  // std::cerr, std::cout, std::ostream

const size_t numResonances = 7;

namespace
{
  /**
   * @brief Construct complex number, given its modulus and phase
   * @param modulus
   * @param phase
   * @return complex number
   */
  PolarimetricVectorTau2a1::cdouble
  get_beta(double modulus, double phase)
  {
    double re = modulus*std::cos(phase);
    double im = modulus*std::sin(phase);
    return PolarimetricVectorTau2a1::cdouble(re, im);
  }

  /**
   * @brief Return element of metric tensor g at index {mu,nu}
   * @param mu first index
   * @param nu second index
   * @return g_{mu,nu}
   */
  double
  get_g(size_t mu, size_t nu)
  {
    if      ( (mu == 0 && nu == 0) || (mu == 1 && nu == 1) || (mu == 2 && nu == 2) ) return -1.;
    else if (  mu == 3 && nu == 3                                                  ) return +1.;
    else                                                                             return  0.;
  }

  /**
   * @brief Auxiliary functions to print-out complex number, real and complex four-vectors 
   *       (to be used only for debugging purposes)
   * @param os reference to std::cout
   * @param val complex number or four-vector to be printed-out
   * @return os
   */
  std::ostream&
  operator<<(std::ostream& os, const PolarimetricVectorTau2a1::cdouble& val)
  {
    os << val.real();
    os << ( val.imag() >= 0. ? " + " : " - " );
    os << std::fabs(val.imag()) << "i";
    return os;
  }

  void
  print(const std::string& label, const PolarimetricVectorTau2a1::cdouble& val)
  {
    std::cout << label << ": " << val << "\n";
  }

  std::ostream&
  operator<<(std::ostream& os, const PolarimetricVectorTau2a1::cLorentzVector& val)
  {
    os << "E = " << val(3) << ", Px = " << val(0) << ", Py = " << val(1) << ", Pz = " << val(2);
    return os;
  }

  void
  print(const std::string& label, const PolarimetricVectorTau2a1::cLorentzVector& val)
  {
    std::cout << label << ": " << val << "\n";
  }

  std::ostream&
  operator<<(std::ostream& os, const PolarimetricVectorTau2a1::LorentzVector& val)
  {
    os << "E = " << val.energy() << ", Px = " << val.px() << ", Py = " << val.py() << ", Px = " << val.pz();
    return os;
  }

  void
  print(const std::string& label, const PolarimetricVectorTau2a1::LorentzVector& val)
  {
    std::cout << label << ": " << val << "\n";
  }
}

PolarimetricVectorTau2a1::PolarimetricVectorTau2a1(int verbosity)
  : verbosity_(verbosity)
{
  verbosity_ = -1; // CV: disable all debug output !!

  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::PolarimetricVectorTau2a1>:\n";
  }

  // define mass of tau lepton, charged and neutral pion;
  // values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01 (PDG)
  m_tau_          = 1.7769;   // [GeV]
  m_chargedPi_    = 0.139570; // [GeV]
  m_neutralPi_    = 0.134977; // [GeV]

  // define mass and width of a1(1260) meson;
  // values are taken from column "nominal fit" in Table VI of Phys.Rev.D 61 (2000) 012002
  //m0_a1_          = 1.331;    // [GeV]
  //Gamma0_a1_      = 0.814;    // [GeV]
  // values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01 (PDG)
  //m0_a1_          = 1.230;    // [GeV]
  //Gamma0_a1_      = 0.420;    // [GeV]
  // values are taken from TAUOLA code
  //m0_a1_          = 1.251;    // [GeV] <-- CV: taken from Setup function
  //Gamma0_a1_      = 0.599;    // [GeV] <-- CV: taken from Setup function
  m0_a1_          = 1.275;    // [GeV] <-- CV: taken from FA1A1P function
  Gamma0_a1_      = 0.700;    // [GeV] <-- CV: taken from FA1A1P function

  // define parameters specifying "running" of a1 width;
  // the values of Gamma_a1 as function of s have been taken from Fig. 9 (b) of Phys.Rev.D 61 (2000) 012002
  Gamma_a1_vs_s_ = { 
    { 0.00, 0.000 }, { 0.20, 0.000 }, { 0.40, 0.000 }, { 0.50, 0.005 }, { 0.60, 0.020 },
    { 0.65, 0.040 }, { 0.70, 0.055 }, { 0.75, 0.075 }, { 0.80, 0.110 }, { 0.85, 0.160 },
    { 0.90, 0.205 }, { 0.95, 0.250 }, { 1.00, 0.295 }, { 1.05, 0.340 }, { 1.10, 0.375 },
    { 1.15, 0.410 }, { 1.20, 0.450 }, { 1.25, 0.475 }, { 1.30, 0.515 }, { 1.35, 0.555 },
    { 1.40, 0.595 }, { 1.45, 0.630 }, { 1.50, 0.660 }, { 1.55, 0.690 }, { 1.60, 0.720 },
    { 1.65, 0.750 }, { 1.70, 0.780 }, { 1.75, 0.815 }, { 1.80, 0.845 }, { 1.85, 0.875 },
    { 1.90, 0.905 }, { 1.93, 0.930 }, { 1.95, 0.980 }, { 2.00, 1.060 }, { 2.05, 1.125 },
    { 2.10, 1.185 }, { 2.15, 1.245 }, { 2.20, 1.300 }, { 2.25, 1.355 }, { 2.30, 1.415 },
    { 2.35, 1.470 }, { 2.37, 1.485 }, { 2.40, 1.520 }, { 2.45, 1.575 }, { 2.50, 1.640 },
    { 2.55, 1.705 }, { 2.60, 1.765 }, { 2.65, 1.835 }, { 2.70, 1.900 }, { 2.75, 1.970 },
    { 2.80, 2.050 }, { 2.85, 2.130 }, { 2.90, 2.205 }, { 2.95, 2.285 }, { 3.00, 2.380 },
    { 3.05, 2.470 }, { 3.10, 2.570 }, { 3.15, 2.690 } 
  };

  // define masses and widths of intermediate rho(770), rho(1450), f2(1270), sigma, f0(1370) resonances;
  // values are taken from Table I of Phys.Rev.D 61 (2000) 012002
  m0_rho770_      = 0.774;    // [GeV]
  Gamma0_rho770_  = 0.149;    // [GeV]
  m0_rho1450_     = 1.370;    // [GeV]
  Gamma0_rho1450_ = 0.386;    // [GeV]
  m0_f2_          = 1.275;    // [GeV]
  Gamma0_f2_      = 0.185;    // [GeV]
  m0_sigma_       = 0.860;    // [GeV]
  Gamma0_sigma_   = 0.880;    // [GeV]
  m0_f0_          = 1.186;    // [GeV]
  Gamma0_f0_      = 0.350;    // [GeV]

  // define coefficients specifying the contribution of meson resonances to the hadronic current J;
  // values are taken from Table III of Phys.Rev.D 61 (2000) 012002
  std::vector<double> beta_moduli = {  1.00,  0.12,  0.37,  0.87,  0.71,  2.10,  0.77 };
  std::vector<double> beta_phases = {  0.00,  0.99, -0.15,  0.53,  0.56,  0.23, -0.54 };
  assert(beta_moduli.size() == numResonances && beta_phases.size() == numResonances);
  for ( size_t idxResonance = 0; idxResonance < numResonances; ++idxResonance )
  {
    cdouble beta = get_beta(beta_moduli[idxResonance], beta_phases[idxResonance]*TMath::Pi());
    if ( verbosity_ >= 2 )
    {
      std::cout << "beta[" << idxResonance << "] = " << beta << "\n";
    }
    beta_.push_back(beta);
  }

  for ( size_t mu = 0; mu < 4; ++mu )
  {
    for ( size_t nu = 0; nu < 4; ++nu )
    {
      g_(mu, nu) = get_g(mu, nu);
    }
  }
  if ( verbosity_ >= 2 )
  {
    std::cout << "g:\n";
    std::cout << g_ << "\n";
  }
}

PolarimetricVectorTau2a1::~PolarimetricVectorTau2a1()
{}

namespace
{
  /**
   * @brief Convert four-vector of type cLorentzVector to type LorentzVector. 
   * Note: The imaginary part of the four-vector given as function argument is discarded in the conversion.
   * @param p four-vector
   * @return LorentzVector(p)
   */
  PolarimetricVectorTau2a1::LorentzVector
  convert_to_LorentzVector(const PolarimetricVectorTau2a1::cLorentzVector& p)
  {
    PolarimetricVectorTau2a1::LorentzVector retVal(p(0).real(), p(1).real(), p(2).real(), p(3).real());
    return retVal;
  } 
}

PolarimetricVectorTau2a1::Vector
PolarimetricVectorTau2a1::operator()(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                                     int charge,
                                     DecayChannel decayChannel) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::operator()>:\n";
  }

  LorentzVector P(0.,0.,0.,m_tau_);
  if ( verbosity_ >= 2 )
  {
    print("P", P);
  }

  LorentzVector N = P - (p1 + p2 + p3);
  if ( verbosity_ >= 2 )
  {
    print("N", N);
  }

  cLorentzVector J = comp_J(p1, p2, p3, decayChannel);
  if ( verbosity_ >= 2 )
  {
    comp_J2(p1, p2, p3, decayChannel);
  }
  // CV: set J to values computed by Vladimir's code in order to check Pi and Pi5 computation
  //J(0) = std::complex(-0.927262, -1.448420);
  //J(1) = std::complex( 9.008140,  2.929270);
  //J(2) = std::complex(-8.217810, -2.278660);
  //J(3) = std::complex( 0.640295, -0.116947);

  LorentzVector Pi = comp_Pi(J, N);
  LorentzVector Pi5 = comp_Pi5(J, N, charge);

  const double gammaVA = 1.; // CV: Standard Model value, cf. Section 3.2 in Comput.Phys.Commun. 64 (1991) 275
  double omega = P.Dot(Pi - gammaVA*Pi5);

  const double M = m_tau_;
  LorentzVector H = (1./(omega*M))*(pow(M, 2)*(Pi5 - gammaVA*Pi) - P.Dot(Pi5 - gammaVA*Pi)*P);
  if ( verbosity_ >= 2 )
  {
    print("H", H);
  }

  Vector retVal = H.Vect().unit();
  if ( verbosity_ >= 2 )
  {
    std::cout << "h: Px = " << retVal.x() << ", Py = " << retVal.y() << ", Pz = " << retVal.z() << "\n";
  }
  return retVal;
}

namespace
{
  /**
   * @brief Return complex conjugate of given four-vector
   * @param p given four-vector
   * @return p*
   */
  PolarimetricVectorTau2a1::cLorentzVector
  star(const PolarimetricVectorTau2a1::cLorentzVector& p)
  {
    PolarimetricVectorTau2a1::cLorentzVector retVal;
    for ( size_t mu = 0; mu < 4; ++mu )
    {
      retVal(mu) = std::complex(p(mu).real(), -p(mu).imag());
    }
    return retVal;
  }

  /**
   * @brief Return certain component of given four-vector
   * @param p given four-vector
   * @param mu component to be returned [1]
   *
   *   [1] use 0 for px
   *       use 1 for py
   *       use 2 for pz
   *       use 3 for energy
   * 
   * @return p[mu]
   */
  double
  get_component(const PolarimetricVectorTau2a1::LorentzVector& p, size_t mu)
  {
         if ( mu == 0 ) return p.px();
    else if ( mu == 1 ) return p.py();
    else if ( mu == 2 ) return p.pz();
    else if ( mu == 3 ) return p.energy();
    else assert(0);
  }

  /**
   * @brief Convert four-vector of type LorentzVector to type cLorentzVector
   * @param p four-vector
   * @return cLorentzVector(p)
   */
  PolarimetricVectorTau2a1::cLorentzVector
  convert_to_cLorentzVector(const PolarimetricVectorTau2a1::LorentzVector& p)
  {
    PolarimetricVectorTau2a1::cLorentzVector retVal;
    for ( size_t mu = 0; mu < 4; ++mu )
    {
      retVal(mu) = get_component(p, mu);
    }
    return retVal;
  } 
}

PolarimetricVectorTau2a1::LorentzVector
PolarimetricVectorTau2a1::comp_Pi(const PolarimetricVectorTau2a1::cLorentzVector& J, const LorentzVector& N) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::comp_Pi>:\n";
  }

  cLorentzVector cN = convert_to_cLorentzVector(N);
  cLorentzVector Jstar = star(J);

  cdouble JstarXN = ROOT::Math::Dot(Jstar, g_*cN);
  cdouble JXN = ROOT::Math::Dot(J, g_*cN);
  cdouble JstarXJ = ROOT::Math::Dot(Jstar, g_*J);
  LorentzVector retVal = convert_to_LorentzVector(2.*(JstarXN*J + JXN*Jstar - JstarXJ*cN));
  if ( verbosity_ >= 2 )
  {
    print("Pi", retVal);
  }
  
  return retVal;
}

namespace
{
  /**
   * @brief Return sign of integer value given as function argument
   * @param x function argument
   * @return sign(x)
   */
  int
  sgn(int x)
  {
    if ( x > 0 ) return +1;
    if ( x < 0 ) return -1;
    return 0;
  }

  /**
   * @brief Compute Levi-Civita symbol epsilon^{mu,nu,rho,sigma},
   *        cf. https://en.wikipedia.org/wiki/Levi-Civita_symbol
   *       (Section "Levi-Civita tensors")
   * @params mu, nu,rho, sigma indices of Levi-Civita symbol
   * @return epsilon^{mu,nu,rho,sigma}
   */
  double
  get_epsilon(unsigned mu, unsigned nu, unsigned rho, unsigned sigma, int verbosity)
  {
    //if ( verbosity >= 3 )
    //{
    //  std::cout << "<get_epsilon>:\n";
    //}
    // CV: formula for computation of four-dimensional Levi-Civita symbol taken from
    //       https://en.wikipedia.org/wiki/Levi-Civita_symbol
    //    (Section "Definition" -> "Generalization to n dimensions")
    int a[4];
    a[0] = mu;
    a[1] = nu;
    a[2] = rho;
    a[3] = sigma;
    double epsilon = 1.;
    for ( int i = 0; i < 4; ++i )
    {
      for ( int j = i + 1; j < 4; ++j )
      {
        epsilon *= sgn(a[j] - a[i]); 
      }
    }
    if ( verbosity >= 3 )
    {
      std::cout << "epsilon(mu = " << mu << ", nu = " << nu << ", rho = " << rho << ", sigma = " << sigma << ") = " << epsilon << "\n";
    }
    return epsilon;
  }
}

PolarimetricVectorTau2a1::LorentzVector
PolarimetricVectorTau2a1::comp_Pi5(const cLorentzVector& J, const LorentzVector& N, int charge) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::comp_Pi5>:\n";
  }

  cLorentzVector cN = convert_to_cLorentzVector(N);
  cLorentzVector Jstar = star(J);

  cLorentzVector vProd;
  for ( size_t mu = 0; mu < 4; ++mu )
  {
    for ( size_t nu = 0; nu < 4; ++nu )
    {
      for ( size_t rho = 0; rho < 4; ++rho )
      {
        for ( size_t sigma = 0; sigma < 4; ++sigma )
        {
          double epsilon = get_epsilon(mu, nu, rho, sigma, verbosity_);
          vProd(mu) += epsilon*Jstar(nu)*J(rho)*cN(sigma);
        }
      }
    }
  }

  // CV: multiply with metric tensor in order to transform Pi5^{mu} into Pi5_{mu},
  //     as required to insert Pi5^{mu} given by Eq. (3.15) into Eq. (3.14) of Comput.Phys.Commun. 64 (1991) 275
  vProd = g_*vProd;

  //double sign = 0.;
  //if      ( charge == +1 ) sign = -1.;
  //else if ( charge == -1 ) sign = +1.;
  //else assert(0);

  //LorentzVector retVal(sign*2.*sum(0).imag(), sign*2.*sum(1).imag(), sign*2.*sum(2).imag(), -sign*2.*sum(3).imag());
  LorentzVector retVal(2.*vProd(0).imag(), 2.*vProd(1).imag(), 2.*vProd(2).imag(), 2.*vProd(3).imag());
  if ( verbosity_ >= 2 )
  {
    print("Pi5", retVal);
  }
  
  return retVal;
}

namespace
{
  /**
   * @brief Compute "decay momentum",
   *        given by bottom line of Eq. (A6) in Phys.Rev.D 61 (2000) 012002
   * @param si (mass of two-pion system that forms resonance)^2
   * @param mj mass of first  pion that forms resonance [1]
   * @param mk mass of second pion that forms resonance [1]
   *
   *  [1] the ordering of the two pions does not matter
   *
   * @return k'_i
   */
  double
  kdash(double si, double mj, double mk)
  {
    double retVal = std::sqrt((si - pow(mj + mk, 2))*(si - pow(mj - mk, 2)))/(2.*std::sqrt(si));
    return retVal;
  }

  /**
   * @brief Compute "running width" of intermediate rho(770), rho(1450), f2(1270), sigma, and f0(1370) resonances,
   *        given by bottom line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
   * @param m0 nominal mass of resonance
   * @param Gamma0 nominal width of resonance
   * @param si (mass of two-pion system that forms resonance)^2
   * @param mj mass of first  pion that forms resonance [1]
   * @param mk mass of second pion that forms resonance [1]
   * @param L angular momentum of resonance (s-wave = 0, p-wave=1, d-wave=2)
   *
   *  [1] the ordering of the two pions does not matter
   *
   * @return Gamma^{Y,L}(s_i)
   */
  double
  Gamma(double m0, double Gamma0, double si, double mj, double mk, unsigned L, int verbosity = -1)
  {
    if ( verbosity >= 2 )
    {
      std::cout << "<Gamma>:\n";
    }

    double kdashi = kdash(si, mj, mk);
    double kdash0 = kdash(pow(m0, 2), mj, mk);
    if ( verbosity >= 2 )
    {
      std::cout << "kdashi = " << kdashi << "\n";
      std::cout << "kdash0 = " << kdash0 << "\n";
    }

    double retVal = Gamma0*pow(kdashi/kdash0, 2*L + 1)*m0/std::sqrt(si);
    return retVal;
  }

  /**
   * @brief Compute Breit-Wigner function of intermediate rho(770), rho(1450), f2(1270), sigma, and f0(1370) resonances,
   *        given by top line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
   * @param m0 nominal mass of resonance
   * @param Gamma0 nominal width of resonance
   * @param si (mass of two-pion system that forms resonance)^2
   * @param mj mass of first  pion that forms resonance [1]
   * @param mk mass of second pion that forms resonance [1]
   * @param L angular momentum of resonance (s-wave = 0, p-wave=1, d-wave=2)
   *
   *  [1] the ordering of the two pions does not matter
   *
   * @return B^{L}_{Y}(s_i)
   */
  PolarimetricVectorTau2a1::cdouble
  BreitWigner(double m0, double Gamma0, double si, double mj, double mk, unsigned L, int verbosity = -1)
  {
    if ( verbosity >= 2 )
    {
      std::cout << "<BreitWigner>:\n";
      std::cout << "m0 = " << m0 << "\n";
      std::cout << "Gamma0 = " << Gamma0 << "\n";
      std::cout << "si = " << si << "\n";
      std::cout << "mj = " << mj << "\n";
      std::cout << "mk = " << mk << "\n";
      std::cout << "L = " << L << "\n";
    }

    double num = -pow(m0, 2);
    PolarimetricVectorTau2a1::cdouble denom(si - pow(m0, 2), m0*Gamma(m0, Gamma0, si, mj, mk, L, verbosity));
    if ( verbosity >= 2 )
    {
      std::cout << "num = " << num << "\n";
      std::cout << "denom = " << denom << "\n";
    }

    PolarimetricVectorTau2a1::cdouble retVal = PolarimetricVectorTau2a1::cdouble(num, 0.)/denom;
    if ( verbosity >= 2 )
    {
      print("B", retVal);
    }
    return retVal;
  }
}

PolarimetricVectorTau2a1::cLorentzVector
PolarimetricVectorTau2a1::comp_J(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                                 DecayChannel decayChannel) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::comp_J>:\n";
  }

  LorentzVector q1 = p2 - p3;
  LorentzVector q2 = p3 - p1;
  LorentzVector q3 = p1 - p2;

  LorentzVector h1 = p2 + p3;
  LorentzVector Q1 = h1 - p1;
  double s1 = h1.mag2();
  LorentzVector h2 = p1 + p3;
  LorentzVector Q2 = h2 - p2;
  double s2 = h2.mag2();
  LorentzVector h3 = p1 + p2;
  LorentzVector Q3 = h3 - p3;
  double s3 = h3.mag2();
  if ( verbosity_ >= 2 )
  {
    std::cout << "s1 = " << s1 << "\n";
    std::cout << "s2 = " << s2 << "\n";
    std::cout << "s3 = " << s3 << "\n";
  }

  double m1, m2, m3;
  if ( decayChannel == k3ChargedPi )
  {
    m1 = m_chargedPi_;
    m2 = m_chargedPi_;
    m3 = m_chargedPi_;
  }
  else if ( decayChannel == kChargedPi2NeutralPi )
  {
    m1 = m_neutralPi_;
    m2 = m_neutralPi_;
    m3 = m_chargedPi_;
  }
  else
  {
    std::cerr << "Error in <PolarimetricVectorTau2a1::J>: Invalid parameter 'decayChannel' = " << decayChannel << " !!\n";
    assert(0);
  }
  
  LorentzVector a = p1 + p2 + p3;
  double s = a.mag2();
  if ( verbosity_ >= 2 )
  {
    std::cout << "s = " << s << "\n";
  }

  cTensor T;
  for ( size_t mu = 0; mu < 4; ++mu )
  {
    for ( size_t nu = 0; nu < 4; ++nu )
    {
      // CV: all vectors without explicite mu indices are assumed to have the mu index as subscript,
      //     hence multiplication with the product of metric tensors g^{mu,mu}*g^{nu,nu} is required to transform
      //     the expression get_component(a, mu)*get_component(a, nu) into a^{mu}*a^{nu}
      T(mu,nu) = g_(mu,nu) - g_(mu,mu)*g_(nu,nu)*get_component(a, mu)*get_component(a, nu)/s;
    }
  }
  if ( verbosity_ >= 2 )
  {
    std::cout << "T:\n";
    std::cout << T << "\n";
  }
  
  // compute amplitudes for individual resonances according to Eq. (A3) in Phys.Rev.D 61 (2000) 012002
  // Note: all the factors F_R in Eq. (A3) are equal to one, as the nominal fit assumes the a1 size parameter R to be zero
  std::vector<cLorentzVector> j(numResonances);
  cLorentzVector cq1 = convert_to_cLorentzVector(q1);
  cLorentzVector cq2 = convert_to_cLorentzVector(q2);
  j[0] = T*(BreitWigner(m0_rho770_,  Gamma0_rho770_,  s1, m2, m3, 1)*cq1 - BreitWigner(m0_rho770_,  Gamma0_rho770_,  s2, m1, m3, 1)*cq2);
  j[1] = T*(BreitWigner(m0_rho1450_, Gamma0_rho1450_, s1, m2, m3, 1)*cq1 - BreitWigner(m0_rho1450_, Gamma0_rho1450_, s2, m1, m3, 1)*cq2);
  double aXq1 = a.Dot(q1);
  cLorentzVector cQ1 = convert_to_cLorentzVector(Q1);
  double aXq2 = a.Dot(q2);
  cLorentzVector cQ2 = convert_to_cLorentzVector(Q2);
  j[2] = T*(aXq1*BreitWigner(m0_rho770_,  Gamma0_rho770_,  s1, m2, m3, 1)*cQ1 - aXq2*BreitWigner(m0_rho770_,  Gamma0_rho770_,  s2, m1, m3, 1)*cQ2);
  j[3] = T*(aXq1*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s1, m2, m3, 1)*cQ1 - aXq2*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s2, m1, m3, 1)*cQ2);
  double aXq3 = a.Dot(q3);
  cLorentzVector cq3 = convert_to_cLorentzVector(q3);
  double q3Xq3 = q3.mag2();
  cLorentzVector ca = convert_to_cLorentzVector(a);
  double h3Xa = h3.Dot(a);
  cLorentzVector ch3 = convert_to_cLorentzVector(h3);
  j[4] = T*(BreitWigner(m0_f2_, Gamma0_f2_, s3, m1, m2, 2)*(aXq3*cq3 - (q3Xq3/3.)*(ca - (h3Xa/s3)*ch3)));
  cLorentzVector cQ3 = convert_to_cLorentzVector(Q3);
  j[5] = T*(BreitWigner(m0_sigma_, Gamma0_sigma_, s3, m1, m2, 0)*cQ3);
  j[6] = T*(BreitWigner(m0_f0_, Gamma0_f0_, s3, m1, m2, 0)*cQ3);
  if ( verbosity_ >= 2 )
  {
    print("B_rho^P(s1)", BreitWigner(m0_rho770_, Gamma0_rho770_, s1, m2, m3, 1, verbosity_));
    print("B_rho'^P(s1)", BreitWigner(m0_rho1450_, Gamma0_rho1450_, s1, m2, m3, 1, verbosity_));
    print("B_rho^P(s2)", BreitWigner(m0_rho770_, Gamma0_rho770_, s2, m1, m3, 1, verbosity_));
    print("B_rho'^P(s2)", BreitWigner(m0_rho1450_, Gamma0_rho1450_, s2, m1, m3, 1, verbosity_));
    print("B_f2^D(s3)", BreitWigner(m0_f2_, Gamma0_f2_, s3, m1, m2, 2, verbosity_));
    print("B_sigma^S(s3)", BreitWigner(m0_sigma_, Gamma0_sigma_, s3, m1, m2, 0, verbosity_));
    print("B_f0^S(s3)", BreitWigner(m0_f0_, Gamma0_f0_, s3, m1, m2, 0, verbosity_));
    for ( size_t idx = 0; idx < numResonances; ++idx )
    {
      print(Form("j[%lu]", idx), j[idx]);
    }
  }

  cLorentzVector retVal;
  for ( size_t idx = 0; idx < numResonances; ++idx )
  {
    retVal += beta_[idx]*j[idx];
  }
  retVal *= BreitWigner_a1(s);
  // CV: multiply with metric tensor in order to transform J^{mu} into J_{mu},
  //     as required to insert J^{mu} given by Eq. (A2) of Phys.Rev.D 61 (2000) 012002 into Eq. (3.15) of Comput.Phys.Commun. 64 (1991) 275
  retVal = g_*retVal;
  if ( verbosity_ >= 2 )
  {
    print("J", retVal);
  }

  return retVal;
}

PolarimetricVectorTau2a1::cLorentzVector
PolarimetricVectorTau2a1::comp_J2(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                                  DecayChannel decayChannel) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::comp_J2>:\n";
  }

  LorentzVector q1 = p2 - p3;
  LorentzVector q2 = p3 - p1;
  LorentzVector q3 = p1 - p2;
  if ( verbosity_ >= 2 )
  {
    print("q1", q1);
    print("q2", q2);
    print("q3", q3);
  }

  LorentzVector a = p1 + p2 + p3;
  double s = a.mag2();
  if ( verbosity_ >= 2 )
  {
    print("a", a);
    std::cout << "s = " << s << "\n";
  }

  cTensor T;
  for ( size_t mu = 0; mu < 4; ++mu )
  {
    for ( size_t nu = 0; nu < 4; ++nu )
    {
      T(mu,nu) = g_(mu,nu) - g_(mu,mu)*g_(nu,nu)*get_component(a, mu)*get_component(a, nu)/s;
      //double delta = ( mu == nu ) ? 1. : 0.;
      //T(mu,nu) = delta - get_component(a, mu)*get_component(a, nu)/s;
    }
  }
  if ( verbosity_ >= 2 )
  {
    std::cout << "T:\n";
    std::cout << T << "\n";
  }

  const double c1 =  1.;
  const double c2 = -1.;
  const double c3 =  1.;
  if ( verbosity_ >= 2 )
  {
    std::cout << "c1 = " << c1 << "\n";
    std::cout << "c2 = " << c2 << "\n";
    std::cout << "c3 = " << c3 << "\n";
  }
  
  cdouble F1 = comp_Fi(1, p1, p2, p3, decayChannel);
  cdouble F2 = comp_Fi(2, p1, p2, p3, decayChannel);
  cdouble F3 = comp_Fi(3, p1, p2, p3, decayChannel);
  if ( verbosity_ >= 2 )
  {
    print("F1", F1);
    print("F2", F2);
    print("F3", F3);
  }

  LorentzVector vec1 = q1 - (a.Dot(q1)/s)*a;
  LorentzVector vec2 = q2 - (a.Dot(q2)/s)*a;
  LorentzVector vec3 = q3 - (a.Dot(q3)/s)*a;
  if ( verbosity_ >= 2 )
  {
    print("vec1", vec1);
    print("vec2", vec2);
    print("vec3", vec3);
  }

  //cdouble F1 = comp_Fi(1, p1, p2, p3, decayChannel);
  //cdouble F2 = comp_Fi(2, p1, p2, p3, decayChannel);
  //cdouble F3 = comp_Fi(3, p1, p2, p3, decayChannel);
  // CV: multiply with metric tensor in order to transform J^{mu} into J_{mu},
  //     as required to insert J^{mu} given by Eq. (A2) of Phys.Rev.D 61 (2000) 012002 into Eq. (3.15) of Comput.Phys.Commun. 64 (1991) 275
  cLorentzVector retValTMP = g_*T*(c1*F1*convert_to_cLorentzVector(q1) + c2*F2*convert_to_cLorentzVector(q2) + c3*F3*convert_to_cLorentzVector(q3));
  if ( verbosity_ >= 2 )
  {    
    print("J2@1", retValTMP);
  }
  cLorentzVector retVal = c1*F1*convert_to_cLorentzVector(vec1) + c2*F2*convert_to_cLorentzVector(vec2) + c3*F3*convert_to_cLorentzVector(vec3);
  if ( verbosity_ >= 2 )
  {
    print("J2@2", retVal);
  }

  return retVal;
}

PolarimetricVectorTau2a1::cdouble
PolarimetricVectorTau2a1::comp_Fi(unsigned i,
                                  const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                                  DecayChannel decayChannel) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::Fi(i = " << i << ")>:\n";
  }

  double s1 = (p2 + p3).mag2();
  double s2 = (p1 + p3).mag2();
  double s3 = (p1 + p2).mag2();

  LorentzVector a = p1 + p2 + p3;
  double s = a.mag2();

  double m1, m2, m3;
  if ( decayChannel == k3ChargedPi )
  {
    m1 = m_chargedPi_;
    m2 = m_chargedPi_;
    m3 = m_chargedPi_;
  }
  else if ( decayChannel == kChargedPi2NeutralPi )
  {
    m1 = m_neutralPi_;
    m2 = m_neutralPi_;
    m3 = m_chargedPi_;
  }
  else
  {
    std::cerr << "Error in <PolarimetricVectorTau2a1::J2>: Invalid parameter 'decayChannel' = " << decayChannel << " !!\n";
    assert(0);
  }

  cdouble retVal;
  if ( i == 1 || i == 2 )
  {
    double sa, sb, ma, mb;
    if ( i == 1 )
    {
      sa = s1;
      sb = s2;
      ma = m1;
      mb = m2;
    }
    else if ( i == 2 )
    {
      sa = s2;
      sb = s1;
      ma = m2;
      mb = m1;
    }

    double s3Msa = (s3 - pow(m3, 2)) - (sa - pow(ma, 2));

    //retVal =  beta_[0]*BreitWigner(m0_rho770_, Gamma0_rho770_, sa, mb, m3, 1)
    //        + beta_[1]*BreitWigner(m0_rho1450_, Gamma0_rho1450_, sa, mb, m3, 1)
    //        - (1./3)*beta_[2]*s3Msa*BreitWigner(m0_rho770_, Gamma0_rho770_, sb, ma, m3, 1)
    //        - (1./3)*beta_[3]*s3Msa*BreitWigner(m0_rho1450_, Gamma0_rho1450_, sb, ma, m3, 1)
    //        + (1./3)*beta_[4]*((s - pow(m3, 2) + s3)*(2.*pow(ma, 2) + 2.*pow(mb, 2) - s3)/(6.*s3))*BreitWigner(m0_f2_, Gamma0_f2_, s3, ma, mb, 2)
    //        + (2./3)*beta_[5]*BreitWigner(m0_sigma_, Gamma0_sigma_, s3, ma, mb, 0)
    //        + (2./3)*beta_[6]*BreitWigner(m0_f0_, Gamma0_f0_, s3, m1, m2, 0, verbosity_);
    cdouble term1 =  beta_[0]*BreitWigner(m0_rho770_, Gamma0_rho770_, sa, mb, m3, 1, verbosity_);
    cdouble term2 =  beta_[1]*BreitWigner(m0_rho1450_, Gamma0_rho1450_, sa, mb, m3, 1, verbosity_);
    cdouble term3 = -(1./3)*beta_[2]*s3Msa*BreitWigner(m0_rho770_, Gamma0_rho770_, sb, ma, m3, 1, verbosity_);
    cdouble term4 = -(1./3)*beta_[3]*s3Msa*BreitWigner(m0_rho1450_, Gamma0_rho1450_, sb, ma, m3, 1, verbosity_);
    cdouble term5 =  (1./3)*beta_[4]*((s - pow(m3, 2) + s3)*(2.*pow(ma, 2) + 2.*pow(mb, 2) - s3)/(6.*s3))*BreitWigner(m0_f2_, Gamma0_f2_, s3, ma, mb, 2, verbosity_);
    cdouble term6 =  (2./3)*beta_[5]*BreitWigner(m0_sigma_, Gamma0_sigma_, s3, ma, mb, 0, verbosity_);
    cdouble term7 =  (2./3)*beta_[6]*BreitWigner(m0_f0_, Gamma0_f0_, s3, m1, m2, 0, verbosity_);
    if ( verbosity_ >= 2 )
    {
      std::cout << "F134 = " << -(1./3)*s3Msa << "\n";
      print("term1", term1);
      print("term2", term2);
      print("term3", term3);
      print("term4", term4);
      print("term5", term5);
      print("term6", term6);
      print("term7", term7);
    }
    retVal = term1 + term2 + term3 + term4 + term5 + term6 + term7;
  } 
  else if ( i == 3 )
  {
    double s1Ms2 = (s1 - pow(m1, 2)) - (s2 - pow(m2, 2));
    double s2Ms3 = (s2 - pow(m2, 2)) - (s3 - pow(m3, 2));
    double s3Ms1 = (s3 - pow(m3, 2)) - (s1 - pow(m1, 2));

    //retVal =  (1./3)*beta_[2]*(s2Ms3*BreitWigner(m0_rho770_, Gamma0_rho770_, s1, m2, m3, 1) + s3Ms1*BreitWigner(m0_rho770_, Gamma0_rho770_, s2, m1, m3, 1))
    //        + (1./3)*beta_[3]*(s2Ms3*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s1, m2, m3, 1) + s3Ms1*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s2, m1, m3, 1))
    //        - (1./2)*beta_[4]*s1Ms2*BreitWigner(m0_f2_, Gamma0_f2_, s3, m1, m2, 2);
    cdouble term3 =  (1./3)*beta_[2]*(s2Ms3*BreitWigner(m0_rho770_, Gamma0_rho770_, s1, m2, m3, 1, verbosity_) + s3Ms1*BreitWigner(m0_rho770_, Gamma0_rho770_, s2, m1, m3, 1, verbosity_));
    cdouble term4 =  (1./3)*beta_[3]*(s2Ms3*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s1, m2, m3, 1, verbosity_) + s3Ms1*BreitWigner(m0_rho1450_, Gamma0_rho1450_, s2, m1, m3, 1, verbosity_));
    cdouble term5 = -(1./2)*beta_[4]*s1Ms2*BreitWigner(m0_f2_, Gamma0_f2_, s3, m1, m2, 2, verbosity_);
    if ( verbosity_ >= 2 )
    {
      print("term3", term3);
      print("term4", term4);
      print("term5", term5);
    }
    retVal = term3 + term4 + term5;
  }
  else 
  {
    std::cerr << "Error in <PolarimetricVectorTau2a1::J2>: Invalid parameter 'i' = " << i << " !!\n";
    assert(0);
  }
  if ( verbosity_ >= 2 )
  {
    print("F3PIFactor", retVal);
  }
  retVal *= BreitWigner_a1(s);
  if ( verbosity_ >= 2 )
  {
    print(Form("F%u", i), retVal);
  }

  return retVal;
}

double
PolarimetricVectorTau2a1::m_a1(double s) const
{
  return m0_a1_;
}

double
PolarimetricVectorTau2a1::Gamma_a1(double s) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::Gamma_a1>:\n";
  } 

  double retVal = 0.;
  bool retVal_isInitialized = false;
  if ( s <= Gamma_a1_vs_s_.front().first )
  {    
    retVal = Gamma_a1_vs_s_.front().second;
    retVal_isInitialized = true;
  }
  else if ( s >= Gamma_a1_vs_s_.back().first )
  {
    retVal = Gamma_a1_vs_s_.back().second;
    retVal_isInitialized = true;
  }
  else
  {
    for ( size_t idx = 0; idx < (Gamma_a1_vs_s_.size() - 1); ++idx )
    {
      double s_lo = Gamma_a1_vs_s_[idx].first;
      double s_hi = Gamma_a1_vs_s_[idx + 1].first;
      if ( s >= s_lo && s <= s_hi )
      {
        double Gamma_lo = Gamma_a1_vs_s_[idx].second;
        double Gamma_hi = Gamma_a1_vs_s_[idx + 1].second;
        assert(Gamma_hi >= Gamma_lo);
        // CV: formula for linear interpolation taken from https://en.wikipedia.org/wiki/Linear_interpolation
        retVal = (Gamma_lo*(s_hi - s) + Gamma_hi*(s - s_lo))/(s_hi - s_lo);
        retVal_isInitialized = true;
        break;
      }
    }
  }
  if ( verbosity_ >= 2 )
  {
    std::cout << "Gamma(s = " << s << ") = " << retVal << "\n";
  }

  assert(retVal_isInitialized);
  return retVal;
}

PolarimetricVectorTau2a1::cdouble
PolarimetricVectorTau2a1::BreitWigner_a1(double s) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<PolarimetricVectorTau2a1::BreitWigner_a1>:\n";
  }

  double m = m_a1(s);
  double Gamma = Gamma_a1(s);
  if ( verbosity_ >= 2 )
  {
    std::cout << "m = " << m << "\n";
    std::cout << "Gamma = " << Gamma << "\n";
  }

  double num = -pow(m, 2);
  PolarimetricVectorTau2a1::cdouble denom(s - pow(m, 2), m0_a1_*Gamma);
  if ( verbosity_ >= 2 )
  {
    std::cout << "num = " << num << "\n";
    std::cout << "denom = " << denom << "\n";
  }

  PolarimetricVectorTau2a1::cdouble retVal = PolarimetricVectorTau2a1::cdouble(num, 0.)/denom;
  if ( verbosity_ >= 2 )
  {
    std::cout << "B_a1(s = " << s << ") = " << retVal << "\n";
  }

  return retVal;
}
