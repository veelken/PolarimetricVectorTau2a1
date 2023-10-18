#include "TauAnalysis/PolarimetricVectorTau2a1/interface/PolarimetricVectorTau2a1.h"

#include "TMath"    // TMath::Pi()

#include <assert.h> // assert()
#include <cmath>    // pow(), std::cos(), std::sin()
#include <iostream> // std::cerr, std::cout

const size_t numResonances = 7;

namespace
{
  /**
   * @brief Construct complex number, given its modulus and phase
   * @param modulus
   * @param phase
   * @return complex number
   */
  cdouble
  get_complex(double modulus, double phase)
  {
    double re = modulus*std::cos(phase);
    double im = modulus*std::sin(phase);
    return cdouble(re, im);
  }
}

PolarimetricVectorTau2a1::PolarimetricVectorTau2a1(int verbosity = -1)
  : verbosity_(verbosity)
{
  // define mass of charged and neutral pion;
  // values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01 (PDG)
  m_chargedPi_    = 0.139570; // [GeV]
  m_neutralPi_    = 0.134977; // [GeV]

  // define mass and width of a1(1260) meson;
  // values are taken from column "nominal fit" in Table VI of Phys.Rev.D 61 (2000) 012002
  m0_a1_;         = 1.331;    // [GeV]
  Gamma0_a1_;     = 0.814;    // [GeV]

  // define parameters specifying "running" of a1 width;
  // the values of Gamma_a1 as function of s have been taken from Fig. 9 (b) of Phys.Rev.D 61 (2000) 012002
  Gamma_a1_vs_s_.push_back(pdouble(0.00, 0.000);
  Gamma_a1_vs_s_.push_back(pdouble(0.20, 0.000);
  Gamma_a1_vs_s_.push_back(pdouble(0.40, 0.000);
  Gamma_a1_vs_s_.push_back(pdouble(0.50, 0.005);
  Gamma_a1_vs_s_.push_back(pdouble(0.60, 0.020);
  Gamma_a1_vs_s_.push_back(pdouble(0.65, 0.040);
  Gamma_a1_vs_s_.push_back(pdouble(0.70, 0.055);
  Gamma_a1_vs_s_.push_back(pdouble(0.75, 0.075);
  Gamma_a1_vs_s_.push_back(pdouble(0.80, 0.110);
  Gamma_a1_vs_s_.push_back(pdouble(0.85, 0.160);
  Gamma_a1_vs_s_.push_back(pdouble(0.90, 0.205);
  Gamma_a1_vs_s_.push_back(pdouble(0.95, 0.250);
  Gamma_a1_vs_s_.push_back(pdouble(1.00, 0.295);
  Gamma_a1_vs_s_.push_back(pdouble(1.05, 0.340);
  Gamma_a1_vs_s_.push_back(pdouble(1.10, 0.375);
  Gamma_a1_vs_s_.push_back(pdouble(1.15, 0.410);
  Gamma_a1_vs_s_.push_back(pdouble(1.20, 0.450);
  Gamma_a1_vs_s_.push_back(pdouble(1.25, 0.475);
  Gamma_a1_vs_s_.push_back(pdouble(1.30, 0.515);
  Gamma_a1_vs_s_.push_back(pdouble(1.35, 0.555);
  Gamma_a1_vs_s_.push_back(pdouble(1.40, 0.595);
  Gamma_a1_vs_s_.push_back(pdouble(1.45, 0.630);
  Gamma_a1_vs_s_.push_back(pdouble(1.50, 0.660);
  Gamma_a1_vs_s_.push_back(pdouble(1.55, 0.690);
  Gamma_a1_vs_s_.push_back(pdouble(1.60, 0.720);
  Gamma_a1_vs_s_.push_back(pdouble(1.65, 0.750);
  Gamma_a1_vs_s_.push_back(pdouble(1.70, 0.780);
  Gamma_a1_vs_s_.push_back(pdouble(1.75, 0.815);
  Gamma_a1_vs_s_.push_back(pdouble(1.80, 0.845);
  Gamma_a1_vs_s_.push_back(pdouble(1.85, 0.875);
  Gamma_a1_vs_s_.push_back(pdouble(1.90, 0.905);
  Gamma_a1_vs_s_.push_back(pdouble(1.93, 0.930);
  Gamma_a1_vs_s_.push_back(pdouble(1.95, 0.980);
  Gamma_a1_vs_s_.push_back(pdouble(2.00, 1.060);
  Gamma_a1_vs_s_.push_back(pdouble(2.05, 1.125);
  Gamma_a1_vs_s_.push_back(pdouble(2.10, 1.185);
  Gamma_a1_vs_s_.push_back(pdouble(2.15, 1.245);
  Gamma_a1_vs_s_.push_back(pdouble(2.20, 1.300);
  Gamma_a1_vs_s_.push_back(pdouble(2.25, 1.355);
  Gamma_a1_vs_s_.push_back(pdouble(2.30, 1.415);
  Gamma_a1_vs_s_.push_back(pdouble(2.35, 1.470);
  Gamma_a1_vs_s_.push_back(pdouble(2.37, 1.485);
  Gamma_a1_vs_s_.push_back(pdouble(2.40, 1.520);
  Gamma_a1_vs_s_.push_back(pdouble(2.45, 1.575);
  Gamma_a1_vs_s_.push_back(pdouble(2.50, 1.640);
  Gamma_a1_vs_s_.push_back(pdouble(2.55, 1.705);
  Gamma_a1_vs_s_.push_back(pdouble(2.60, 1.765);
  Gamma_a1_vs_s_.push_back(pdouble(2.65, 1.835);
  Gamma_a1_vs_s_.push_back(pdouble(2.70, 1.900);
  Gamma_a1_vs_s_.push_back(pdouble(2.75, 1.970);
  Gamma_a1_vs_s_.push_back(pdouble(2.80, 2.050);
  Gamma_a1_vs_s_.push_back(pdouble(2.85, 2.130);
  Gamma_a1_vs_s_.push_back(pdouble(2.90, 2.205);
  Gamma_a1_vs_s_.push_back(pdouble(2.95, 2.285);
  Gamma_a1_vs_s_.push_back(pdouble(3.00, 2.380);
  Gamma_a1_vs_s_.push_back(pdouble(3.05, 2.470);
  Gamma_a1_vs_s_.push_back(pdouble(3.10, 2.570);
  Gamma_a1_vs_s_.push_back(pdouble(3.15, 2.690);

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
    beta_.push_back(comp_beta(beta_moduli[idxResonance], beta_phases[idxResonance]));    
  }
}

PolarimetricVectorTau2a1::~PolarimetricVectorTau2a1()
{}

/**
 * @brief Compute polarimeter vector H according to Eq. (3.14) in Comput.Phys.Commun. 64 (1991) 275
 * @param p1 four-vector of first  same-sign pion [1,2]
 * @param p2 four-vector of second same-sign pion [1,2]
 * @param p3 four-vector of opposite-sign pion    [3]
 * @param charge tau lepton charge [4]
 * @param decayChannel flag indicating whether a1 decays into three charged pions or into one charged pion plus two neutral pions [5]
 *
 * The three pion four-vectors are expected to be given in the tau lepton rest-frame
 *
 *  [1] the two pi- in tau- -> pi- pi+ pi- nu decays
 *          two pi+ in tau+ -> pi+ pi- pi+ nu decays
 *          two neutral pions in tau+/- -> pi+/- pi0 pi0 nu decays
 *  [2] the ordering of the two same-sign pions does not matter
 *  [3] the pi+ in tau- -> pi- pi+ pi- nu decays
 *          pi- in tau+ -> pi+ pi- pi+ nu decays
 *          charged pion in tau+/- -> pi+/- pi0 pi0 nu decays
 *  [4] use -1 for tau-
 *          +1 for tau+
 *  [5] use k3ChargedPi for tau- -> pi- pi+ pi- nu and tau+ -> pi+ pi- pi+ nu decays
 *          kChargedPi2NeutralPi for tau+/- -> pi+/- pi0 pi0 nu decays
 *
 * @return H
 */
Vector
PolarimetricVectorTau2a1::operator()(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                                     int charge,
                                     DecayChannel decayChannel)
{

}

namespace
{
  /**
   * @brief Return element of metric tensor g at index {mu,nu}
   * @param mu first index
   * @param nu second index
   * @return g_{mu,nu}
   */
  double
  g(size_t mu, size_t nu)
  {
    if      ( (mu == 0 && nu == 0) || (mu == 1 && nu == 1) || (mu == 2 && nu == 2) ) return -1.;
    else if (  mu == 3 && nu == 3                                                  ) return +1.;
    else                                                                             return  0.;
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
  get_component(const LorentzVector& p, size_t mu)
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
  cLorentzVector
  convert_to_cLorentzVector(const LorentzVector& p)
  {
    cLorentzVector retVal;
    for ( size_t mu = 0; mu < 4; ++mu )
    {
      retVal(mu) = get_component(p, size_t mu);
    }
    return retVal;
  }
}

/**
 * @brief Compute hadronic current J according to Eq. (3) in Phys.Rev.D 61 (2000) 012002
 * @param p1 four-vector of first  same-sign pion [1,2]
 * @param p2 four-vector of second same-sign pion [1,2]
 * @param p3 four-vector of opposite-sign pion    [3] 
 * @param decayChannel flag indicating whether a1 decays into three charged pions or into one charged pion plus two neutral pions [4]
 *
 * The three pion four-vectors are expected to be given in the tau lepton rest-frame
 *
 *  [1] the two pi- in tau- -> pi- pi+ pi- nu decays
 *          two pi+ in tau+ -> pi+ pi- pi+ nu decays
 *          two neutral pions in tau+/- -> pi+/- pi0 pi0 nu decays
 *  [2] the ordering of the two same-sign pions does not matter
 *  [3] the pi+ in tau- -> pi- pi+ pi- nu decays
 *          pi- in tau+ -> pi+ pi- pi+ nu decays
 *          charged pion in tau+/- -> pi+/- pi0 pi0 nu decays
 *  [4] use k3ChargedPi for tau- -> pi- pi+ pi- nu and tau+ -> pi+ pi- pi+ nu decays
 *          kChargedPi2NeutralPi for tau+/- -> pi+/- pi0 pi0 nu decays
 *
 * @return J
 */
cLorentzVector
PolarimetricVectorTau2a1::J(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3,
                            DecayChannel decayChannel)
{
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

  cTensor T;
  for ( size_t mu = 0; mu < 4; ++mu )
  {
    for ( size_t nu = 0; nu < 4; ++nu )
    {
      T(mu,nu) = g(mu,nu) - get_component(a, mu)*get_component(a, nu)/s;
    }
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

  cLorentzVectorretVal retVal;
  for ( size_t idx = 0; idx < numResonances; ++idx )
  {
    retVal += beta_[idx]*j[idx];
  }
  retVal *= BreitWigner_a1(s);
  return retVal;
}

/**
 * @brief Compute "running mass" of a1 meson.
 *        The "running mass" of the a1 meson can only be computed for non-zero values of the a1 size parameter R.
 *        Since the nominal fit sets R equal to zero, we take the a1 mass to be constant,
 *        following Section VI.C.2 in Phys.Rev.D 61 (2000) 012002.
 * @param s (mass of three-pion system)^2
 * @return m_a1(s)
 */
double
PolarimetricVectorTau2a1::m_a1(double s) const
{
  return m0_a1_;
}

/**
 * @brief Compute "running width" of a1 meson.
 *        Instead of implementing the "running" of the a1 width according to Eq. (11),
 *        we interpolate (linearly) between the Gamma_tot values shown in Fig. 9 (b) of Phys.Rev.D 61 (2000) 012002
 * @param s (mass of three-pion system)^2
 * @return Gamma_a1(s)
 */
double
PolarimetricVectorTau2a1::Gamma_a1(double s) const
{
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
        assert(Gamma_hi > Gamma_lo);
        // CV: formula for linear interpolation taken from https://en.wikipedia.org/wiki/Linear_interpolation
        retVal = (Gamma_lo*(s_hi - s) + Gamma_hi*(s - s_lo))/(s_hi - s_lo);
        retVal_isInitialized = true;
        break;
      }
    }
  }
  assert(retVal_isInitialized);
  return retVal;
}

/**
 * @brief Compute Breit-Wigner function of a1 meson, given by first term on the right-hand-side of Eq. (7) in Phys.Rev.D 61 (2000) 012002
 * @param s (mass of three-pion system)^2
 * @return B_a1(s)
 */
cdouble
PolarimetricVectorTau2a1::BreitWigner_a1(double s) const
{
  double retVal = cdouble(s - pow(m_a1(s), 2), m0_a1_*Gamma_a1(s));
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
cdouble
PolarimetricVectorTau2a1::BreitWigner(double m0, double Gamma0, double si, double mj, double mk, unsigned L) const;
{
  cdouble retVal = cdouble(pow(m0, 2), 0.)/cdouble(pow(m0, 2) - si, -m0*Gamma(m0, Gamma0, si, mj, mk, L));
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
PolarimetricVectorTau2a1::Gamma(double m0, double Gamma0, double si, double mj, double mk, unsigned L) const
{
  // compute "decay momentum" and "decay momentum at s_i=m0^2" according to Eq. (A6) in Phys.Rev.D 61 (2000) 012002
  double kdashi = kdash(si, mj, mk);
  double kdash0 = kdash(pow(m0, 2), mj, mk);

  double retVal = Gamma0*pow(kdashi/kdash0, 2*L + 1)*m0/std::sqrt(si);
  return retVal;
}


