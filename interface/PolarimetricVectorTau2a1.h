#ifndef PolarimetricVectorTau2a1_h
#define PolarimetricVectorTau2a1_h

/**
 * PolarimetricVectorTau2a1
 *
 * Authors: Vladimir Cherepanov (vladimir.cherepanov@cern.ch)
 *          Christian Veelken   (christian.veelken@cern.ch)
 *
 * Compute polarimetric vector for tau->a1->pi- pi0 pi0 nu and tau->a1->pi- pi+ pi- nu decays.
 *
 * The computation of the polarimetric vector is based on Section 3.2 of the paper:
 *   Comput.Phys.Commun. 64 (1991) 275,
 * with the hadronic current J taken from the measurement of tau->a1->pi- pi0 pi0 nu decays 
 * by the CLEO collaboration, published in the paper:
 *   Phys.Rev.D 61 (2000) 012002
 *
 * If you use this code, please cite: 
 *   XXXXXX
 *
 */

#include "Math/Cartesian3D.h"   // ROOT::Math::Cartesian3D<>
#include "Math/LorentzVector.h" // ROOT::Math::LorentzVector<>
#include "Math/PtEtaPhiE4D.h"   // ROOT::Math::PxPyPzE4D<>
#include "Math/SMatrix.h"       // ROOT::Math::SMatrix<>
#include "Math/SVector.h"       // ROOT::Math::SVector<>
#include "Math/Vector3D.h"      // ROOT::Math::DisplacementVector3D<>

#include <complex>              // std::complex<>
#include <utility>              // std::pair<>
#include <vector>               // std::vector<>

class PolarimetricVectorTau2a1
{ 
 public:
  typedef std::complex<double> cdouble;
  typedef std::pair<double,double> pdouble;

  // in cTensor and cLorentzVector objects, the index 
  //   0 refers to px
  //   1 refers to py
  //   2 refers to pz
  //   3 refers to energy
  typedef ROOT::Math::SMatrix<cdouble,4,4> cTensor;
  typedef ROOT::Math::SVector<cdouble,4> cLorentzVector;

  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> Vector;

  PolarimetricVectorTau2a1(int verbosity = -1);
  ~PolarimetricVectorTau2a1();

  enum DecayChannel { k3ChargedPi, kChargedPi2NeutralPi };

  // 
  Vector
  operator()(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3, 
             int charge,
             DecayChannel decayChannel);

 private:
  // hadronic current J,
  // computed according to Eq. (3) in Phys.Rev.D 61 (2000) 012002
  cLorentzVector
  J(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3, 
    DecayChannel decayChannel) const;

  // "running mass" of the a1 meson
  double
  m_a1(double s) const;

  // "running width" of the a1 meson,
  // computed according to Eq. (11) in Phys.Rev.D 61 (2000) 012002
  double
  Gamma_a1(double s) const;

  // Breit-Wigner function of a1 meson,
  // computed according to Eq. (7) in Phys.Rev.D 61 (2000) 012002
  cdouble
  BreitWigner_a1(double s) const;

  // Breit-Wigner function of intermediate rho(770), rho(1450), f2(1270), sigma, f0(1370) resonances,
  // computed according to top line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
  cdouble
  BreitWigner(double m0, double Gamma0, double si, double mj, double mk, unsigned L) const;

  // "running width" of intermediate rho(770), rho(1450), f2(1270), sigma, f0(1370) resonances,
  // computed according to bottom line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
  double
  Gamma(double m0, double Gamma0, double si, double mj, double mk, unsigned L) const;

  // mass of charged and neutral pion
  double m_chargedPi_;    // mass of charged pion     [GeV]
  double m_neutralPi_;    // mass of neutral pion     [GeV]

  // mass and width of a1(1260) meson
  double m0_a1_;          // mass of a1(1260) meson   [GeV]
  double Gamma0_a1_;      // width of a1(1260) meson  [GeV]

  // parameters specifying "running" of a1 width
  std::vector<pdouble> Gamma_a1_vs_s_;

  // masses and widths of intermediate rho(770), rho(1450), f2(1270), sigma, f0(1370) resonances
  double m0_rho770_;      // mass of rho(770) meson   [GeV]
  double Gamma0_rho770_;  // width of rho(770) meson  [GeV]
  double m0_rho1450_;     // mass of rho(1450) meson  [GeV]
  double Gamma0_rho1450_; // width of rho(1450) meson [GeV]
  double m0_f2_;          // mass of f2(1270) meson   [GeV]
  double Gamma0_f2_;      // width of f2(1270) meson  [GeV]
  double m0_sigma_;       // mass of sigma meson      [GeV]
  double Gamma0_sigma_;   // width of sigma meson     [GeV]
  double m0_f0_;          // mass of f0(1370) meson   [GeV]
  double Gamma0_f0_;      // width of f0(1370) meson  [GeV]

  // coefficients specifying the contribution of meson resonances to the hadronic current J,
  // cf. Eq. (3) of Phys.Rev.D 61 (2000) 012002
  std::vector<cdouble> beta_;

  // verbosity:
  //   -1 no output
  //    1 minimal output
  //  > 1 extensive output (to be used only for debugging purposes)
  int verbosity_;
};

#endif
