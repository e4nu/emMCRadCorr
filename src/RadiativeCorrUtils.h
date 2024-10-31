#ifndef _RADUTILS_H_
#define _RADUTILS_H_

#include <string>
#include <TMath.h>
#include "NuHepMC/EventUtils.hxx"

namespace e4nu {
  namespace utils
  {
    // Constants
    static const double kAem = 1./137.03599976; // EM coupling const, dimensionless 
    static const double kAem2  = TMath::Power(kAem,2);
    static const double kPi = TMath::Pi(); 
    static const double kElectronMass = 0.000510998 ;
    static const double kMuonMass = 0.105658367 ;
    static const double kTauMass = 1.77682 ;
    static const double ka1 = 2.7299;
    static const double ka2 = 1.1115;
    static const double km1 = 0.209026;
    static const double km2 = 1.70554;

    // General functions
    unsigned int GetTargetNProtons( const unsigned int target_pdg ) ;
    HepMC3::FourVector GetEmittedHardPhoton( const HepMC3::FourVector electron, double eloss ) ; 
    double SIMCBFactor( const double tgt_pdg ) ;
    double GetLambda( const HepMC3::FourVector particle ) ;
    double GetGFactor( const HepMC3::FourVector particle, const double tgt_pdg, const double thickness ) ;
    double GetUFactor( const double mass, const double Q2 );

    // Energy loss probability functions
    double SIMCEnergyLoss(const HepMC3::FourVector particle, const double tgt_pdg, const double thickness, const double max_Ephoton, const double Delta_Em ) ;
    double VanderhagenELoss( const double Q2 , const double Ee, const double resolution ) ;

    // Weight calculators for cros section
    double RadCorrWeight( const HepMC3::GenEvent & evt, const double true_Q2, const double thickness, const double max_Ephoton, const double Delta_E, const double Integral_peak, const double Integral_tail, const std::string model );
    double VacuumPolarization( const double mass, const double Q2 );
  }
}

#endif
