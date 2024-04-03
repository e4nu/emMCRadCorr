#ifndef _RADUTILS_H_
#define _RADUTILS_H_

#include "TLorentzVector.h"
#include <string>

namespace e4nu {
  namespace utils
  {
    // Constants
    static const double kAem = 1./137.03599976; // EM coupling const, dimensionless 
    static const double kAem2  = TMath::Power(kAem,2);
    static const double kPi = TMath::Pi(); 
    static const double kElectronMass = 0.000510998 ;
    
    // QEL
    unsigned int GetTargetNProtons( const unsigned int target_pdg ) ;
    double RadCorrQELVertex( const double Q2 ) ; 
    double RadCorrQELVacumm( const double Q2 ) ; 
    double RadCorrQELRealRad( const double Q2, const double E, const double Ep, const double theta) ; 

    // Corrects outgoing electron for external radiation and also returns photon 4momenta
    TLorentzVector RadOutElectron( const TLorentzVector electron_vertex, TLorentzVector & out_electron, 
				   const int tgt, const double thickness, const double max_Ephoton, const std::string model ) ;
   
    // Computes total correction weight
    //    double SIMCRadCorrWeight( const e4nu::Event & event, const double thickness, const double max_Ephoton, const std::string model );

    // Probability funcitons
    double SIMCBFactor( const double tgt_pdg ) ;
    double SIMCEnergyLoss( const TLorentzVector particle, const int p_pdg, const double tgt_pdg, 
			   const double thickness, const double max_Ephoton ) ;
    double VanderhagenELoss( const double Q2 , const double Ee ) ;

    // Compute TLorentzVector for emited photon
    TLorentzVector GetEmittedHardPhoton( TLorentzVector electron, double eloss ) ; 
  }
}

#endif
