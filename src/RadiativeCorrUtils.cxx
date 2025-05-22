/**
 * This file contains utils which aim to study the effect of radiative corrections
 * qualitively before included in the GENIE event generator
 * It computes the correction factors to the cross section due to vertex, vacum or radiative effects
 * The methods used depend on the exclusive final state measured
 * References are provided for each case
 * \author Julia Tena Vidal \at Tel Aviv University                                                          
 * \date Nov 2023                                                                                                                   
 **/
#include <iostream>
#include <TF1.h>
#include "RadiativeCorrUtils.h"
#include "NuHepMC/HepMC3Features.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/WriterUtils.hxx"
#include "NuHepMC/make_writer.hxx"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"

using namespace std;
using namespace e4nu;

unsigned int utils::GetTargetNProtons( const unsigned int target_pdg ) {

  static const unsigned int kPdgH   = 1000010010; 
  static const unsigned int kPdgD   = 1000010020; 
  static const unsigned int kPdgHe3 = 1000020030; 
  static const unsigned int kPdgHe4 = 1000020040; 
  static const unsigned int kPdgC12 = 1000060120 ; 
  static const unsigned int kPdgFe56 = 1000260560 ;

  if ( target_pdg == kPdgH || target_pdg == 2212 || target_pdg == kPdgD ) return 1 ;
  else if ( target_pdg == kPdgHe3 || target_pdg == kPdgHe4 ) return 2 ; 
  else if ( target_pdg == kPdgC12 ) return 6 ; 
  else if ( target_pdg == kPdgFe56) return 26 ; 
  return 0;
}

HepMC3::FourVector utils::GetEmittedHardPhoton( const HepMC3::FourVector electron, double eloss ) {

  double ptLoss;
  double pzLoss;
  // for the z direction going probe theta = -nan 
  if (electron.pz()==electron.e()) {
    ptLoss = 0.;
    pzLoss = eloss;
  } else {
    ptLoss = eloss*sin(electron.theta()); 
    pzLoss = eloss*cos(electron.theta());
  }
  
  HepMC3::FourVector p4RadGamma(ptLoss*cos(electron.phi()),ptLoss*sin(electron.phi()),pzLoss,eloss);
  return p4RadGamma ;
}

double utils::SIMCBFactor( const double tgt_pdg ) { 
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double Z = utils::GetTargetNProtons(tgt_pdg);
  double L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
  double L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
  if( Z ==1 ) { 
    L1 = 5.31;
    L2 = 6.144;
  }
  double b = (1./9)*(12 + (Z+1)/(Z*L1 + L2));
  return b ;
}

double utils::GetLambda( const HepMC3::FourVector particle ) { 
  double lambda = TMath::Log(4*pow(particle.p3mod(),2)/pow(kElectronMass,2)) - 1 ;
  // The factor two is applied to account for the total radiative strenght for lambda 2, according to https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
  // Eq 61
  if( particle.pz() != particle.e() ) lambda += 2 * TMath::Log(0.5*(1-cos(particle.theta()))) ;
  lambda *= (kAem/kPi) ;
  if( lambda < 0 ) return 0;
  return lambda ; 
}


double utils::GetGFactor( const HepMC3::FourVector particle, const double tgt_pdg, const double thickness ) { 
  double b = SIMCBFactor( tgt_pdg );
  double lambda = GetLambda(particle);
  return lambda ; 
}

double utils::SIMCEnergyLoss(const HepMC3::FourVector particle, const double tgt_pdg, const double thickness, const double max_Ephoton, const double Delta_Em ) {
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double gfactor = GetGFactor( particle, tgt_pdg, thickness ); 
  double e_gamma_max = max_Ephoton;
  double e_gamma_min = 1E-25 ; // We generate the full distribution then cut it off
  double power_hi = pow(e_gamma_max,gfactor);
  double power_lo  = pow(Delta_Em,gfactor);
  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,gfactor);
  f->SetParameter(1,power_hi - power_lo);
  double energyLoss = f->GetRandom();
  if( energyLoss < 0 || energyLoss < Delta_Em ) energyLoss = 0 ; // cut off and set it to 0 -> unradiated event

  delete f;
  return energyLoss ; 
}

double utils::GetUFactor( const double mass, const double Q2 ){
  // computing u factor used for vacuum polarization calculation
  if ( Q2 == 0 ) return 0 ;
  return 4 * pow( mass, 2 ) / Q2 ;
}

double utils::VacuumPolarization( const double mass, const double Q2 ){
  // From Phys. Rev. C, 64:054610, Eq 28
  // In the limit Q2 >> m2, the vacumm polarization becomes: 
  //double delta_vac = 1./3. * ( - 5./3. + lnQm )  ;
  
  double u = GetUFactor( mass, Q2 ) ;
  double su = sqrt( 1 + u ) ;
  double pigg = (( 1 - u * 0.5 ) * su * TMath::Log( (su + 1)/(su - 1)) + u - 5./3.)/3.;
  pigg *= kAem ;
  return pigg; 
}

double utils::VanderhagenELoss( const double Q2 , const double Ee, const double resolution ) {
  // http://dx.doi.org/10.1103/PhysRevC.62.025501
  // Depends on event kinematics - not used
  // For QEL interactions, we can estiamte Q2. Similar results are obtaiend when using SIMC in this case
  double e_gamma_min = resolution ;
  double e_gamma_max = 0.2*Ee ;
  TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",e_gamma_min,e_gamma_max);
  double a = (kAem/kPi)*(TMath::Log(Q2)/pow(kElectronMass,2) - 1.);
  f->SetParameter(0,a);
  f->SetParameter(1,Ee);
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}

double utils::RadCorrWeight( const HepMC3::GenEvent & evt, const double true_Q2, const double thickness, const double max_Ephoton, const double Delta_E, const double Integral_peak, const double Integral_tail, const std::string model ){
  double weight = 1;
  // Grab beam and outgoing electron (after radiation - detected)
  auto beampt = NuHepMC::Event::GetBeamParticle(evt);
  auto primary_leptons = NuHepMC::Event::GetParticles_All(evt, NuHepMC::ParticleStatus::UndecayedPhysical, {beampt->pid()} );
  auto fslep = primary_leptons.back();
  // Grab intermediate leptons (vertex leptons)
  int id_radcorr = 123;
  auto corr_leptons = NuHepMC::Event::GetParticles_All(evt, id_radcorr, {beampt->pid()} );	

  // Grab event kinematic variables
  auto q = beampt->momentum() - fslep->momentum();
  double Q2 = true_Q2;
  unsigned int tgt_pdg = NuHepMC::Event::GetTargetParticle(evt)->pid();
  double Emax = max_Ephoton ; 
  double Emin = 1E-15;
  double radcorr_cutoff = 0.0001; 
  double delta = 0 ;
  
  if ( model == "vanderhaghen" ) { 
    // 10.1103/physrevc.62.025501
    double e_gamma_min = 1E-25;
    double lnQm  = TMath::Log(true_Q2/pow(kElectronMass,2)) ;
    double delta_vac = 2./3. * ( 5./3. - lnQm )  ;
    double delta_vertex = - 3./2. * lnQm + 2 + 0.5 * pow( lnQm, 2. ) + pow(kPi,2)/6. ; 
    double delta_soft = 0;//TMath::Log( beampt->momentum().e()*fslep->momentum().e()/ pow(Delta_E,2)) * ( lnQm - 1 ) ;
    //delta_soft += 0.5 * pow( TMath::Log( beampt->momentum().e()/fslep->momentum().e() ),2 ) ;
    //delta_soft -= 0.5 * pow( lnQm,2) - pow(kPi,2)/3. + TMath::DiLog(pow(TMath::Cos(fslep->momentum().theta())/2,2.));
    delta = (kAem/kPi) * ( delta_vac + delta_vertex + delta_soft ) ; 
  } else if ( model == "motsai" ) {
    // 10.1103/RevModPhys.41.205
    // https://inspirehep.net/files/1fcaa81f63f50d7bf56a22ce2c6b8b58 Equation II.2
    double SP = TMath::DiLog(-pow(TMath::Sin(fslep->momentum().theta())/2,2.));
    double Fth = TMath::Log(pow(TMath::Sin(fslep->momentum().theta()/2),2.))*TMath::Log(pow(TMath::Cos(fslep->momentum().theta()/2),2.))-SP;
    delta = ( TMath::Log( beampt->momentum().e() / Delta_E ) - (13./12.)) * (TMath::Log(Q2/pow(kElectronMass,2)) - 1) + 17./36. + 0.5*Fth;
    delta *= 2*(kAem/kPi);
  } 

  weight = 1 - delta ; 

  if ( model == "simc" ){
    // 10.1103/PhysRevC.64.054610
    // Here we have to be careful. Two weights are needed, depending on the event having hard bremstrahalung radiation or only virtual corrections
    // For the virtual correction, the weighting factor is (1-delta_hard)e^{-delta_soft(Delta_E)}
    // For the hard part, we only apply (1-delta_hard) * phi_ext_i * phi_ext_f * w_i * w_f 

    bool tail = false ;
    bool peak = false ;
    // Energy of the emmitted photons
    auto gamma_i = beampt->momentum() - corr_leptons[0]->momentum() ; // should be energy of decayed electron
    auto gamma_f =  fslep->momentum() - corr_leptons[1]->momentum() ; // should be energy of decayed electron

    if( gamma_i.e() < Delta_E && gamma_f.e() < Delta_E ) peak = true ; 
    else tail = true ;
    
    // Cross section correction when Egamma < Delta_E
    // Compute delta_hard:
    double lnQm  = TMath::Log(true_Q2/pow(kElectronMass,2)) ;

    // The contribution to the vacumm polarization from electron pairs is dominant. 
    // However, muon and hadronic pairs can contribute to a similar order for high Q2. 
    // We include the electron and muon polarization using the prenscription from Phys. Rev. C, 64:054610
    // The hadronic polarization is given by the Jegerlehner parameterization
    // use fit from the dispersion analysis of Jegerlehner (on-shell scheme).
    // Expression for polarization function for QED taken from Eq. (3.18) of Maximon & Tjon, PRC 62, 054320 (2000).
    double dvp_e   = VacuumPolarization( kElectronMass, true_Q2 ) ;
    double dvp_mu  = VacuumPolarization( kMuonMass, true_Q2 ) ;
    double dvp_tau = VacuumPolarization( kTauMass, true_Q2 ) ;
    double dvp_m1  = VacuumPolarization( km1, true_Q2 ) ;
    double dvp_m2  = VacuumPolarization( km2, true_Q2 ) ;
    double dvp_hadron = ka1 * dvp_m1 + ka2 * dvp_m2 ;
    double delta_vac = dvp_e + dvp_mu + dvp_tau + dvp_hadron ; 
    double delta_hard = 2 * kAem * ( -3./4. * lnQm + 1. - delta_vac ) /kPi ; 

    // Compute delta_soft(Delta_E) 
    double delta_soft = kAem / kPi * TMath::Log( beampt->momentum().e()*fslep->momentum().e()/ pow(Delta_E,2))*( lnQm - 1 ); 

    if( peak ) { 
      delta = delta_hard + delta_soft ;
      weight = ( 1 - delta_hard ) * exp( -delta_soft ) ;
      weight /= Integral_peak ; // renormalze peak to integral of flux. supression accounted for in weights 
      return weight ; 
    }

    // External weights
    double b = SIMCBFactor( tgt_pdg );
    double g_i = GetLambda( beampt->momentum() );
    double g_f = GetLambda( fslep->momentum() ) ;
    double phi_i = 1 ;

    if( beampt->momentum().p3mod() != 0 ) phi_i -= b * thickness * gamma_i.e() / beampt->momentum().p3mod() / g_i ; 
    double phi_f = 1 ;
    if( fslep->momentum().p3mod() != 0 ) phi_f -= b * thickness * gamma_f.e() / fslep->momentum().p3mod() / g_f ;

    // Soft weight 
    double gamma = tgamma( 1 + b * thickness ) ;
    double C_i = g_i / pow( beampt->momentum().p3mod(), b * thickness) ;
    C_i /= gamma ;
    C_i /= pow( sqrt( beampt->momentum().p3mod() * fslep->momentum().p3mod() ), GetLambda( beampt->momentum() ) ) ; 

    double C_f = g_i / pow( beampt->momentum().p3mod(), b * thickness) ;
    C_f /= gamma ;
    C_f /= pow( sqrt( beampt->momentum().p3mod() * fslep->momentum().p3mod() ), GetLambda( corr_leptons[1]->momentum() ) ) ; 

    double wght_simple_i = 1 ; 
    if( gamma_i.e() != 0 ) wght_simple_i = C_i * ( pow( max_Ephoton, g_i ) - pow( Delta_E, g_i ) ) / g_i ; 
    double wght_simple_f = 1 ; 

    weight = ( 1 - delta_hard ) * wght_simple_i * wght_simple_f * phi_i * phi_f ;

    // Here the theoretical calculation assumes that the tail has been generated using a p.d.f of only the tail which integral is 1
    // To account for the fact that the tail integral is not one, we need to add the wght_simple_i components and also the Phi_i for
    // the external radiation
    // However, note that the way we generated the events assumes that the whole p.d.f, which includes the peak at the beam energy, 
    // is equal to one. Hence, we need to reweigth also by the integral of our tail, so we re-normalize it to be 1
    weight /= Integral_tail ; 
    return weight ; 
  }

  return weight ; 
}
