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
  static const unsigned int kPdgHe3 = 1000020030; 
  static const unsigned int kPdgHe4 = 1000020040; 
  static const unsigned int kPdgC12 = 1000060120 ; 
  static const unsigned int kPdgFe56 = 1000260560 ;

  if ( target_pdg == kPdgH || target_pdg == 2212 ) return 1 ;
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

double utils::SIMCEnergyLoss(const HepMC3::FourVector particle, const double tgt_pdg, const double thickness, const double max_Ephoton, const double resolution ) {
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double b = SIMCBFactor( tgt_pdg );
  double lambda = TMath::Log(4*pow(particle.p3mod(),2)/pow(kElectronMass,2)) - 1 ;
  if( particle.pz() != particle.e() ) lambda += TMath::Log(0.5*(1-cos(particle.theta()))) ;
  lambda *= (kAem/kPi) ;
  lambda += b*thickness;
  if( lambda < 0 ) return 0; 

  double e_gamma_max = max_Ephoton;
  double e_gamma_min = 1E-25;;//resolution ;
  double power_hi = pow(e_gamma_max,lambda);
  double power_lo  = pow(e_gamma_min,lambda);
  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,lambda);
  f->SetParameter(1,power_hi - power_lo);
  double energyLoss = f->GetRandom();

  delete f;
  return energyLoss ; 
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

double utils::RadCorrWeight( const HepMC3::GenEvent & evt, const double true_Q2, const double thickness, const double max_Ephoton, const double Delta_E, const std::string model ){
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
  } else if ( model == "simc" ){
    // 10.1103/PhysRevC.64.054610
    delta = 2.*(kAem/kPi)*( -13./12.*TMath::Log(true_Q2/pow(kElectronMass,2))+14./9.);
    double lnQm  = TMath::Log(true_Q2/pow(kElectronMass,2)) ;
    double delta_vac = 1./3./kPi * ( - 5./3. + lnQm )  ;
    double delta_hard = 2 * kAem * ( -3./4./kPi * lnQm + 1./kPi - delta_vac ) ; 
    delta = delta_hard ;
  }
  weight = 1 - delta ; 
  
  return weight ; 
}
