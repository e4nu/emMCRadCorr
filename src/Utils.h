/**
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date March 2024
 **/

#ifndef _UTILS_H_
#define _UTILS_H_
#include <iostream>
#include <string> 
#include <TMath.h>
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
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

// Pdg codes 
static const int kPdgProton = 2212 ;
static const int kPdgNeutron = 2112 ;
static const int kPdgPiP = 211 ; 
static const int kPdgPiM = -211 ; 
static const int kPdgPi0 = 111 ; 
static const int kPdgKP = 321 ;
static const int kPdgKM = -321 ;
static const int kPdgK0 = 311 ;
static const int kPdgElectron = 11 ; 
static const int kPdgPositron = -11 ; 
static const int kPdgPhoton = 22 ; 

// Define branch variables as defined in GENIE
//
int    Iev         = 0;      // Event number 
int    Neutrino    = 0;      // Neutrino pdg code
int    FSPrimLept  = 0;      // Final state primary lepton pdg code
int    Target      = 0;      // Nuclear target pdg code (10LZZZAAAI)
int    TargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
int    TargetA     = 0;      // Nuclear target A (extracted from pdg code above)
int    HitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
int    HitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
bool   FromSea     = false;  // Hit quark is from sea     (set for DIS events only)
int    ResId       = 0;      // Produced baryon resonance (set for resonance events only)
bool   IsQel       = false;  // Is QEL?
bool   IsRes       = false;  // Is RES?
bool   IsDis       = false;  // Is DIS?
bool   IsCoh       = false;  // Is Coherent?
bool   IsMec       = false;  // Is MEC?
bool   IsDfr       = false;  // Is Diffractive?
bool   IsImd       = false;  // Is IMD?
bool   IsNrm       = false;  // Is Norm?
bool   IsSingleK   = false;  // Is single kaon?  
bool   IsImdAnh    = false;  // Is IMD annihilation?
bool   IsNuEL      = false;  // Is ve elastic?
bool   IsEM        = false;  // Is EM process?
bool   IsCC        = false;  // Is Weak CC process?
bool   IsNC        = false;  // Is Weak NC process?
bool   IsCharmPro  = false;  // Produces charm?
bool   IsAMNuGamma = false;  // is anomaly mediated nu gamma
bool   IsHNL       = false;  // is HNL decay?
int    CodeNeut    = 0;      // The equivalent NEUT reaction code (if any)
int    CodeNuance  = 0;      // The equivalent NUANCE reaction code (if any)
double Weight      = 0;      // Event weight
double KineXs      = 0;      // Bjorken x as was generated during kinematical selection; takes fermi momentum / off-shellness into account
double KineYs      = 0;      // Inelasticity y as was generated during kinematical selection; takes fermi momentum / off-shellness into account
double KineTs      = 0;      // Energy transfer to nucleus at COH events as was generated during kinematical selection
double KineQ2s     = 0;      // Momentum transfer Q^2 as was generated during kinematical selection; takes fermi momentum / off-shellness into account
double KineWs      = 0;      // Hadronic invariant mass W as was generated during kinematical selection; takes fermi momentum / off-shellness into account
double KineX       = 0;      // Experimental-like Bjorken x; neglects fermi momentum / off-shellness 
double KineY       = 0;      // Experimental-like inelasticity y; neglects fermi momentum / off-shellness 
double KineT       = 0;      // Experimental-like energy transfer to nucleus at COH events 
double KineQ2      = 0;      // Experimental-like momentum transfer Q^2; neglects fermi momentum / off-shellness
double KineW       = 0;      // Experimental-like hadronic invariant mass W; neglects fermi momentum / off-shellness 
double EvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
double Ev          = 0;      // Neutrino energy @ LAB
double Pxv         = 0;      // Neutrino px @ LAB
double Pyv         = 0;      // Neutrino py @ LAB
double Pzv         = 0;      // Neutrino pz @ LAB
double Evcorr      = 0;      // Neutrino energy @ LAB after rad corr (vertex)
double Pxvcorr     = 0;      // Neutrino px @ LAB after rad corr (vertex)
double Pyvcorr     = 0;      // Neutrino py @ LAB after rad corr (vertex)
double Pzvcorr     = 0;      // Neutrino pz @ LAB after rad corr (vertex)
double En          = 0;      // Initial state hit nucleon energy @ LAB
double Pxn         = 0;      // Initial state hit nucleon px @ LAB
double Pyn         = 0;      // Initial state hit nucleon py @ LAB
double Pzn         = 0;      // Initial state hit nucleon pz @ LAB
double El          = 0;      // Final state primary lepton energy @ LAB
double Pxl         = 0;      // Final state primary lepton px @ LAB
double Pyl         = 0;      // Final state primary lepton py @ LAB
double Pzl         = 0;      // Final state primary lepton pz @ LAB
double Pl          = 0;      // Final state primary lepton p  @ LAB
double Costhl      = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
double Elcorr      = 0;      // Final state primary lepton energy @ LAB before rad corr (vertex)
double Pxlcorr     = 0;      // Final state primary lepton px @ LAB before rad corr (vertex)
double Pylcorr     = 0;      // Final state primary lepton py @ LAB before rad corr (vertex)
double Pzlcorr     = 0;      // Final state primary lepton pz @ LAB before rad corr (vertex)
double Plcorr      = 0;      // Final state primary lepton p  @ LBA before rad corr (vertex)
double Costhlcorr  = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
int    NfP         = 0;      // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
int    NfN         = 0;      // Nu. of final state n's + \bar{n}'s
int    NfPip       = 0;      // Nu. of final state pi+'s
int    NfPim       = 0;      // Nu. of final state pi-'s
int    NfPi0       = 0;      // Nu. of final state pi0's (
int    NfKp        = 0;      // Nu. of final state K+'s
int    NfKm        = 0;      // Nu. of final state K-'s
int    NfK0        = 0;      // Nu. of final state K0's + \bar{K0}'s
int    NfEM        = 0;      // Nu. of final state gammas and e-/e+ 
int    NfOther     = 0;      // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
int    NiP         = 0;      // Nu. of `primary' (: before intranuclear rescattering) p's + \bar{p}'s  
int    NiN         = 0;      // Nu. of `primary' n's + \bar{n}'s  
int    NiPip       = 0;      // Nu. of `primary' pi+'s 
int    NiPim       = 0;      // Nu. of `primary' pi-'s 
int    NiPi0       = 0;      // Nu. of `primary' pi0's 
int    NiKp        = 0;      // Nu. of `primary' K+'s  
int    NiKm        = 0;      // Nu. of `primary' K-'s  
int    NiK0        = 0;      // Nu. of `primary' K0's + \bar{K0}'s 
int    NiEM        = 0;      // Nu. of `primary' gammas and e-/e+ 
int    NiOther     = 0;      // Nu. of other `primary' hadron shower particles
int    Nf          = 0;      // Nu. of final state particles in hadronic system
int    Pdgf  [500];       // Pdg code of k^th final state particle in hadronic system
double Ef    [500];       // Energy     of k^th final state particle in hadronic system @ LAB
double Pxf   [500];       // Px         of k^th final state particle in hadronic system @ LAB
double Pyf   [500];       // Py         of k^th final state particle in hadronic system @ LAB
double Pzf   [500];       // Pz         of k^th final state particle in hadronic system @ LAB
double Pf    [500];       // P          of k^th final state particle in hadronic system @ LAB
double Costhf[500];       // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
int    Ni = 0;               // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
int    Pdgi[500];         // Pdg code of k^th particle in 'primary' hadronic system 
int    Resc[500];         // FSI code of k^th particle in 'primary' hadronic system 
double Ei  [500];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
double Pxi [500];         // Px       of k^th particle in 'primary' hadronic system @ LAB
double Pyi [500];         // Py       of k^th particle in 'primary' hadronic system @ LAB
double Pzi [500];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
double VtxX;                 // Vertex x in detector coord system (SI)
double VtxY;                 // Vertex y in detector coord system (SI)
double VtxZ;                 // Vertex z in detector coord system (SI)
double VtxT;                 // Vertex t in detector coord system (SI)
double SumKEf;               // Sum of kinetic energies of all final state particles
double CalResp0;             // Approximate calorimetric response to the hadronic system computed as sum of
Double_t  XSec;              // the event cross section in 1E-38cm^2
Double_t  DXSec;             // is the differential cross section for the selected in 1E-38cm^2/{K^n}
UInt_t    KPS;               // phase space that the xsec has been evaluated into

namespace e4nu { 
  namespace utils
  {
    void SetGSTBranchAddress( TTree * output_tree ){
      // Create tree branches
      //
      output_tree->Branch("iev", &Iev, "iev/I"  );
      output_tree->Branch("neu", &Neutrino, "neu/I" );
      output_tree->Branch("fspl",  &FSPrimLept, "fspl/I");
      output_tree->Branch("tgt", &Target, "tgt/I" );
      output_tree->Branch("Z",   &TargetZ,  "Z/I" );
      output_tree->Branch("A",   &TargetA,  "A/I" );
      output_tree->Branch("hitnuc", &HitNuc, "hitnuc/I" );
      output_tree->Branch("hitqrk", &HitQrk, "hitqrk/I" );
      output_tree->Branch("resid",  &ResId, "resid/I" );
      output_tree->Branch("sea", &FromSea,  "sea/O" );
      output_tree->Branch("qel", &IsQel, "qel/O" );
      output_tree->Branch("mec", &IsMec, "mec/O" );
      output_tree->Branch("res", &IsRes, "res/O" );
      output_tree->Branch("dis", &IsDis, "dis/O" );
      output_tree->Branch("coh", &IsCoh,  "coh/O" );
      output_tree->Branch("dfr", &IsDfr,  "dfr/O" );
      output_tree->Branch("imd", &IsImd, "imd/O" );
      output_tree->Branch("norm", &IsNrm,  "norm/O" );
      output_tree->Branch("imdanh", &IsImdAnh, "imdanh/O" );
      output_tree->Branch("singlek",  &IsSingleK,  "singlek/O"  );  
      output_tree->Branch("nuel", &IsNuEL, "nuel/O" );
      output_tree->Branch("em", &IsEM, "em/O" );
      output_tree->Branch("cc", &IsCC, "cc/O" );
      output_tree->Branch("nc", &IsNC, "nc/O" );
      output_tree->Branch("charm",  &IsCharmPro, "charm/O" );
      output_tree->Branch("amnugamma",  &IsAMNuGamma,   "amnugamma/O"   );
      output_tree->Branch("hnl", &IsHNL,  "hnl/O"  );
      output_tree->Branch("neut_code",  &CodeNeut, "neut_code/I"   );
      output_tree->Branch("nuance_code",   &CodeNuance, "nuance_code/I" );
      output_tree->Branch("wght", &Weight, "wght/D" );
      output_tree->Branch("xs", &KineXs, "xs/D" );
      output_tree->Branch("ys", &KineYs, "ys/D" );
      output_tree->Branch("ts", &KineTs, "ts/D" );
      output_tree->Branch("Q2s", &KineQ2s,  "Q2s/D" );
      output_tree->Branch("Ws", &KineWs, "Ws/D" );
      output_tree->Branch("x", &KineX, "x/D" );
      output_tree->Branch("y", &KineY, "y/D" );
      output_tree->Branch("t", &KineT, "t/D" );
      output_tree->Branch("Q2", &KineQ2, "Q2/D" );
      output_tree->Branch("W", &KineW, "W/D" );
      output_tree->Branch("EvRF", &EvRF, "EvRF/D" );
      output_tree->Branch("Ev", &Ev, "Ev/D" );
      output_tree->Branch("pxv", &Pxv, "pxv/D" );
      output_tree->Branch("pyv", &Pyv, "pyv/D" );
      output_tree->Branch("pzv", &Pzv, "pzv/D" );
      output_tree->Branch("Evcorr", &Evcorr, "Evcorr/D" );
      output_tree->Branch("pxvcorr", &Pxvcorr, "pxvcorr/D" );
      output_tree->Branch("pyvcorr", &Pyvcorr, "pyvcorr/D" );
      output_tree->Branch("pzvcorr", &Pzvcorr, "pzvcorr/D" );
      output_tree->Branch("En", &En, "En/D" );
      output_tree->Branch("pxn", &Pxn, "pxn/D" );
      output_tree->Branch("pyn", &Pyn, "pyn/D" );
      output_tree->Branch("pzn", &Pzn, "pzn/D" );
      output_tree->Branch("El", &El, "El/D" );
      output_tree->Branch("pxl", &Pxl, "pxl/D" );
      output_tree->Branch("pyl", &Pyl, "pyl/D" );
      output_tree->Branch("pzl", &Pzl, "pzl/D" );
      output_tree->Branch("pl",  &Pl,  "pl/D" );
      output_tree->Branch("cthl", &Costhl, "cthl/D" );
      output_tree->Branch("Elcorr", &Elcorr, "Elcorr/D" );
      output_tree->Branch("pxlcorr", &Pxlcorr, "pxlcorr/D" );
      output_tree->Branch("pylcorr", &Pylcorr, "pylcorr/D" );
      output_tree->Branch("pzlcorr", &Pzlcorr, "pzlcorr/D" );
      output_tree->Branch("plcorr",  &Plcorr,  "plcorr/D" );
      output_tree->Branch("cthlcorr", &Costhlcorr, "cthlcorr/D" );
      output_tree->Branch("nfp", &NfP, "nfp/I" );
      output_tree->Branch("nfn", &NfN, "nfn/I" );
      output_tree->Branch("nfpip",  &NfPip, "nfpip/I" );
      output_tree->Branch("nfpim",  &NfPim, "nfpim/I" );
      output_tree->Branch("nfpi0",  &NfPi0, "nfpi0/I" );
      output_tree->Branch("nfkp", &NfKp, "nfkp/I" );
      output_tree->Branch("nfkm", &NfKm, "nfkm/I" );
      output_tree->Branch("nfk0", &NfK0, "nfk0/I" );
      output_tree->Branch("nfem", &NfEM, "nfem/I" );
      output_tree->Branch("nfother",  &NfOther,  "nfother/I"  );
      output_tree->Branch("nip", &NiP, "nip/I" );
      output_tree->Branch("nin", &NiN, "nin/I" );
      output_tree->Branch("nipip",  &NiPip, "nipip/I" );
      output_tree->Branch("nipim",  &NiPim, "nipim/I" );
      output_tree->Branch("nipi0",  &NiPi0, "nipi0/I" );
      output_tree->Branch("nikp", &NiKp, "nikp/I" );
      output_tree->Branch("nikm", &NiKm, "nikm/I" );
      output_tree->Branch("nik0", &NiK0, "nik0/I" );
      output_tree->Branch("niem", &NiEM, "niem/I" );
      output_tree->Branch("niother",  &NiOther,  "niother/I"  );
      output_tree->Branch("ni",  &Ni,  "ni/I" );
      output_tree->Branch("pdgi", Pdgi, "pdgi[ni]/I"   );
      output_tree->Branch("resc", Resc, "resc[ni]/I"   );
      output_tree->Branch("Ei", Ei,  "Ei[ni]/D" );
      output_tree->Branch("pxi", Pxi, "pxi[ni]/D"  );
      output_tree->Branch("pyi", Pyi, "pyi[ni]/D"  );
      output_tree->Branch("pzi", Pzi, "pzi[ni]/D"  );
      output_tree->Branch("nf",  &Nf,  "nf/I" );
      output_tree->Branch("pdgf", Pdgf, "pdgf[nf]/I"   );
      output_tree->Branch("Ef", Ef,  "Ef[nf]/D" );
      output_tree->Branch("pxf", Pxf, "pxf[nf]/D"  );
      output_tree->Branch("pyf", Pyf, "pyf[nf]/D"  );
      output_tree->Branch("pzf", Pzf, "pzf[nf]/D"  );
      output_tree->Branch("pf", Pf, "pf[nf]/D" );
      output_tree->Branch("cthf", Costhf,  "cthf[nf]/D" );
      output_tree->Branch("vtxx",  &VtxX, "vtxx/D" );
      output_tree->Branch("vtxy",  &VtxY, "vtxy/D" );
      output_tree->Branch("vtxz",  &VtxZ, "vtxz/D" );
      output_tree->Branch("vtxt",  &VtxT, "vtxt/D" );
      output_tree->Branch("sumKEf",  &SumKEf, "sumKEf/D" );
      output_tree->Branch("calresp0",  &CalResp0, "calresp0/D" );
      output_tree->Branch("XSec",  &XSec, "XSec/D" );
      output_tree->Branch("DXSec",  &DXSec, "DXSec/D" );
      output_tree->Branch("KPS", &KPS, "KPS/i" );
    }

    void StoreHepMCToGST( const HepMC3::GenEvent evt, TTree * output_tree ){
      Iev = evt.event_number();
      auto in_gen_run_info = evt.run_info();
      bool is_GENIE = false ;
      for(auto const &tool : in_gen_run_info->tools()){
	if( tool.name == "GENIE") is_GENIE = true ; 
      }
      //HepMC3::Print::content(evt);

      auto beampt = NuHepMC::Event::GetBeamParticle(evt);
      Ev = beampt->momentum().e();
      Pxv = beampt->momentum().px();
      Pyv = beampt->momentum().py();
      Pzv = beampt->momentum().pz();
      Neutrino = beampt->pid();

      auto primary_leptons = NuHepMC::Event::GetParticles_All(evt, NuHepMC::ParticleStatus::UndecayedPhysical, {beampt->pid()} );	
      auto fslep = primary_leptons.back();  
      FSPrimLept = fslep->pid();
      El = fslep->momentum().e();
      Pxl = fslep->momentum().px();
      Pyl = fslep->momentum().py();
      Pzl = fslep->momentum().pz();
      Pl = fslep->momentum().p3mod();
      Costhl = cos(fslep->momentum().theta());

      if( beampt->pid() == kPdgElectron ) IsEM=true ; 
      else if ( TMath::Abs(beampt->pid()) == 12 || TMath::Abs(beampt->pid()) == 14 || TMath::Abs(beampt->pid()) == 16 ) {
	if( FSPrimLept == beampt->pid() ) IsNC=true;
	else IsCC=true;
      }

      if(NuHepMC::GC1::SignalsConvention(evt.run_info(),"E.C.1")){
	//the generator has promised that 300 <= process_id < 349 == IsCCRes and 350 <= process_id < 399 == IsNCRes
	auto process_id = NuHepMC::ER3::ReadProcessID(evt);
	if( process_id >=200 && process_id <= 299 ) IsQel=true;
	else if ( process_id >=300 &&process_id <= 399 ) IsMec=true;
	else if ( process_id >=400 &&process_id <= 499 ) IsRes=true;
	else if ( process_id >=500 &&process_id <= 699 ) IsDis=true;
	
	// If missing it gets it from the GENIE format directly
	if( process_id == 0 && is_GENIE ) { 
	  int scatt_type = NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.ScatteringType");
	  if( scatt_type == 1 ) IsQel = true ; 
	  else if( scatt_type == 2 ) IsSingleK = true ; 
	  else if( scatt_type == 3 ) IsDis = true ; 
	  else if( scatt_type == 4 ) IsRes = true ; 
	  else if( scatt_type == 5 ) IsCoh = true ;
	  else if( scatt_type == 6 ) IsDfr = true ;
	  else if( scatt_type == 7 ) IsNuEL = true ;
	  else if( scatt_type == 8 ) IsImd = true ;
	  else if( scatt_type == 9 ) IsAMNuGamma = true ;
	  else if( scatt_type == 10 ) IsMec = true ;
	  else if( scatt_type == 11 ) IsCharmPro = true ;
	  
	  XSec = = NuHepMC::CheckedAttributeValue<double>(&evt,"GENIE.XSec");
	  DiffXSec = = NuHepMC::CheckedAttributeValue<double>(&evt,"GENIE.DiffXSec");
	}
      }
      
      auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);
      Target = tgtpt->pid();
      if( Target == 2212 ) Target = 1000010010 ;
      TargetZ=(Target/10000) - 1000*(Target/10000000);
      TargetA=(Target/10) - 1000*(Target/10000);;

      // Get electrons with id radcorr
      int id_radcorr = 123;
      auto corr_leptons = NuHepMC::Event::GetParticles_All(evt, id_radcorr, {beampt->pid()} );	
      Evcorr = corr_leptons[0]->momentum().e();
      Pxvcorr = corr_leptons[0]->momentum().px();
      Pyvcorr = corr_leptons[0]->momentum().py();
      Pzvcorr = corr_leptons[0]->momentum().pz();

      Elcorr = corr_leptons[1]->momentum().e();
      Pxlcorr = corr_leptons[1]->momentum().px();
      Pylcorr = corr_leptons[1]->momentum().py();
      Pzlcorr = corr_leptons[1]->momentum().pz();
      Plcorr = corr_leptons[1]->momentum().p3mod();
      Plcorr = cos( corr_leptons[1]->momentum().theta());
      
      // Get sturck nucleon information
      for(auto & part : evt.particles()){
	if(part->status() == NuHepMC::ParticleStatus::StruckNucleon) { 
	  //this is the struck nucleon
	  En = part->momentum().e();
	  Pxn = part->momentum().px() ; 
	  Pyn = part->momentum().py() ; 
	  Pzn = part->momentum().pz() ; 
	  HitNuc=part->pid();
	} else if (part->status() == NuHepMC::ParticleStatus::UndecayedPhysical || part->status() == NuHepMC::ParticleStatus::DecayedPhysical || part->status() == NuHepMC::ParticleStatus::IncomingBeam || part->status() == NuHepMC::ParticleStatus::Target ) { continue ; }
	break;
      }

      // Storing particles produced besides the outgoing lepton
      std::vector<HepMC3::ConstGenParticlePtr> final_parts = NuHepMC::Event::GetParticles_AllRealFinalState(evt,{});

      Nf = final_parts.size() ;
      NfEM = 0 ;
      for( unsigned int i = 0 ; i < Nf ; ++i ) { 
	Pdgf[i] = final_parts[i]->pid();
	if( Pdgf[i] != beampt->pid() ) { 
	  if( Pdgf[i] == kPdgProton ) ++NfP ; 
	  else if ( Pdgf[i] == kPdgNeutron ) ++NfN;
	  else if ( Pdgf[i] == kPdgPiP ) ++NfPip ; 
	  else if ( Pdgf[i] == kPdgPiM ) ++NfPim ; 
	  else if ( Pdgf[i] == kPdgPi0 ) ++NfPi0 ;
	  else if ( Pdgf[i] == kPdgKP ) ++NfKp ;
	  else if ( Pdgf[i] == kPdgKM ) ++NfKm ;
	  else if ( Pdgf[i] == kPdgK0 ) ++NfK0 ; 
	  else if ( Pdgf[i] == kPdgPhoton ) ++NfEM ; 
	  else ++NfOther;

	  Ef[i] = final_parts[i]->momentum().e();
	  Pxf[i] = final_parts[i]->momentum().px();
	  Pyf[i] = final_parts[i]->momentum().py();
	  Pzf[i] = final_parts[i]->momentum().pz();
	  Pf[i] = final_parts[i]->momentum().p3mod();
	  Costhf[i] = cos(final_parts[i]->momentum().theta());
	  //	  if( is_GENIE ) { 
	  //  Resc[i]=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.RescatterCode");
	  // }
	}
      }
      
      static const double kProtonMass = 0.9382720813 ;
      static const double kNeutronMass = 0.939565 ;
      double M = (kProtonMass+kNeutronMass)/2.;
      // Compute kinematics with true vertex kinematics
      // We use the leptons at the vertex to compute the variable W
      // It corresponds to the true W used for event generation
      auto q = corr_leptons[0]->momentum() - corr_leptons[1]->momentum();
      KineQ2 = -q.m2();
      KineX = 0.5*KineQ2/(M+q.e()) ;
      KineY = q.e()/beampt->momentum().e() ; 
      KineW = sqrt( pow(M,2) + 2*M*q.e() - KineQ2 ) ; 
      
      /*
      if( is_GENIE ) { 
	// True event information
	HitQrk=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.HitQuarkPDG");
	FromSea=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.HitSeaQuark");
	ResId=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.Resonance");
	XSec = NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.XSec") ;
	DXSec = NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.DiffXSec") ; 
	KineQ2=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.Q2");
	KineW=NuHepMC::CheckedAttributeValue<int>(&evt, "GENIE.Interaction.W");

	//	std::cout << KineQ2 - -q.m2()<<std::endl;
	// Not implemented: 
	VtxX = 0 ; 
	VtxY = 0 ; 
	VtxZ = 0 ; 
	SumKEf = 0 ; 
	CalResp0 = 0 ; 
	KPS = 0 ; 
      }
      */      
      Weight = evt.weights()[0]; // get event weight
 
      output_tree->Fill();
      return ;
    }
    
    HepMC3::GenParticlePtr GetPartFromId(HepMC3::GenEvent &evt, int id){
      for(auto &part : evt.particles()){
	if(part->id() == id){ return part; }
      }
      return nullptr;
    }

    std::string GetArg(std::string op, int argc, char ** argv )
      {
	const int buf_size = 2048*128;
	char *  argument   = new char[buf_size];
	strcpy(argument, "");

	while(argc>2)
	  {
	    if (argv[1][0] == '-' && argv[1][1] == '-') {

	      char op_cur[buf_size];
	      strcpy(op_cur,&argv[1][2]);

	      if (strcmp(op.c_str(),op_cur)==0) {
		if (strlen(&argv[2][0]) ) {
		  strcpy(argument,&argv[2][0]);
		}
	      }
	    }
	    argc--;
	    argv++;

	  }

	std::string value = std::string(argument);
	delete [] argument;
	return value ;
      }


    bool ExistArg(std::string op, int argc, char ** argv )
    {
      const int buf_size = 2048*128;
      char *  argument   = new char[buf_size];
      strcpy(argument, "");

      while(argc>2)
	{
	  if (argv[1][0] == '-' && argv[1][1] == '-') {

	    char op_cur[buf_size];
	    strcpy(op_cur,&argv[1][2]);

	    if (strcmp(op.c_str(),op_cur)==0) {
	      return true ;
	    }
	  }
	  argc--;
	  argv++;

	}
      delete [] argument ;
      return false;
    }

    double GetCLAS6TargetThickness( double tgt ) { 
      static const unsigned int kPdgHe3 = 1000020030; 
      static const unsigned int kPdgC12 = 1000060120 ; 
      static const unsigned int kPdgFe56 = 1000260560 ;

      const double tH = 0.027760 ;
      const double tHe = 0.005772;
      const double tC = 0.004183;
      const double tFe = 0.008532 ;
  
      if ( tgt == kPdgHe3 ) return tHe ; 
      else if ( tgt == kPdgC12 ) return tC ; 
      else if ( tgt == kPdgFe56 ) return tFe ; 
      return tH ; 
    }

  }
}

#endif 
