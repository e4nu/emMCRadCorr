// _________________________________________________________________________________
/* This app is used to generate the radiated flux to be used for event generation */
// _________________________________________________________________________________
#include <iostream>
#include <vector>
#include "RadiativeCorrUtils.h"
#include "Utils.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
// --output-file : name and path of file where to store output //
// --ebeam : beam energy of your experiment                    //
// --target : target pdg                                       //
// --thickness : thickness of your experiment target           //
// --resolution : electron energy resolution (bin width)       //
// --Emin : minimum energy for your histogram axis             //
// --Emax : maximum energy for your histogram axis             //
// --rad-model : model used for external radiation of in e-    //
//               simc or simple                                //
// --max-Ephoton : defaulted to 0.2 (of beam energy)           //
//                                                             //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  cout << "Generating radiated flux..." << endl;
  string output_file = "radiatedflux.root";
  double EBeam = 1 ; 
  int tgt = 1000060120 ;
  double resolution = 0.01;
  double Emin = 0.75 ;
  double Emax = EBeam+0.02 ;
  double thickness = utils::GetCLAS6TargetThickness(tgt); // Defaulted to CLAS6
  string rad_model = "simc";
  double MaxEPhoton = 0.2 ;

  if( argc > 1 ) { // configure rest of analysis
    if( utils::ExistArg("output-file",argc,argv)) {
      output_file = utils::GetArg("output-file",argc,argv); 
    }
    if( utils::ExistArg("rad-model",argc,argv)) {
      rad_model = utils::GetArg("rad-model",argc,argv); 
    }
    if( utils::ExistArg("ebeam",argc,argv)) {
      EBeam = stod(utils::GetArg("ebeam",argc,argv)); 
      Emax = EBeam+0.02 ; 
    }
    if( utils::ExistArg("target",argc,argv)) {
      tgt = stoi(utils::GetArg("target",argc,argv)); 
      thickness = utils::GetCLAS6TargetThickness(tgt); // Defaulted to CLAS6
    }
    if( utils::ExistArg("thickness",argc,argv)) {
      thickness = stoi(utils::GetArg("thickness",argc,argv)); 
    }
    if( utils::ExistArg("resolution",argc,argv)) {
      resolution = stod(utils::GetArg("resolution",argc,argv)); 
    }
    if( utils::ExistArg("Emin",argc,argv)) {
      Emin = stod(utils::GetArg("Emin",argc,argv)); 
    }
    if( utils::ExistArg("Emax",argc,argv)) {
      Emax = stod(utils::GetArg("Emax",argc,argv)); 
    }
    if( utils::ExistArg("max-Ephoton",argc,argv)) {
      MaxEPhoton = stod(utils::GetArg("max-Ephoton",argc,argv)); 
    }
  }
  // Calculate number of bins given resolution
  if( resolution == 0 ) return 0 ; 
  int nbins = (Emax - Emin)/resolution ; 

  std::unique_ptr<TFile> myFile( TFile::Open(output_file.c_str(), "RECREATE") );
  TH1D * hradflux = new TH1D( "hradflux", "Radiated Flux", nbins, Emin, Emax) ;   

  HepMC3::FourVector V4_beam(0,0,EBeam,EBeam);
  unsigned int nentries = 100000; 
  for( unsigned int i = 0 ; i < nentries ; ++i ) { 
    double egamma = SIMCEnergyLoss( V4_beam, 11, tgt, thickness, MaxEPhoton ) ;
    double Ee = EBeam - egamma ;
    hradflux -> Fill( Ee ) ; 
  }
  hradflux->Scale(1./hradflux->GetEntries());
  myFile->WriteObject(hradflux,"hradflux");

  std::cout << " Flux generated with " << rad_model << " with a resolution of " << resolution << " was stored in " << output_file << std::endl;

  return 0 ;
}
