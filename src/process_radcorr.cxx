// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include "NuHepMC/HepMC3Features.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/WriterUtils.hxx"
#include "NuHepMC/make_writer.hxx"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "RadiativeCorrUtils.h"
#include "Utils.h"

///////////////////////////////////////////////////////////////////////////////
// process_radcorr.cxx input_file outputfile.hepmc3                          //
///////////////////////////////////////////////////////////////////////////////
using namespace e4nu;
using namespace utils;

int main(int, char const *argv[]) {

  auto rdr = HepMC3::deduce_reader(argv[1]);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << argv[1]
              << std::endl;
    return 1;
  }

  HepMC3::GenEvent evt;
  rdr->read_event(evt);
  if (rdr->failed()) {
    std::cout << "Failed to read first event from " << argv[1] << "."
              << std::endl;
    return 1;
  }
  
  TFile * output_gst = new TFile("genie.gst.root","RECREATE");
  TTree * output_tree = new TTree("gst","GENIE Summary Event Tree");
  SetGSTBranchAddress( output_tree );

  auto in_gen_run_info = evt.run_info();
  auto vtx_statuses = NuHepMC::GR5::ReadVertexStatusIdDefinitions(in_gen_run_info);
  auto part_statuses = NuHepMC::GR6::ReadParticleStatusIdDefinitions(in_gen_run_info);
  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // Define ID for radiative corrections that is not used
  int MyRadVertexStatus = 123;
  part_statuses[MyRadVertexStatus] = {"rad_lepton", "Radiated corrections"};

  out_gen_run_info->tools().push_back(HepMC3::GenRunInfo::ToolInfo{ "emMCRadCorr", "version 1", "Adding radiative corrections to EM interactions"});
  
  NuHepMC::GR5::WriteVertexStatusIDDefinitions(out_gen_run_info, vtx_statuses);
  NuHepMC::GR6::WriteParticleStatusIDDefinitions(out_gen_run_info, part_statuses);
  
  // add link to your paper describing this model to the citation metadata
  NuHepMC::GC6::AddGeneratorCitation(out_gen_run_info, "arxiv", {"2404.12345v3",});
  
  auto wrtr = std::unique_ptr<HepMC3::Writer>(NuHepMC::Writer::make_writer(argv[2], out_gen_run_info));
  
  // re-open the file so that you start at the beginning
  rdr = HepMC3::deduce_reader(argv[1]);
  size_t nprocessed = 0;
  while (true) { // loop while there are events
    
    rdr->read_event(evt);
    if (rdr->failed()) {
      std::cout << "Reached the end of the file after " << nprocessed
                << " events." << std::endl;
      break;
    }
    evt.set_run_info(out_gen_run_info);
    evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

    auto beampt = NuHepMC::Event::GetBeamParticle(evt);
    auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);

    if (!beampt || !tgtpt) {  // this event didn't have a beam particle or
                              // target, its an odd one
      wrtr->write_event(evt); // write out events that we don't modify
      continue;
    }
    
    // Events are generated with the radiated flux
    // The original beam is a monocromatic beam of a specific energy, EBeam
    // The correct calculation requires to add back true beam in the event rate
    auto beampt_corr = GetPartFromId(evt, beampt->id());
    beampt_corr->set_status( MyRadVertexStatus ); // This corresponds to the radiated electron (radiated corrected electron)
    
    // Store generated photon
    auto beampt_mono = std::make_shared<HepMC3::GenParticle>(beampt_corr->data()); // Monocromatic beam
    const HepMC3::FourVector true_beam ( 0,0,4.325,4.325); // from configuration in z direction
    beampt_mono->set_pid(beampt_corr->pid());
    beampt_mono->set_momentum( true_beam );
    beampt_mono->set_status( NuHepMC::ParticleStatus::IncomingBeam ) ; // This is the true beam 

    // We compute the kinematics of the emited photon using the information from the beam and radiated electron:
    auto beam_photon = std::make_shared<HepMC3::GenParticle>(beampt_corr->data());
    beam_photon->set_momentum( true_beam - beampt->momentum() ) ; 
    beam_photon->set_pid( kPdgPhoton ) ;
    beam_photon->set_status( NuHepMC::ParticleStatus::UndecayedPhysical ) ;
    
    auto primary_vtx = NuHepMC::Event::GetPrimaryVertex(evt);
    auto emfslep_pid = beampt->pid();
    
    // From the primary vertex, we can get the outgoing electron kinematics before it radiates:
    auto primary_leptons = NuHepMC::Vertex::GetParticlesOut_All(primary_vtx, NuHepMC::ParticleStatus::UndecayedPhysical, {emfslep_pid} );
    if ( primary_leptons.size()!=1 ) { // this event had no primary leptons.
      std::cout << "More than one outgoing electron exits. Not stored.."<<std::endl;
      continue;
    }

    auto fslep = primary_leptons.back();
    auto fslep_corr = GetPartFromId(evt, fslep->id());
    // As it is the outgoing electron at the vertex (before radiation), we give it a different status
    fslep_corr->set_status( MyRadVertexStatus ); 

    // Now we account for the fact that the outgoing electron might also radiate
    // Modify outgoing electron kinematics
    auto fslep_detected = std::make_shared<HepMC3::GenParticle>(fslep_corr->data());
    fslep_detected->set_status(NuHepMC::ParticleStatus::UndecayedPhysical) ;

    // Store generated photon and alter out.electron kinematics
    auto out_photon = std::make_shared<HepMC3::GenParticle>(fslep_corr->data());
    out_photon->set_status(NuHepMC::ParticleStatus::UndecayedPhysical) ;
    out_photon->set_pid(kPdgPhoton);

    // This should all be configurable
    double thickness = utils::GetCLAS6TargetThickness(tgtpt->pid());
    double Emin = 0.75 ;
    double Emax = beampt->momentum().e()+0.02 ;
    double max_Ephoton = 0.2*beampt_mono->momentum().e();
    string model = "simc";//"vanderhaeghen";//"motsai";//"vanderhaeghen"; 

    // Compute true detected outgoing electron kinematics with energy loss method
    double egamma = utils::SIMCEnergyLoss( fslep_corr->momentum(), fslep_corr->pid(), tgtpt->pid(), thickness, max_Ephoton ) ;
    if( egamma < 0 ) egamma = 0 ;
    HepMC3::FourVector OutGamma = utils::GetEmittedHardPhoton( fslep_corr->momentum(), egamma ) ;
    if( OutGamma.e() < 0 )  OutGamma.set(0,0,0,0);
    out_photon->set_momentum(OutGamma);
    fslep_detected->set_momentum( fslep_corr->momentum() - OutGamma ) ;

    // Add all particles to event record
    evt.add_particle(beampt_mono);
    evt.add_particle(fslep_detected);
    evt.add_particle(beam_photon);
    evt.add_particle(out_photon);

    // For the weight calculation, we need the true Q2 used for event generation
    // We compute it with vertex kinematics
    auto q = beampt_corr->momentum() - fslep_corr->momentum();
    double vertex_Q2 = -q.m2();

    // Alter event weigth to account for vertex and vacumm effects
    evt.weights()[0] = utils::RadCorrWeight( evt, vertex_Q2, thickness, max_Ephoton, "simc" );

    // Store back in hepmc3 format:
    wrtr->write_event(evt); 

    // Store in GENIE gst output
    StoreHepMCToGST( evt, output_tree );
    ++nprocessed;
  }
  wrtr->close();
  output_tree->Write();
  output_gst->Close();
}
