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
#include "Utils.h"

///////////////////////////////////////////////////////////////////////////////
// ./processor_skeleton input_file outputfile                                //
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

  // for files that you know can be opened multiple times without failing (i.e.
  // not a stream) it is easiest to read one event, and grab the
  // HepMC3::GenRunInfo here.
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

  // modify gen_run_info here to add to the file provenance that your code has
  // run on the file

  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // Define ID for radiative corrections that is not used
  int MyRadVertexStatus = 123;
  while ( vtx_statuses.count(MyRadVertexStatus) > 0 ) {
    ++MyRadVertexStatus ; 
  }
  vtx_statuses[MyRadVertexStatus] = {"rad_vcorr", "Radiated corrections"};
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
    // ensure units are in MeV will perform conversions if the previous
    // generation step used GeV
    evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

    auto beampt = NuHepMC::Event::GetBeamParticle(evt);
    auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);

    if (!beampt || !tgtpt) {  // this event didn't have a beam particle or
                              // target, its an odd one
      wrtr->write_event(evt); // write out events that we don't modify
      continue;
    }


    // Add back true beam
    auto rad_beam = evt.particles()[beampt->id()];
    rad_beam->set_status( MyRadVertexStatus );

    // Store generated photon
    auto beampt_preRad = std::make_shared<HepMC3::GenParticle>(rad_beam->data());
    const HepMC3::FourVector true_beam ( 0,0,4.325,4.325); // from configuration
    beampt_preRad->set_momentum( true_beam );
    beampt_preRad->set_status( NuHepMC::ParticleStatus::IncomingBeam ) ; // This is the true beam 
    
    auto beam_photon = std::make_shared<HepMC3::GenParticle>(rad_beam->data());
    beam_photon->set_momentum( true_beam - beampt->momentum() ) ; 
    beam_photon->set_status( NuHepMC::ParticleStatus::UndecayedPhysical ) ;

    auto primary_vtx = NuHepMC::Event::GetPrimaryVertex(evt);
    auto emfslep_pid = beampt->pid();
    
    // grab all primary leptons that were considered 'final state' by the
    // previous simulation
    //  might have to adjust for simulations that include lepton FSI already
    auto primary_leptons = NuHepMC::Vertex::GetParticlesOut_All(primary_vtx, NuHepMC::ParticleStatus::UndecayedPhysical, {emfslep_pid} );

    if ( primary_leptons.size()!=1 ) { // this event had no primary leptons.
      wrtr->write_event(evt);      // write out events that we don't modify
      continue;
    }

    auto fslep = primary_leptons.back();
    auto fslep_preRad = evt.particles()[fslep->id()];
    fslep_preRad->set_status( MyRadVertexStatus );

    // Store generated photon
    auto out_photon = std::make_shared<HepMC3::GenParticle>(fslep_preRad->data());
    HepMC3::FourVector Delta_Photon ( 0,0,0.1,0.1); // random for now 
    out_photon->set_momentum(Delta_Photon);
    out_photon->set_status( NuHepMC::ParticleStatus::UndecayedPhysical ) ;

    // Modify outgoing electron kinematics
    auto fslep_postRad = std::make_shared<HepMC3::GenParticle>( fslep_preRad->data() ); 
    fslep_postRad->set_momentum(fslep_preRad->momentum()-Delta_Photon);
    fslep_postRad->set_status(NuHepMC::ParticleStatus::UndecayedPhysical) ; // detected electron

    // make a new vertex to represent the radiated event
    auto lepRadvtx = std::make_shared<HepMC3::GenVertex>();
    lepRadvtx->set_status(MyRadVertexStatus);
    evt.add_vertex(lepRadvtx); 
    lepRadvtx->add_particle_in(fslep_preRad);
    lepRadvtx->add_particle_out(fslep_postRad);
    lepRadvtx->add_particle_out(out_photon);

    wrtr->write_event(evt); // write out events your modified event

    // Store in gst output
    StoreHepMCToGST( evt, output_tree );
    ++nprocessed;
    if( nprocessed == 40 ) break;
  }
  wrtr->close();
  output_tree->Write();
  output_gst->Close();
}
