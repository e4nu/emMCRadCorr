// __________________________________________________________________________
/* This app is used to plot quantities related to radiative corrections    */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <TMath.h>
#include <TRandom3.h>

/////////////////////////////////////////////////////////////////////////////////////////////
// Files can be accessed from the following link                                           //
// https://drive.google.com/drive/folders/1kkkKo8SUmdNYJYdm4Vg8XkD0yNne7Hej?usp=share_link //
/////////////////////////////////////////////////////////////////////////////////////////////

int plot_radFactor_Hydrogen( ) {
  bool plot_diff = true;
  TCanvas * c1 = new TCanvas("c1","c1",800,600);

  // This file is a SIMC simulation without radiative effects
  TFile * simc_unrad = TFile::Open("cafe_heep_singles_kin0_norad.root");
  if (!simc_unrad || simc_unrad->IsZombie()) {
      std::cerr << "Error opening file!" << std::endl;
      return 1;
  }

  // Get the tree from the file
  TTree *tree_simc_unrad = (TTree*)simc_unrad->Get("SNT");
  if (!tree_simc_unrad) {
      std::cerr << "Error: Tree not found!" << std::endl;
      simc_unrad->Close();
      return 1;
  }

  // This one is radiated
  TFile * simc_rad = TFile::Open("cafe_heep_singles_kin0_rad.root");
  if (!simc_rad || simc_rad->IsZombie()) {
      std::cerr << "Error opening file!" << std::endl;
      return 1;
  }

  // Get the tree from the file
  TTree *treerad = (TTree*)simc_rad->Get("SNT");
  if (!treerad) {
      std::cerr << "Error: Tree not found!" << std::endl;
      simc_rad->Close();
      return 1;
  }


  // This is the genie simulation
  // Open GENIE GST FILES
  TFile * genie_born = TFile::Open("GENIE_Heep_unradiated.gst.root");
  TTree * gst  = (TTree*)genie_born->Get("gst");
  if (!gst ) {
      std::cerr << "Error: Tree not found!" << std::endl;
      return 1;
  }

  TFile * genie_rad2 = TFile::Open("GENIE_Heep_radiated.gst.root");
  TTree * gst_rad2  = (TTree*)genie_rad2->Get("gst");
  if (!gst_rad2 ) {
      std::cerr << "Error: Tree not found!" << std::endl;
      return 1;
  }

  // This contains data from cafe
  TFile * data_file = TFile::Open("cafe_sample_LH2_heep_singles_16968_5000000_skimmed.root");
  if (!data_file || data_file->IsZombie()) {
      std::cerr << "Error opening file!" << std::endl;
      return 1;
  }
  TTree * data_tree = (TTree*) data_file->Get("T");


  double max = 3;
  int bins = 400;
  //Define histogram
  TH1D * hist_simc_unrad = new TH1D("hist_simc_unrad","",bins,0.9,max);
  TH1D * histRad = new TH1D("HistRAD","",bins,0.9,max);
  TH1D * histGENIE = new TH1D("HistGENIE","",bins,0.9,max);
  TH1D * histGENIERad2 = new TH1D("HistGENIERAD2","",bins,0.9,max);
  TH1D * data = new TH1D("data","",bins,0.9,max);

  // Variables to hold the data
  double W, Weight, theta_e, sig;
  // Set branch addresses
  tree_simc_unrad->SetBranchAddress("W", &W);
  tree_simc_unrad->SetBranchAddress("Weight", &Weight);
  tree_simc_unrad->SetBranchAddress("theta_e", &theta_e);
  tree_simc_unrad->SetBranchAddress("sig", &sig);

  double angle = 16;
  // Loop over all entries in the tree
  Long64_t nentries = tree_simc_unrad->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
      tree_simc_unrad->GetEntry(i);
      if( theta_e < 0.13 ) continue ;
      if( theta_e > 0.17 ) continue ;
      //if( theta != angle ) continue ;
      hist_simc_unrad->Fill(W,Weight);
  }

  treerad->SetBranchAddress("W", &W);
  treerad->SetBranchAddress("Weight", &Weight);
  treerad->SetBranchAddress("theta_e", &theta_e);
  nentries = treerad->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
      treerad->GetEntry(i);
      if( theta_e < 0.13 ) continue ;
      if( theta_e > 0.17 ) continue ;
      //if( theta != angle ) continue ;
      histRad->Fill(W,Weight);
  }

  // Variables to hold the data
  double W_genie ;
  double cthl;
  double Q2_genie;
  double pxv, pyv, pzv, Ev;
  double pxl, pyl, pzl, El;
  bool qel;
  // particle mass values
  static const double kElectronMass = 0.000510998;
  static const double kProtonMass = 0.9382720813;
  static const double kNeutronMass = 0.939565;
  static const double kNucleonMass = (kProtonMass + kNeutronMass) / 2;
  static const double kPiMMass = 0.139570;
  double M = (kProtonMass+kNeutronMass)/2.;
  double wght;
  // Set branch addresses
  gst->SetBranchAddress("Ev", &Ev);
  gst->SetBranchAddress("El", &El);
  gst->SetBranchAddress("cthl", &cthl);
  gst->SetBranchAddress("Q2", &Q2_genie);
  gst->SetBranchAddress("W", &W_genie);
  gst->SetBranchAddress("pxl", &pxl);
  gst->SetBranchAddress("pyl", &pyl);
  gst->SetBranchAddress("pzl", &pzl);
  gst->SetBranchAddress("El", &El);
  gst->SetBranchAddress("Ev", &Ev);
  gst->SetBranchAddress("pxv", &pxv);
  gst->SetBranchAddress("pyv", &pyv);
  gst->SetBranchAddress("pzv", &pzv);
  gst->SetBranchAddress("qel", &qel);
  gst->SetBranchAddress("wght", &wght);
  // We want to smear the electron momentum by its resolution
  double resolution = 0.001;//068;
  double Wshift = 0.003;
  // Loop over all entries in the tree
  nentries = gst->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
      gst->GetEntry(i);
      TLorentzVector beam ( pxv,pyv,pzv,Ev );
      TLorentzVector outl ( pxl,pyl,pzl,El );

      // we smear the outgoing electron only
      double outl_true_p = outl.P();
      double SmearedPOut = gRandom->Gaus(outl.P(),resolution*outl.P());
      double SmearedEOut = sqrt( pow( SmearedPOut,2 ) + pow( kElectronMass,2 ) ) ;
      outl.SetPxPyPzE( SmearedPOut/outl_true_p * outl.Px(), SmearedPOut/outl_true_p * outl.Py(), SmearedPOut/outl_true_p * outl.Pz(), SmearedEOut ) ;

      TLorentzVector q = (beam - outl) ;
      double Q2_reco = -q.M2();
      double W_reco = sqrt(pow(M,2) + 2*M*q.E() - Q2_reco );
      if( !qel ) continue ;
      if( cthl > TMath::Cos(0.13) ) continue ;
      if( cthl < TMath::Cos(0.17) ) continue ;
      //if( Q2_genie < 1.5 ) continue ;
      histGENIE->Fill(W_reco+Wshift,wght);
  }
  std::cout << " nentries def " << nentries<<std::endl;
  gst_rad2->SetBranchAddress("Ev", &Ev);
  gst_rad2->SetBranchAddress("El", &El);
  gst_rad2->SetBranchAddress("cthl", &cthl);
  gst_rad2->SetBranchAddress("Q2", &Q2_genie);
  gst_rad2->SetBranchAddress("W", &W_genie);
  gst_rad2->SetBranchAddress("pxl", &pxl);
  gst_rad2->SetBranchAddress("pyl", &pyl);
  gst_rad2->SetBranchAddress("pzl", &pzl);
  gst_rad2->SetBranchAddress("El", &El);
  gst_rad2->SetBranchAddress("Ev", &Ev);
  gst_rad2->SetBranchAddress("pxv", &pxv);
  gst_rad2->SetBranchAddress("pyv", &pyv);
  gst_rad2->SetBranchAddress("pzv", &pzv);
  gst_rad2->SetBranchAddress("qel", &qel);
  gst_rad2->SetBranchAddress("wght", &wght);
  //nentries = gst_rad2->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
      gst_rad2->GetEntry(i);
      TLorentzVector beam ( pxv,pyv,pzv,Ev );
      TLorentzVector outl ( pxl,pyl,pzl,El );

      double beam_true_p = beam.P();
      double SmearedPBeam = gRandom->Gaus(beam_true_p,resolution*beam_true_p);
      double SmearedEBeam = sqrt( pow( SmearedPBeam,2 ) + pow( kElectronMass,2 ) ) ;
      //beam.SetPxPyPzE( SmearedPBeam/beam_true_p * beam.Px(), SmearedPBeam/beam_true_p * beam.Py(), SmearedPBeam/beam_true_p * beam.Pz(), SmearedEBeam ) ;

      double outl_true_p = outl.P();
      double SmearedPOut = gRandom->Gaus(outl.P(),resolution*outl.P());
      double SmearedEOut = sqrt( pow( SmearedPOut,2 ) + pow( kElectronMass,2 ) ) ;
      outl.SetPxPyPzE( SmearedPOut/outl_true_p * outl.Px(), SmearedPOut/outl_true_p * outl.Py(), SmearedPOut/outl_true_p * outl.Pz(), SmearedEOut ) ;


      TLorentzVector q = (beam - outl) ;
      double Q2_reco = -q.M2();
      double W_reco = sqrt(pow(M,2) + 2*M*q.E() - Q2_reco );
      if( !qel || W_genie > 1 ) continue ;
      if( cthl > TMath::Cos(0.13) ) { continue ;}
      if( cthl < TMath::Cos(0.17) ) continue ;
      if(wght>0){
        histGENIERad2->Fill(W_reco+Wshift,wght);
      }
  }

  double W_data, theta_e_data ;
  data_tree->SetBranchAddress("P.kin.primary.W",&W_data);
  data_tree->SetBranchAddress("P.kin.primary.scat_ang_rad",&theta_e_data);
  nentries = data_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    data_tree->GetEntry(i);
    if( theta_e_data < 0.13 ) continue ;
    if( theta_e_data > 0.17 ) continue ;
    data->Fill(W_data-0.002);
  }
  data->Scale(1./data->GetEntries());

  histGENIE -> SetLineColor(kAzure+7 );
  histGENIERad2 -> SetLineColor(kAzure+7);
  hist_simc_unrad -> SetLineColor(kBlack);
  histRad -> SetLineColor(kBlack);
  hist_simc_unrad -> SetLineStyle(2);
  histGENIE -> SetLineStyle(2);
  histGENIERad2 -> SetLineWidth(4);
  histGENIE -> SetLineWidth(4);
  hist_simc_unrad -> SetLineWidth(4);
  histRad -> SetLineWidth(4);

  double simcintegral = hist_simc_unrad->Integral();
  double radintegral = histRad->Integral();
  hist_simc_unrad->Scale(1./simcintegral);
  double max_simc = hist_simc_unrad->GetMaximum() ;
  histRad->Scale(1./simcintegral);
  double histGENIE_Integral = histGENIE->Integral();
  double GENIEradintegral = histGENIERad2->Integral();
  histGENIE->Scale(1./histGENIE_Integral);
  double max_GENIE = histGENIE->GetMaximum() ;
  histGENIE->Scale(max_simc/max_GENIE);
  histGENIERad2->Scale(1./histGENIE_Integral*max_simc/max_GENIE);

  data->Sumw2();
  hist_simc_unrad->Sumw2();
  histRad->Sumw2();
  histGENIE->Sumw2();
  histGENIERad2->Sumw2();

  TH1D * hist_data_minus = (TH1D*)data->Clone();
  TH1D * hist_data = (TH1D*)data->Clone();
  if( plot_diff ){
    // substract data
    hist_data_minus->Scale(-1);
    hist_simc_unrad->Add(hist_data_minus);
    histRad->Add(hist_data_minus);
    histGENIE->Add(hist_data_minus);
    histGENIERad2->Add(hist_data_minus);
    data->Add(hist_data_minus);
    // Devide by original data
    hist_simc_unrad->Divide(hist_data);
    histRad->Divide(hist_data);
    histGENIE->Divide(hist_data);
    histGENIERad2->Divide(hist_data);
  }

  data -> SetMarkerStyle(8);
  data -> SetLineWidth(0);

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(data,"H(e,e')@10.6 GeV, Hall C data");
  legend->AddEntry(hist_simc_unrad,"SIMC, Unradiated");
  legend->AddEntry(histRad,"SIMC, Radiated");
  legend->AddEntry(histGENIE,"GENIE, Unradiated");
  legend->AddEntry(histGENIERad2,"GENIE, Radiated");

  if( !plot_diff ){
    hist_simc_unrad->GetXaxis()->SetRangeUser(0,1.1);
    hist_simc_unrad->GetXaxis()->SetTitle("W_{reco}[GeV]");
    hist_simc_unrad->GetYaxis()->SetTitle("Normalized event rate");
    hist_simc_unrad->Draw("hist err");
    histRad->Draw("hist same err");
    histGENIE->Draw("hist err same");
    histGENIERad2->Draw("hist err same");
    data->Draw("same err");
    legend->Draw();
  } else {
    data -> SetLineWidth(4);
    data -> SetMarkerSize(0);
    hist_simc_unrad->GetXaxis()->SetRangeUser(0,1.1);
    hist_simc_unrad->GetXaxis()->SetTitle("W_{reco}[GeV]");
    hist_simc_unrad->GetYaxis()->SetTitle("Pred-Data/Data");
    hist_simc_unrad->Draw("hist err");
    histRad->Draw("hist same err");
    histGENIE->Draw("hist err same");
    histGENIERad2->Draw("hist err same");
    data->Draw("same hist err");
    legend->Draw();
  }
  return 0 ;
}
