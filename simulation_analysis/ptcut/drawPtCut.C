// plots the effects of detector simulation, see readme

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "TF1.h"
#include <string>

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod12.C"


int drawPtCut(std::string fileName) {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile(fileName.c_str(), "READ");
  TProfile* ptProfile = (TProfile*) inFile->Get( "pT bins" ); 
  TProfile* ecalProfile = (TProfile*) inFile->Get( "ecal bins" );
  TProfile* hcalProfile = (TProfile*) inFile->Get( "hcal bins" );
  TProfile* gev3Profile = (TProfile*) inFile->Get( "3gev bins" );
  TCanvas *canv = new TCanvas("c1","c1",1200,1200);
  setTDRStyle();
  canv->UseCurrentStyle();
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  canv->SetSelected(canv);
  canv->Divide(2,2);
  canv->cd(1);
  // Show histogram  
  ptProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  ptProfile->GetYaxis()->SetTitle("simulation for detected pt/full pt (1)");
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  ptProfile->SetStats(0);
  ptProfile->GetXaxis()->SetMoreLogLabels();
  ptProfile->GetXaxis()->SetNoExponent();
  ptProfile->Draw();
  pythiaFinal();
  canv->cd(2);
  ecalProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  ecalProfile->GetYaxis()->SetTitle("pt with simulated false ecal photons/full pt (1)");
  ecalProfile->GetYaxis()->SetTitleOffset(1.6);
  ecalProfile->SetStats(0);
  ecalProfile->GetXaxis()->SetMoreLogLabels();
  ecalProfile->GetXaxis()->SetNoExponent();
  ecalProfile->Draw();
  pythiaFinal();
  canv->cd(3);
  hcalProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  hcalProfile->GetYaxis()->SetTitle("pt with simulated hcal response/full pt (1)");
  hcalProfile->GetYaxis()->SetTitleOffset(1.6);
  hcalProfile->SetStats(0);
  hcalProfile->GetXaxis()->SetMoreLogLabels();
  hcalProfile->GetXaxis()->SetNoExponent();
  hcalProfile->Draw();
  pythiaFinal();
  canv->cd(4);
  gev3Profile->GetXaxis()->SetTitle("p_{T} (GeV)");
  gev3Profile->GetYaxis()->SetTitle("pt with simulated 3GeV cut for neutral hadrons/full pt (1)");
  gev3Profile->GetYaxis()->SetTitleOffset(1.6);
  gev3Profile->SetStats(0);
  gev3Profile->GetXaxis()->SetMoreLogLabels();
  gev3Profile->GetXaxis()->SetNoExponent();
  gev3Profile->Draw();

  pythiaFinal();
  canv->Modified();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  
  // Done.
  return 0;
}
