// This file uses the ptcut.root file to plot fractions of observed pT in a 
// detector simulation. 

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
#include "tdrstyle_mod3.C"


int drawProgressivePtCut(std::string fileName) {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile(fileName.c_str(), "READ");
  TH1D* ptProfile = ((TProfile*) inFile->Get( "pT bins" ))->ProjectionX(""); 
  TH1D* hcalProfile = ((TProfile*) inFile->Get( "hcal bins" ))->ProjectionX("");
  TH1D* gev3Profile = ((TProfile*) inFile->Get( "3gev bins" ))->ProjectionX("");
  TH1D* combProfile = (TH1D*) hcalProfile->Clone("hnew");
  combProfile->Multiply( gev3Profile );
  TH1D* hiEProfile = ((TProfile*) inFile->Get( "hie bins" ))->ProjectionX("");

  TCanvas *canv = new TCanvas("c1","c1",1200,1200);
  setTDRStyle();
  canv->UseCurrentStyle();
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  canv->SetSelected(canv);
  canv->cd();
  // Show histogram  
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  TPad *pad3 = new TPad("pad3","",0,0,1,1);
  TPad *pad4 = new TPad("pad4","",0,0,1,1);
  TPad *pad5 = new TPad("pad5","",0,0,1,1);

  pad2->SetFillStyle(4000);
  pad3->SetFillStyle(4000);
  pad4->SetFillStyle(4000);
  pad5->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  ptProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  ptProfile->GetYaxis()->SetTitle("simulation for detected pt/full pt");
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  ptProfile->SetStats(0);
  ptProfile->GetXaxis()->SetMoreLogLabels();
  ptProfile->GetXaxis()->SetNoExponent();
  ptProfile->Draw();
  TLegend *leg = tdrLeg(0.4,0.73,0.6,0.93);
  leg->AddEntry( ptProfile, "All effects" );
  leg->AddEntry( combProfile, "HCal+3GeV" );
  leg->AddEntry( hcalProfile, "HCal" );
  leg->AddEntry( gev3Profile, "3GeV" );
  leg->AddEntry( hiEProfile, "High pt ineff." );
  leg->Draw();
  pythiaFinal();
  canv->Modified();
  pad1->UseCurrentStyle();
  ptProfile->SetMarkerStyle(kCircle);
  ptProfile->SetMaximum(1.04);
  ptProfile->SetMinimum(0.89);
  ptProfile->SetMarkerColor(7);
  pad1->Update();
  
  canv->cd();
  pad2->Draw();
  pad2->cd();
  combProfile->SetStats(0);
  combProfile->GetXaxis()->SetMoreLogLabels();
  combProfile->GetXaxis()->SetNoExponent();
  combProfile->Draw();
  pad2->UseCurrentStyle();
  combProfile->SetMaximum(1.04);
  combProfile->SetMinimum(0.89);
  combProfile->SetMarkerColor(3);
  pad2->Update();

  canv->cd();
  pad3->Draw();
  pad3->cd();
  hcalProfile->SetStats(0);
  hcalProfile->GetXaxis()->SetMoreLogLabels();
  hcalProfile->GetXaxis()->SetNoExponent();
  hcalProfile->Draw();
  pad3->UseCurrentStyle();
  hcalProfile->SetMaximum(1.04);
  hcalProfile->SetMinimum(0.89);
  hcalProfile->SetMarkerStyle(kStar);
  hcalProfile->SetMarkerColor(5);
  pad3->Update();
  
  canv->cd();
  pad4->Draw();
  pad4->cd();
  gev3Profile->SetStats(0);
  gev3Profile->GetXaxis()->SetMoreLogLabels();
  gev3Profile->GetXaxis()->SetNoExponent();
  gev3Profile->Draw();
  pad4->UseCurrentStyle();
  gev3Profile->SetMaximum(1.04);
  gev3Profile->SetMinimum(0.89);
  gev3Profile->SetMarkerStyle(kMultiply);
  pad4->Update();

  canv->cd();
  pad5->Draw();
  pad5->cd();
  hiEProfile->SetStats(0);
  hiEProfile->GetXaxis()->SetMoreLogLabels();
  hiEProfile->GetXaxis()->SetNoExponent();
  hiEProfile->Draw();
  pad5->UseCurrentStyle();
  hiEProfile->SetMaximum(1.04);
  hiEProfile->SetMinimum(0.89);
  hiEProfile->SetMarkerStyle(kCircle);
  hiEProfile->SetMarkerColor(kRed);
  pad5->Update();

  std::cout << "\nDouble click on the histogram window to quit.\n";
  
  // Done.
  return 0;
}
