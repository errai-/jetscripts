// This file uses the ptcut.root file to plot fractions of observed pT in a 
// detector simulation. 

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "TF1.h"

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


int drawSamePtCut() {
  int oldFlag = 0;
  char oldName[] = "ptcut_olds.root";
  char newName[] = "ptcut_news.root";
  double minH = 0.94; double maxH = 1.01;
  if (oldFlag){
    maxH = 1.01; minH = 0.94;
  } else {
    maxH= 1.005; minH = 0.95;
  }  
  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile(oldFlag ? oldName : newName, "READ");
  TH1D* ptProfile = ((TProfile*) inFile->Get( "pT bins" ))->ProjectionX(""); 
  TH1D* ecalProfile; 
  if (oldFlag)
    ecalProfile = ((TProfile*) inFile->Get( "ecal bins" ))->ProjectionX("");
  TH1D* hcalProfile = ((TProfile*) inFile->Get( "hcal bins" ))->ProjectionX("");
  TH1D* gev3Profile = ((TProfile*) inFile->Get( "3gev bins" ))->ProjectionX("");
  TCanvas *canv = new TCanvas("c1","c1",1200,1200);
  setTDRStyle();
  canv->UseCurrentStyle();
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  canv->SetSelected(canv);
  canv->cd();
  // Show histogram  
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2;
  if (oldFlag)
    pad2 = new TPad("pad2","",0,0,1,1);
  TPad *pad3 = new TPad("pad3","",0,0,1,1);
  TPad *pad4 = new TPad("pad4","",0,0,1,1);
  if (oldFlag)
    pad2->SetFillStyle(4000);
  pad3->SetFillStyle(4000);
  pad4->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  ptProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  ptProfile->GetYaxis()->SetTitle("simulation for detected pt/full pt");
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  ptProfile->SetStats(0);
  ptProfile->GetXaxis()->SetMoreLogLabels();
  ptProfile->GetXaxis()->SetNoExponent();
  ptProfile->Draw();
  TLegend *leg = tdrLeg(0.2,0.7,0.4,0.9);
  leg->AddEntry( ptProfile, "total" );
  if (oldFlag)
    leg->AddEntry( ecalProfile, "ecal" );
  leg->AddEntry( hcalProfile, "hcal" );
  leg->AddEntry( gev3Profile, "3GeV" );
  leg->Draw();
  pythiaFinal();

  canv->Modified();
  pad1->UseCurrentStyle();
  ptProfile->SetMarkerStyle(kMultiply);
  ptProfile->SetMaximum(maxH);
  ptProfile->SetMinimum(minH);
  ptProfile->SetMarkerColor(kBlack);
  pad1->Update();
  
  if (oldFlag){
    canv->cd();
    pad2->Draw();
    pad2->cd();
    ecalProfile->SetStats(0);
    ecalProfile->GetXaxis()->SetMoreLogLabels();
    ecalProfile->GetXaxis()->SetNoExponent();
    ecalProfile->Draw();
    pad2->UseCurrentStyle();
    ecalProfile->SetMaximum(maxH);
    ecalProfile->SetMinimum(minH);
    ecalProfile->SetMarkerColor(kMagenta);
    pad2->Update();
  }

  canv->cd();
  pad3->Draw();
  pad3->cd();
  hcalProfile->SetStats(0);
  hcalProfile->GetXaxis()->SetMoreLogLabels();
  hcalProfile->GetXaxis()->SetNoExponent();
  hcalProfile->Draw();
  pad3->UseCurrentStyle();
  hcalProfile->SetMaximum(maxH);
  hcalProfile->SetMinimum(minH);
  hcalProfile->SetMarkerStyle(kStar);
  hcalProfile->SetMarkerColor(kGreen);
  pad3->Update();
  
  canv->cd();
  pad4->Draw();
  pad4->cd();
  gev3Profile->SetStats(0);
  gev3Profile->GetXaxis()->SetMoreLogLabels();
  gev3Profile->GetXaxis()->SetNoExponent();
  gev3Profile->Draw();
  pad4->UseCurrentStyle();
  gev3Profile->SetMaximum(maxH);
  gev3Profile->SetMinimum(minH);
  gev3Profile->SetMarkerStyle(kCircle);
  gev3Profile->SetMarkerColor(kCyan);
  pad4->Update();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  
  // Done.
  return 0;
}
