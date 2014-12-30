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


int drawCaloPtCut() {
  double minH = 0.94; double maxH = 1.01;
  int caseFlag = 0;
  if (caseFlag == 0){
    minH = 0.9; maxH = 1;
  } else if (caseFlag == 1){
    minH = 0.9; maxH = 1.01;
  } else {
    minH = 0.91; maxH = 1.01;
  }
  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile("ptcut_hists_full.root", "READ");
  TH1D* normalDummy; TH1D* caloDummy;
  if (caseFlag == 0){
    normalDummy = ((TProfile*) inFile->Get( "pT bins" ))->ProjectionX("");  
    caloDummy = ((TProfile*) inFile->Get( "pT calo bins" ))->ProjectionX("");
  } else if (caseFlag == 1){
    normalDummy = ((TProfile*) inFile->Get( "hcal bins" ))->ProjectionX("");
    caloDummy = ((TProfile*) inFile->Get( "hcal calo bins" ))->ProjectionX(""); 
  } else {
    normalDummy = ((TProfile*) inFile->Get( "3gev bins" ))->ProjectionX(""); 
    caloDummy = ((TProfile*) inFile->Get( "3gev calo bins" ))->ProjectionX(""); 
  }
  TCanvas *canv = new TCanvas("c1","c1",1200,1200);
  setTDRStyle();
  canv->UseCurrentStyle();
  normalDummy->GetYaxis()->SetTitleOffset(1.6);
  canv->SetSelected(canv);
  canv->cd();
  // Show histogram  
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad1->cd();
  normalDummy->GetXaxis()->SetTitle("p_{T} (GeV)");
  normalDummy->GetYaxis()->SetTitle("simulation for detected pt/full pt");
  normalDummy->GetYaxis()->SetTitleOffset(1.6);
  normalDummy->SetStats(0);
  normalDummy->GetXaxis()->SetMoreLogLabels();
  normalDummy->GetXaxis()->SetNoExponent();
  normalDummy->Draw();
  TLegend *leg = tdrLeg(0.2,0.75,0.4,0.95);
  leg->AddEntry( normalDummy, "nh" );
  leg->AddEntry( caloDummy, "nh and ch" );
  leg->Draw();
  pythiaFinal();

  canv->Modified();
  pad1->UseCurrentStyle();
  normalDummy->SetMarkerStyle(kCircle);
  normalDummy->SetMaximum(maxH);
  normalDummy->SetMinimum(minH);
  normalDummy->SetMarkerColor(kBlack);
  pad1->Update();
  
  canv->cd();
  pad2->Draw();
  pad2->cd();
  caloDummy->SetStats(0);
  caloDummy->GetXaxis()->SetMoreLogLabels();
  caloDummy->GetXaxis()->SetNoExponent();
  caloDummy->Draw();
  pad2->UseCurrentStyle();
  caloDummy->SetMaximum(maxH);
  caloDummy->SetMinimum(minH);
  caloDummy->SetMarkerColor(kCyan);
  pad2->Update();

  std::cout << "\nDouble click on the histogram window to quit.\n";
  
  // Done.
  return 0;
}
