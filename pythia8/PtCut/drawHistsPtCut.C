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


int drawHistsPtCut() {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile("ptcut_hists_full.root", "READ");
  TH1D* ptProfile = ((TProfile*) inFile->Get( "hist bins" ))->ProjectionX(""); 

  TCanvas *canv = new TCanvas("c1","c1",1200,1200);
  setTDRStyle();
  canv->UseCurrentStyle();
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  canv->SetSelected(canv);
  canv->cd();
  // Show histogram  
  TPad *pad1 = new TPad("pad1","",0,0,1,1);

  pad1->Draw();
  pad1->cd();
  ptProfile->GetXaxis()->SetTitle("p_{T} (GeV)");
  ptProfile->GetYaxis()->SetTitle("simulation for detected pt/full pt");
  ptProfile->GetYaxis()->SetTitleOffset(1.6);
  ptProfile->SetStats(0);
  ptProfile->GetXaxis()->SetMoreLogLabels();
  ptProfile->GetXaxis()->SetNoExponent();
  ptProfile->Draw();

  std::cout << "\nDouble click on the histogram window to quit.\n";
  
  // Done.
  return 0;
}
