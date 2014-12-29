// plots jet energy fractions, see README_ScriptInfo

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// ROOT, for histogramming.
#include "TROOT.h"
#include "TF1.h"
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
#include "tdrstyle_mod1.C"

int runSortedJets() {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile("sortedjets.root", "READ");
  THStack* partFracs = new THStack("particlestacks", "");
  THStack* errorFracs = new THStack("errorstacks", "" );
  
  // Read histograms.
  vector<TProfile*> fractionProfiles;
  for (int idx = 0; idx != 16; ++idx){
    std::stringstream tmpString("");
    tmpString << "a" << idx; 
    // a = all, g = gluonjets, q = quarkjets, lq = low pt quarks, hq = high pt quarks
    fractionProfiles.push_back( (TProfile*) inFile->Get( tmpString.str().c_str()) );
  }
  cout << fractionProfiles.size();

  // Create assosiated histograms for the profiles
  vector<TH1D*> fractionHists;
  for (int idx = 0; idx != 16; ++idx){  
    fractionHists.push_back( fractionProfiles[idx]->ProjectionX("","") );
  }

  // Set fill colors
  fractionHists[0]->SetFillColor(kYellow);
  fractionHists[1]->SetFillColor(kYellow+1);
  fractionHists[2]->SetFillColor(kYellow+2);
  fractionHists[3]->SetFillColor(kBlue);
  fractionHists[4]->SetFillColor(kBlue+1);
  fractionHists[5]->SetFillColor(kBlue+2);
  fractionHists[6]->SetFillColor(kBlue-2);
  fractionHists[7]->SetFillColor(kRed);
  fractionHists[8]->SetFillColor(kRed+1);
  fractionHists[9]->SetFillColor(kRed+2);
  fractionHists[10]->SetFillColor(kRed+3);
  fractionHists[11]->SetFillColor(kGreen);
  fractionHists[12]->SetFillColor(kGreen+1);
  fractionHists[13]->SetFillColor(kGreen+2);
  fractionHists[14]->SetFillColor(kGreen+3);
  fractionHists[15]->SetFillColor(kCyan);

  for (unsigned int i = 0; i != fractionHists.size(); ++i){
    partFracs->Add( fractionHists[i] );
    errorFracs->Add( fractionProfiles[i]->ProjectionX() );
  }

  TCanvas *canv = new TCanvas("c1","c1",600,600);
  setTDRStyle();
  canv->SetSelected(canv);
  canv->cd();
  TLegend *leg = new TLegend(0.84,0.1,0.98,0.95);
  leg->AddEntry( fractionHists[15], "Rest" );
  leg->AddEntry( fractionHists[14], "e^{+},e^{-},#mu^{+},#mu^{-} " );
  leg->AddEntry( fractionHists[13], "#Sigma^{+},#Sigma^{-}" );
  leg->AddEntry( fractionHists[12], "#Lambda^{0}" );
  leg->AddEntry( fractionHists[11], "#gamma (rest)" );
  leg->AddEntry( fractionHists[10], "#bar{n}^{0}" );
  leg->AddEntry( fractionHists[9], "n^{0}" );
  leg->AddEntry( fractionHists[8], "#bar{p}^{-}" );
  leg->AddEntry( fractionHists[7], "p^{+}" );
  leg->AddEntry( fractionHists[6], "K_{L}^{0}" );
  leg->AddEntry( fractionHists[5], "K_{S}^{0}" );
  leg->AddEntry( fractionHists[4], "K^{-}" );
  leg->AddEntry( fractionHists[3], "K^{+}" );
  leg->AddEntry( fractionHists[2], "#gamma (#pi^{0})" );
  leg->AddEntry( fractionHists[1], "#pi^{-}" );
  leg->AddEntry( fractionHists[0], "#pi^{+}" );
  // Show histogram  
  fractionHists[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  fractionHists[0]->GetYaxis()->SetTitle("Jet energy fraction (1)");
  fractionHists[0]->SetStats(0);
  fractionHists[0]->GetXaxis()->SetMoreLogLabels();
  fractionHists[0]->GetXaxis()->SetNoExponent();
  fractionHists[0]->GetYaxis()->SetTitleOffset(1.25);
  fractionHists[0]->GetXaxis()->SetTitleOffset(1.1);

  partFracs->SetHistogram( fractionHists[0] );
  errorFracs->SetHistogram( fractionHists[0] );
  partFracs->Draw();
  errorFracs->Draw("same");
  canv->UseCurrentStyle();
  pythiaFinal();
  leg->Draw();
  canv->Modified();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();

  // Save histogram on file and close file.
  canv->Write();
  
  // Done.
  return 0;
}
