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
#include "TColor.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TLatex.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod14_pythia8.C"


void stackModify(TH1D *setter){
  setter->GetXaxis()->SetTitle("p_{T} (GeV)");
  setter->GetYaxis()->SetTitle("Jet energy fraction");
  setter->SetStats(0);
  setter->GetXaxis()->SetMoreLogLabels();
  setter->GetXaxis()->SetNoExponent();
  setter->GetYaxis()->SetTitleOffset(1.25);
  setter->GetXaxis()->SetTitleOffset(1.0);
  setter->GetXaxis()->SetLabelSize(0.045);
  setter->GetYaxis()->SetLabelSize(0.045);
  setter->GetXaxis()->SetTitleSize(0.045);
  setter->GetYaxis()->SetTitleSize(0.045);
}

int runPFSortedJets() {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile("jet_storage.root", "READ");
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

  vector<int> order(16,0);
  order[1]=1; order[2]=3; order[3]=4; order[4]=13; order[5]=7; order[6]=8; order[7]=9;
  order[8]=10; order[9]=5; order[10]=6; order[11]=12; order[12]=2; order[13]=11;
  order[14]=14; order[15]=15;
  // Create assosiated histograms for the profiles
  vector<TH1D*> fractionHists;
  for (int idx = 0; idx != 16; ++idx){
    fractionHists.push_back( fractionProfiles[idx]->ProjectionX("","") );
  }

  // Set fill colors
  fractionHists[0]->SetFillColor(kRed-10);
  fractionHists[1]->SetFillColor(kRed-9);
  fractionHists[3]->SetFillColor(kRed-7);
  fractionHists[4]->SetFillColor(kRed-6);
  fractionHists[13]->SetFillColor(kMagenta-10);
  fractionHists[7]->SetFillColor(kMagenta-9);
  fractionHists[8]->SetFillColor(kMagenta-7);
  fractionHists[9]->SetFillColor(kGreen-10);
  fractionHists[10]->SetFillColor(kGreen-9);
  fractionHists[5]->SetFillColor(kGreen-7);
  fractionHists[6]->SetFillColor(kGreen-6);
  fractionHists[12]->SetFillColor(kTeal-9);
  fractionHists[2]->SetFillColor(kBlue-10);
  fractionHists[11]->SetFillColor(kBlue-9);
  fractionHists[14]->SetFillColor(kCyan-9);
  fractionHists[15]->SetFillColor(kOrange-9);

  for (unsigned int i = 0; i != fractionHists.size(); ++i){
    partFracs->Add( fractionHists[order[i]] );
    errorFracs->Add( fractionProfiles[order[i]]->ProjectionX() );
  }
 
  TCanvas *canv = tdrCanvas("c1",fractionHists[0],12,0,1);//new TCanvas("c1","c1",600,600);

  stackModify(fractionHists[0]);
  partFracs->SetHistogram( fractionHists[0] );
  errorFracs->SetHistogram( fractionHists[0] );
  partFracs->GetHistogram()->SetMaximum(1.01);
  partFracs->GetHistogram()->SetMinimum(-0.01);
  errorFracs->GetHistogram()->SetMaximum(1.0);
  errorFracs->GetHistogram()->SetMinimum(0.0);
  //partFracs->SetFillStyle(4000);

  TLegend *leg = tdrLeg(0.83,0.1,0.98,0.93);
  leg->AddEntry( fractionHists[15], "rest" );
  leg->AddEntry( fractionHists[14], "e,#mu" );
  leg->AddEntry( fractionHists[11], "#gamma (rest)" );
  leg->AddEntry( fractionHists[2], "#gamma (#pi^{0})" );
  leg->AddEntry( fractionHists[12], "#Lambda^{0}" );
  leg->AddEntry( fractionHists[6], "K_{L}^{0}" );
  leg->AddEntry( fractionHists[5], "K_{S}^{0}" );
  leg->AddEntry( fractionHists[10], "#bar{n}^{0}" );
  leg->AddEntry( fractionHists[9], "n^{0}" );
  leg->AddEntry( fractionHists[8], "#bar{p}^{-}" );
  leg->AddEntry( fractionHists[7], "p^{+}" );
  leg->AddEntry( fractionHists[13], "#Sigma^{#pm}" );
  leg->AddEntry( fractionHists[4], "K^{-}" );
  leg->AddEntry( fractionHists[3], "K^{+}" );
  leg->AddEntry( fractionHists[1], "#pi^{-}" );
  leg->AddEntry( fractionHists[0], "#pi^{+}" );
  leg->SetTextSize(0.045);
 
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(0.6*0.07);    

  partFracs->Draw("same");
  latex.DrawLatex(0.8,0.85,"Pythia 8");
  //errorFracs->Draw("same");
  
  fixOverlay();
  canv->Print("efracs.pdf"); 
  // Done.
  return 0;
}
