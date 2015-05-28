// plots jet energy fractions, see README_ScriptInfo

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

#include <assert.h>

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

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod4.C"

using std::string;
using std::vector;
using std::cout;
using std::endl;

void drawEfracsPtcut(std::string fileName) {

  // Create file on which histogram(s) can be saved.
  TFile *inFile = new TFile(fileName.c_str(), "READ");
  THStack* partFracs = new THStack("particlestacks", "");
  THStack* errorFracs = new THStack("errorstacks", "" );
  
  // Read histograms.
  vector<TProfile*> fractionProfiles;
  fractionProfiles.push_back( (TProfile*) inFile->Get( "ch bins" ));
  fractionProfiles.push_back( (TProfile*) inFile->Get( "ph bins" ));
  fractionProfiles.push_back( (TProfile*) inFile->Get( "nh bins" ));
  fractionProfiles.push_back( (TProfile*) inFile->Get( "e bins" ));
  fractionProfiles.push_back( (TProfile*) inFile->Get( "mu bins" ));

  // Create assosiated histograms for the profiles
  vector<TH1D*> fractionHists;
  for (int idx = 0; idx != 5; ++idx){
    fractionHists.push_back( fractionProfiles[idx]->ProjectionX("","") );
  }

  // Set fill colors
  fractionHists[0]->SetFillColor(kRed-7);
  fractionHists[1]->SetFillColor(kBlue-7);
  fractionHists[2]->SetFillColor(kGreen-6);
  fractionHists[3]->SetFillColor(kCyan-6);
  fractionHists[4]->SetFillColor(kMagenta-6);

  for (unsigned int i = 0; i != fractionHists.size(); ++i){
    partFracs->Add( fractionHists[i] );
    errorFracs->Add( fractionProfiles[i]->ProjectionX() );
  }

  TCanvas *canv = new TCanvas("c1","c1",600,600);
  canv->SetSelected(canv);
  canv->cd();
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.35);
  leg->AddEntry( fractionHists[4], "Muons" );
  leg->AddEntry( fractionHists[3], "Electrons" );
  leg->AddEntry( fractionHists[2], "Neutral hadrons" );
  leg->AddEntry( fractionHists[1], "Photons" );
  leg->AddEntry( fractionHists[0], "Charged hadrons" );
  leg->SetTextSize(0.03);
  // Show histogram  
  setTDRStyle();
  canv->UseCurrentStyle();
  stackModify(fractionHists[0]);
 
  partFracs->SetHistogram( fractionHists[0] );
  errorFracs->SetHistogram( fractionHists[0] );
  partFracs->Draw();
  //errorFracs->Draw("same");
  pythiaFinal();
  leg->Draw();
  canv->Modified();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();

  // Save histogram on file and close file.
  canv->Write();
}
