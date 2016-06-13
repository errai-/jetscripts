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


int mergeBins() {

  // Read histograms 
  TFile *inFile = new TFile("sortedjets.root", "READ");

  vector<TProfile*> fractionProfiles;
  for (int idx = 0; idx != 16; ++idx){
    std::stringstream tmpString("");
    tmpString << "a" << idx; 
    // a = all, g = gluonjets, q = quarkjets, lq = low pt quarks, hq = high pt quarks
    fractionProfiles.push_back( (TProfile*) inFile->Get( tmpString.str().c_str()) );
  }
  cout << fractionProfiles.size();

  const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 967,
    1101, 1248, 1410, 1684, 2000};//,

  vector<TProfile*> fractionHists;
  TFile *outFile = new TFile("sortednets.root", "RECREATE");
  for (int idx = 0; idx != 16; ++idx){
    std::stringstream tmpString("");
    tmpString << "a" << idx; 
    fractionHists.push_back( (TProfile*) fractionProfiles[idx]->Rebin(40,tmpString.str().c_str(),ptRange));
    fractionHists[idx]->Write();
  }

  return 0;
}
