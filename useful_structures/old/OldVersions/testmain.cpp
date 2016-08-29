#include <iostream>
#include <vector>
#include <string>
#include <sstream>
// Nice libraries from C
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"
// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
// ROOT, for saving a file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
// ROOT, Trees
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TKey.h"
#include "TList.h"
#include "TIterator.h"
#include "TString.h"

// tdrStyle
//#include "tdrstyle_mod1.C"
// scripts
#include "../../JetSorter/jetsorter_auxiliary.h"
#include "../include/Hope.h"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[]) {

  string treePath = "Pythia8Tree";
  string name ="particle_storage.root";
  TChain *forest = new TChain(treePath.c_str());

  Hope suu;
  suu.Loop();
}