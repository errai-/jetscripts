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
#include "RootJetSort.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;



int main(int argc, char* argv[]) {

  // Create a chain that includes the necessary tree
  // There are a couple of different data options available.
  string treePath = "Pythia8Tree";
  string name ="pythia8_particles.root";
  if (argc > 1) {
    if ( argv[2] == "1" ){
      treePath = "HerwigTree";
      name = "herwig_particles.root";
    }
  }
  
  TChain *forest = new TChain(treePath.c_str());
  
  // This opens the tree with the highest key value with the given treePath
  forest->AddFile(name.c_str());
  
  RootJetSort treeHandle(forest);
  
  treeHandle.EventLoop();
  
  treeHandle.WriteResults();
  
  delete forest;
}

