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
#include "../include/RootJetSort.h"

using std::cout;
using std::endl;
using std::string;

string IntToString (int a)
{
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}

int main(int argc, char* argv[]) {

  // Create a chain that includes all the necessary trees
  string treePath = "Pythia8Tree";
  string name ="particle_storage.root";
  TChain *forest = new TChain(treePath.c_str());
    
  // Check that all the subtrees are included in the chain
  TFile *probe = new TFile("particle_storage.root");
  TIter *iter = new TIter(probe->GetListOfKeys());
  TKey *tmpKey = 0;
  string alterName;
  while ( tmpKey = (TKey *) iter->Next() ){
    if ( strcmp( tmpKey->GetName(), treePath.c_str() ) !=  0) continue;
    alterName = treePath;
    alterName += ";";
    alterName += IntToString(tmpKey->GetCycle());
    forest->AddFile(name.c_str(),forest->kBigNumber,alterName.c_str());
  }
  delete iter;
  delete probe;
  
  RootJetSort treeHandle(forest);
  
  treeHandle.EventLoop();
  
  treeHandle.WriteResults();
  
  delete forest;
}