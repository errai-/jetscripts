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
#include "../include/MinimalEvent.h"
#include "../include/JetStore.h"
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;

string IntToString (int a)
{
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}

int main(int argc, char* argv[]) {

  string treePath = "Pythia8Tree";
  string name ="particle_storage.root";
  TChain *forest = new TChain(treePath.c_str());
  //stringstream converter;
    
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
  
  MinimalEvent *event;
  //vector<double> *px;
  forest->SetBranchAddress("event",&event);
  //forest->SetBranchAddress("px",px);
  Long64_t nentries = forest->GetEntries();
//   JetStore storaging;

  for (Long64_t jentry=0; jentry!=nentries;++jentry) {
    forest->Show(jentry); 
    forest->LoadTree(jentry);
    Int_t ientry = forest->GetEntry(jentry);
    if (ientry < 0 || !forest) break;
    size_t i = 0;
    i = 1;
//     cout << event->id[i] << " " << event->px[i] << " " << event->px.size() << endl;
    for (size_t i = 0; i != 2; ++i){ 
//       cout << event->px.size() << " " << event->px[i] << endl;
      double d = event->px[i];
      cout << d << endl;
      // " " << event->py[i] << " " << event->pz[i] << " " << event->e[i] << endl;
    }
    //storaging.Insert(1);//,event->py,event->pz,event->e);//,event->status,
      //event->id, event->mother1, event->mother2, event->daughter1, event->daughter2);
  }
  delete forest;
}