#ifndef JETSTORE_H
#define JETSTORE_H


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
#include "../JetSorter/jetsorter_auxiliary.h"
#include "MinimalEvent.h"
#include <boost/lexical_cast.hpp>

using namespace Pythia8;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;


class JetStore{
private:
  vector<fastjet::PseudoJet> Momentum;
  vector<int> status;
  vector<int> id;
  vector<int> mother1;
  vector<int> mother2;
  vector<int> daughter1;
  vector<int> daughter2;
public:
  JetStore(){}
  ~JetStore(){}
  
  void Insert(double px, double py=0, double pz=0, double e=0, int _status=0, int _id=0,
    int _mother1=0, int _mother2=0, int _daughter1=0, int _daughter2=0){
    fastjet::PseudoJet tmpJet(px,py,pz,e);
    Momentum.push_back(tmpJet);
    status.push_back(_status);
    id.push_back(_id);
    mother1.push_back(_mother1);
    mother2.push_back(_mother2);
    daughter1.push_back(_daughter1);
    daughter2.push_back(_daughter2);
  }
  
  void Clear(){
    Momentum.clear();
    status.clear();
    id.clear();
    mother1.clear();
    mother2.clear();
    daughter1.clear();
    daughter2.clear();
  }
};











#endif JETSTORE_H