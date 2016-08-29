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
// fastjet
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// tdrStyle
//#include "tdrstyle_mod1.C"
// scripts
#include "../../JetSorter/jetsorter_auxiliary.h"
#include "MinimalEvent.h"
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;


class JetStore{
private:
  vector<fastjet::PseudoJet> Momentum;
  vector<int> status;
  vector<int> id;
public:
  JetStore(){}
  ~JetStore(){}
  
  void Insert(double,double=0,double=0,double=0,int=0,int=0);
  
  void Clear();
  

};











#endif JETSTORE_H