#ifndef AnalyzeHelper_h
#define AnalyzeHelper_h

#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>
#include <assert.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <ProfileBuilder.h>

#define JET40 (56)
#define JET80 (97)
#define JET140 (174)
#define JET200 (245)
#define JET260 (300)
#define JET320 (362)
#define JET400 (362)
#define JETMAX (4000)

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;


// A listing of different triggers. The numerical values and the amount of used
// triggers changes in time (values 0-5: trigger index, -1 : bad jet)
Int_t TriggerType( Double_t pt )
{
  if (pt < JET40 ) return -1;
  if (pt < JET80 ) return 0;
  if (pt < JET140) return 1;
  if (pt < JET200) return 2;
  if (pt < JET260) return 3;
  if (pt < JET320) return 4;
  if (pt < JETMAX) return 5;
  return -1; // For very large pt, discards the jet
}


// Tags for an event
struct EventID
{
  Long64_t run;
  Long64_t lumi;
  Long64_t event;
  EventID( Long64_t r, Long64_t l, Long64_t e): run(r), lumi(l), event(e){};

  bool operator<(const EventID& rhs) const {
    return rhs.run < run || (rhs.run == run && (rhs.lumi < lumi || 
      (rhs.lumi == lumi && rhs.event < event)));
  }

  bool operator==(const EventID& rhs) const {
    return rhs.run == run && (rhs.lumi == lumi && rhs.event == event);
  }
};


// A function to shorten obtaining the path of desired files
string GivePath( Int_t correction, Int_t isMC )
{
  string fullPath("CondFormats/JetMETObjects/data/Winter14_V1_");
  isMC ? fullPath += "MC_" : fullPath += "DATA_";
  switch( correction ){
  case 1:
    fullPath += "L1FastJet_";
    break;
  case 2:
    fullPath += "L2Relative_";
    break;
  case 3:
    fullPath += "L3Absolute_";
    break;
  case 23:
    fullPath += "L2L3Residual_";
    break;
  default:
    fullPath += "L1FastJet_";
    break;
  }

  fullPath += "AK5PF.txt";
  return fullPath;
}


FactorizedJetCorrector *EnergyCorrSetup( vector<JetCorrectorParameters> *corParams, Int_t isMC ){
  corParams->push_back( JetCorrectorParameters(GivePath(1, isMC)) );
  corParams->push_back( JetCorrectorParameters(GivePath(2, isMC)) );
  corParams->push_back( JetCorrectorParameters(GivePath(3, isMC)) );
  if (!isMC){
    corParams->push_back( JetCorrectorParameters(GivePath(23, isMC)) );
  }
  FactorizedJetCorrector *jetECor = new FactorizedJetCorrector(*corParams);
  return jetECor;
}


#endif
