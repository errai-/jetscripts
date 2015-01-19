// This holds all the auxiliary functions and classes of jetsorter
#ifndef JETSORTER_AUXILIARY
#define JETSORTER_AUXILIARY

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

#include <math.h>
// FastJet interface
#include "Pythia8Plugins/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod1.C"

using namespace Pythia8;

// From CMSSW
class PtHatReweightUserHook : public UserHooks
{ 
  public:
  PtHatReweightUserHook(double _pt = 15, double _power = 4.5) :
  pt(_pt), power(_power) {}
  virtual ~PtHatReweightUserHook() {}
  
  virtual bool canBiasSelection() { return true; }
  
  virtual double biasSelectionBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent)
  { 
    //the variable selBias of the base class should be used;
    if ((sigmaProcessPtr->nFinal() == 2)) {
    selBias = pow(phaseSpacePtr->pTHat() / pt, power);
    return selBias;
    }
    selBias = 1.;
    return selBias;
  }
  
  private:
  double pt, power;
};

// A function that checks whether a photon is originated from a pi0 and that
// the energy of the photon-pair corresponds to the pion. returns 0 if
// the origin is not a pion with good energy and 1 if it is
int gammaChecker( Event &event, int idx ){
  int mother = event[idx].mother1();
  if ( event[mother].id() != 111 ) return 0;
  double eDifference = abs( event[mother].e() -
    event[event[mother].daughter1()].e() - event[event[mother].daughter2()].e() );
  if ( eDifference < 0.001 ) return 1;
  return 0;
}

double deltaR( double phi1, double phi2, double eta1, double eta2 ){
  double dPhi = phi1 - phi2;
  double dEta = eta1 - eta2;
  return pow( pow( dPhi, 2 ) + pow( dEta, 2 ) , 0.5 );
}

void histFiller( vector<TProfile*> &hists, double pt, double eTot, double piPlus,
  double piMinus, double pi0Gamma, double kaPlus, double kaMinus, double kSZero,
  double kLZero, double proton, double aproton, double neutron, double aneutron,
  double gamma, double lambda0, double sigma, double elecmuon, double others ){
  hists[0]->Fill( pt, piPlus/eTot ); hists[1]->Fill( pt, piMinus/eTot );
  hists[2]->Fill( pt, pi0Gamma/eTot ); hists[3]->Fill( pt, kaPlus/eTot );
  hists[4]->Fill( pt, kaMinus/eTot ); hists[5]->Fill( pt, kSZero/eTot );
  hists[6]->Fill( pt, kLZero/eTot ); hists[7]->Fill( pt, proton/eTot );
  hists[8]->Fill( pt, aproton/eTot ); hists[9]->Fill( pt, neutron/eTot );
  hists[10]->Fill( pt, aneutron/eTot ); hists[11]->Fill( pt, gamma/eTot );
  hists[12]->Fill( pt, lambda0/eTot ); hists[13]->Fill( pt, sigma/eTot );
  hists[14]->Fill( pt, elecmuon/eTot ); hists[15]->Fill( pt, others/eTot );
}

int isBottom( int id ) {
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

int isCharm( int id ) {
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

int isStrange( int id ) {
  int code1;
  int code2;
  bool tmpHasStrange = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 3 || code2 == 3) tmpHasStrange = true;
  return tmpHasStrange;
}

int isDown( int id ) {
  int code1;
  int code2;
  bool tmpHasDown = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 2 || code2 == 2) tmpHasDown = true;
  return tmpHasDown;
}

int isUp( int id ) {
  int code1;
  int code2;
  bool tmpHasUp = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 1 || code2 == 1) tmpHasUp = true;
  return tmpHasUp;
}

int statusCheck( int id1, int id2 ){
  if ( id1 == 5 && isBottom( id2 ) ) return 1;
  if ( id1 == 4 && isCharm( id2 ) ) return 1;
  if ( id1 == 3 && isStrange( id2 ) ) return 1;
  if ( id1 == 2 && isDown( id2 ) ) return 1;
  if ( id1 == 1 && isUp( id2 ) ) return 1;
  return 0;
}

int isExcitedState( Event &event, int idx, int id ) {
  int d1 = event[idx].daughter1(), d2 = event[idx].daughter2();
  if (d2!=0){
    if (d1 < d2){
      for (int i = d1; i <= d2; i++){
        if ( statusCheck( id, event[i].id() ) ) return 1;
      }
    } else {
      if ( statusCheck( id, event[d1].id() ) ) return 1;
      if ( statusCheck( id, event[d2].id() ) ) return 1;
    }
  } else if (d1!=0){
    if ( statusCheck( id, event[d1].id() ) ) return 1;
  }
  return 0;
}

int ChargeSign( int id ){
  if ( id == 1 ) return 1;
  if ( id == -2 ) return 1;
  if ( id == -3 ) return 1;
  if ( id == 4 ) return 1;
  if ( id == -5 ) return 1;
  if ( id == 6 ) return 1;
  if ( id == -1 ) return -1;
  if ( id == 2 ) return -1;
  if ( id == 3 ) return -1;
  if ( id == -4 ) return -1;
  if ( id == 5 ) return -1;
  if ( id == -6 ) return -1;
  return 1;
}


#endif
