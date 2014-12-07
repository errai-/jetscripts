// File: hist.cc
// This is a simple test program.
// It studies the charged piPlusiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.
// Copyright (C) 2014 Torbjorn Sjostrand

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

// FastJet interface
#include "Pythia8/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace Pythia8;

// From CMSSW
class PtHatReweightUserHook : public UserHooks
{
  public:
  PtHatReweightUserHook(double _pt = 15, double _power = 4.5) :
  pt(_pt), power(_power) {}
  ~PtHatReweightUserHook() {}

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

//class Event { //TObject is not required by this example
class RootEvent {
public:
  // The four momentum of a jet
  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t e;
  // The direction of a jet
  Double_t eta;
  Double_t phi;
  
  Double_t chf;
  Double_t nhf;
  Double_t phf;
  Double_t elf;
  Double_t muf;

  Long64_t eventNo;
  Int_t jetsInEvt;

  RootEvent(){ }
  ~RootEvent(){ }

  void SetFourMom( double _px, double _py, double _pz, double _e ){
    px = _px; py = _py; pz = _pz; e = _e;
  }

  void SetDirections( double _eta, double _phi ){
    eta = _eta; phi = _phi;
  }

  void SetFractions( double c, double n, double p, double e, double m ){
    chf = c; nhf = n; phf = p; elf = e; muf = m;
  }

};

void EnergyIdentifier( int id, double tmpEt, double &chf, double &nhf, double &phf,
  double &elf, double &muf ){
  if ( abs( id ) == 211 ){ chf += tmpEt; // pi+-
  }else if ( abs( id ) == 2212 ){ chf += tmpEt; // p+-
  }else if ( abs( id ) == 321 ){ chf += tmpEt; // K+-
  }else if ( abs( id ) == 22 ){ phf += tmpEt; // gamma
  }else if ( abs( id ) == 11 ){ elf += tmpEt; // e
  }else if ( abs( id ) == 13 ){ muf += tmpEt; // mu
  }else if ( abs( id ) == 2112){ nhf += tmpEt; // n0
  }else if ( abs( id ) == 310 ){ nhf += tmpEt; // ks0
  }else if ( abs( id ) == 130 ){ nhf += tmpEt; // kl0
  }else if ( abs( id ) == 3122 ){ nhf += tmpEt; // lambda0
  }else if ( abs( id ) == 3112 ){ chf += tmpEt; // Sigma-
  }else if ( abs( id ) == 3222 ){ chf += tmpEt; // Sigma+ 
  }else if ( abs( id ) == 3212 ){ nhf += tmpEt; // Sigma0
  }else if ( abs( id ) == 3312 ){ chf += tmpEt; // Xi-
  }else if ( abs( id ) == 3322 ){ nhf += tmpEt; // Xi0
  }else if ( abs( id ) == 3334 ){ chf += tmpEt; // Omega-
  }else{ nhf += tmpEt; // For now other particles are interpreted to nhf 
  }
  return;
}

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("event_generation", &argc, argv);
  int weightedPt = 0;

  // Settings
  int  nEvent = 100000;

  double etaMax = 3.0;
  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  // Reweighting for event generation
  PtHatReweightUserHook* ptGenReweight;

  if (weightedPt){
    ptGenReweight = new PtHatReweightUserHook();
    pythia.setUserHooksPtr( ptGenReweight );
  }

  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");

  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init( 2212, 2212, 14000.);
  pythia.settings.listChanged();

  // Create file on which histogram(s) can be saved.
  TFile outFile("pythiajets.root", "RECREATE");

  // Fastjet analysis - select algorithm and parameters
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

  // Creating a tree
  TTree pythiaTree("pythia8","pythia 8 fastjet simulation data");
  RootEvent rEvent;

  // Create a branch
  pythiaTree.Branch("px", &rEvent.px);
  pythiaTree.Branch("py", &rEvent.py);
  pythiaTree.Branch("pz", &rEvent.pz);
  pythiaTree.Branch("e", &rEvent.e);
  pythiaTree.Branch("eta", &rEvent.eta);
  pythiaTree.Branch("phi", &rEvent.phi);
  pythiaTree.Branch("jetsinevent", &rEvent.jetsInEvt);
  pythiaTree.Branch("chf", &rEvent.chf);
  pythiaTree.Branch("nhf", &rEvent.nhf);
  pythiaTree.Branch("phf", &rEvent.phf);
  pythiaTree.Branch("elf", &rEvent.elf);
  pythiaTree.Branch("muf", &rEvent.muf);
  pythiaTree.Branch("events", &rEvent.eventNo);

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Reset Fastjet input
    fjInputs.resize(0);
    // Particle loop
    for (int i = 0; i != event.size(); ++i) if (event[i].isFinal()) {
      if ( !event[i].isVisible() ) continue;
      fastjet::PseudoJet particleTemp = event[i];
      particleTemp.set_user_index( i );
      fjInputs.push_back( particleTemp );
    }

    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveJets = clustSeq.inclusive_jets( pTMin );
    sortedJets    = sorted_by_pt(inclusiveJets);
    
    rEvent.jetsInEvt = sortedJets.size();
    rEvent.eventNo = nEvent;
    
    for (unsigned int i = 0; i != sortedJets.size(); ++i ){
      // Only count jets that have |eta| < etaMax
      //if (abs(sortedJets[i].pseudorapidity()) > etaMax ) continue;

      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      double chf = 0, nhf = 0, phf = 0, elf = 0, muf = 0, etSum = 0;

      for (unsigned int idx = 0; idx != jetParts.size(); ++idx ){
        double tmpEt = event[ jetParts[idx].user_index() ].e();
        etSum += tmpEt;
        EnergyIdentifier( event[ jetParts[idx].user_index() ].id(), tmpEt, chf,
          nhf, phf, elf, muf );
      }

      chf /= etSum; nhf /= etSum; phf /= etSum; elf /= etSum; muf /= etSum;
      rEvent.SetFourMom( sortedJets[i].px(), sortedJets[i].py(),
        sortedJets[i].pz(), sortedJets[i].e() );
      rEvent.SetDirections( sortedJets[i].eta(), sortedJets[i].phi() );
      rEvent.SetFractions( chf, nhf, phf, elf, muf );
      pythiaTree.Fill();
    }
  }

  //gPad->WaitPrimitive();

  pythiaTree.Write();

  if (weightedPt){
    delete ptGenReweight;
  }
  // Done.
  return 0;
}

