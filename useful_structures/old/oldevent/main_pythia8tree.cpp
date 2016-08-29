// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
// Nice libraries from C
#include <cmath>
#include <ctime>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

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

// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
#include "jetsorter_auxiliary.h"
#include "SimEvent.h"
#include "help_functions.h"

using namespace Pythia8;
using std::cout;
using std::endl;
using std::vector;

int main(int argc, char* argv[]) {

  // Settings
  TApplication theApp("event_generation", &argc, argv);
  int  nEvent = 100;
  if (argc > 1){ nEvent = atoi(argv[1]); }
  int weightedPt = 1;
  bool everything = false;
  if (argc > 2 && atoi(argv[2]) != 0){ everything = true; }

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  // Reweighting for event generation
  PtHatReweightUserHook ptGenReweight;

  if (weightedPt){ pythia.setUserHooksPtr( &ptGenReweight ); }
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PartonLevel:ISR=on");
  pythia.readString("Beams:eCM = 7000.");
  pythia.init();
  pythia.settings.listChanged();
  //pythia.particleData.listAll();
  
  // Create file on which a particle data tree is saved (before sampling to jets)
  TFile outFile("particle_storage.root", "RECREATE");
  // An object to store the events
  SimEvent *sEvent;
  // A tree is created
  TTree *eventStorage = new TTree("Pythia8Tree","event storage");
  eventStorage->Branch("event","SimEvent",&sEvent,16000,2);

  Timer timer;
  timer.set_params(nEvent,100);
  timer.start_timing();  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // event.list();
    if (iEvent!=0&&iEvent%100==0){
      timer.print_time();
    }
    // Particle loop
    for (int i = 0; i != event.size(); ++i) {
      double status = abs( event[i].status() );
      if (everything || (event[i].isFinal() && event[i].isVisible()) ) {  
	      int tmpId = event[i].id();
	      tmpId -= ( ((tmpId == 22) && gammaChecker( event, i )) ? 2 : 0 ); // Indicate pi0 photons
        sEvent->Build(event[i].px(),event[i].py(),event[i].pz(),event[i].e(), tmpId);
      }
      //cout << endl;
    }
    eventStorage->Fill();
    sEvent->Clear();
  }

  eventStorage->Print();
  eventStorage->AutoSave("Overwrite"); // Using Write makes the data display in two trees.
  
  outFile.Close();
  
  // Done.
  return 0;
}
