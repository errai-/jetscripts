// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
// Nice libraries from C
#include <cmath>
#include <ctime>
#include <cstdint>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
// FastJet interface
#include "Pythia8Plugins/FastJet3.h"

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
#include "../../JetSorter/jetsorter_auxiliary.h"
#include "../include/MinimalEvent.h"
#include "../include/help_functions.h"

using namespace Pythia8;
using std::cout;
using std::endl;
using std::vector;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("event_generation", &argc, argv);
  int weightedPt = 1;

  // Settings
  int  nEvent = 100;
  // 0 for gluon jet, 1 for all quarks, 2 for light quarks, 3 for heavy quarks
  if (argc > 1){
    nEvent = atoi(argv[1]);
  }

  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range

  //// Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  // Reweighting for event generation
  PtHatReweightUserHook ptGenReweight;

  if (weightedPt){
    pythia.setUserHooksPtr( &ptGenReweight );
  }
  
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PartonLevel:ISR=on");
  //pythia.particleData.listAll();

  pythia.readString("Beams:eCM = 14000.");
  pythia.init();
  pythia.settings.listChanged();

  // An object to store the events
  MinimalEvent tmpEvent;
  // A timer to count ETA's
  Timer timer(nEvent,100);
  // File to write into
  std::ofstream output;
  output.open("pythia8data.txt", std::ofstream::out | std::ofstream::app);
  output << std::fixed << std::setprecision(60);
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
      if (event[i].isFinal() && event[i].isVisible()) {
	int tmpId = event[i].id();
	tmpId -= ( ((tmpId == 22) && gammaChecker( event, i )) ? 2 : 0 ); // Indicate pi0 photons
        tmpEvent.SetVals(event[i].px(),event[i].py(),event[i].pz(),event[i].e(), tmpId );
      }
    }
    tmpEvent.Write(&output);
    tmpEvent.Nullify();
  }

  output.close();
  // Done.
  return 0;
}
