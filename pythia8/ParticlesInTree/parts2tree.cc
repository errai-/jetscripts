// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
// Nice libraries from C
#include <cmath>
#include <ctime>

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
#include "../JetSorter/jetsorter_auxiliary.h"
#include "MinimalEvent.h"

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
  int ptBins = 48.;
  const double ptRange[]=
 //{1, 5, 6, 8, 10, 12, 15
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};//,
    //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
    //3637, 3832, 4037};//

  int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 20.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range

  // Create Pythia instance and set it up to generate hard QCD processes
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

  // Create file on which a particle data tree is saved (before sampling to jets)
  TFile outFile("particle_storage.root", "RECREATE");
  // An object to store the events
  MinimalEvent *mEvent = new MinimalEvent;
  // A tree is created
  TTree *eventStorage = new TTree("Pythia8Tree","event storage");
  eventStorage->Branch("event","MinimalEvent",&mEvent);
  
  std::clock_t start = std::clock();
  double time_processor = 0; int hours; int minutes; int seconds; 

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // event.list();
    if (iEvent!=0&&iEvent%100==0){
      time_processor = (std::clock() - start)/(( (double) CLOCKS_PER_SEC ) );
      time_processor = time_processor*( ((double) nEvent)/iEvent-1); 
      minutes =  time_processor/60; hours = minutes/60;
      seconds = time_processor-60*minutes;
      minutes = minutes - hours*60;
      cout << iEvent << " events created, ETA : " << hours << "h" <<
        minutes << "m" << seconds << "s." << endl;
    }  
    // Particle loop
    for (int i = 0; i != event.size(); ++i) {
      double status = abs( event[i].status() );
      if (event[i].isFinal() && event[i].isVisible()) {  
        mEvent->SetVals(event[i].px(),event[i].py(),event[i].pz(),event[i].e(), 
          event[i].status(), event[i].id(), event[i].mother1(), 
          event[i].mother2(), event[i].daughter1(), event[i].daughter2() );
      }
    }
    eventStorage->Fill();
    mEvent->Nullify();
  }

  eventStorage->Print();
  eventStorage->Write();
  
  delete eventStorage;
  
  // Done.
  return 0;
}
