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

// ROOT
#include "Riostream.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
#include "help_functions.h"
#include "pythia8_functions.h"
#include "SimEvent.h"

using namespace Pythia8;
using std::cout;
using std::endl;
using std::vector;

int main(int argc, char **argv)
{
  // Settings
  TApplication theApp("event_generation", &argc, argv);
  size_t nevent = 400;     // by default create 400 events
  Int_t comp   = 1;       // by default file is compressed
  Int_t arg5   = 600;     //default number of tracks per event
  if (argc > 1)  nevent = atoi(argv[1]);
  Int_t branchStyle = 1; //new style by default 
  bool saveEverything = false;
  if (argc > 2 && atoi(argv[2]) != 0){ saveEverything = true; }
  
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;

  pythia.readString("HardQCD:all = on");
  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PartonLevel:ISR=on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.5");
  pythia.readString("PhaseSpace:bias2SelectionRef = 15");
  pythia.readString("Beams:eCM = 7000.");
  pythia.init();
  pythia.settings.listChanged();
  //pythia.particleData.listAll();
   
  // Create file on which a particle data tree is saved (before sampling to jets)
  TFile *outFile = new TFile("pythia8_particles.root", "RECREATE");
  TTree *tree = new TTree("Pythia8Tree","Tree filled with pythia8 data.");
  SimEvent *sEvent = new SimEvent();
   
  Int_t bufsize = 64000/4;
  outFile->SetCompressionLevel(comp);
   
  // Create a ROOT Tree and one superbranch
  tree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
  tree->SetCacheSize(10000000);  // set a 10 MBytes cache (useless when writing local files)
  TTree::SetBranchStyle(branchStyle);
  TBranch *branch = tree->Branch("event", &sEvent, bufsize,2);
  branch->SetAutoDelete(kFALSE);
  if(branchStyle) tree->BranchRef();
  Float_t ptmin = 1;

  Timer timer;
  timer.set_params(nevent,100);
  timer.start_timing();  
  for (size_t ev = 0; ev < nevent; ev++) {
    if (!pythia.next()) continue;
    // event.list();
    if (ev!=0&&ev%100==0){
      timer.print_time();
    }
    
    for (size_t prt = 0; prt!=event.size(); ++prt){
      double status = abs( event[prt].status() );
      int tmpId = event[prt].id();
      if (saveEverything || (event[prt].isFinal() && event[prt].isVisible()) ) {  
	int pi0Gamma = 0;
	if ((tmpId == 22) && gammaChecker( event, prt )) pi0Gamma = 1; // Indicate pi0 photons
	sEvent->Build(event[prt].px(),event[prt].py(),event[prt].pz(),event[prt].e(), tmpId, event[prt].charge(), pi0Gamma);
      } else if ( status == 71 || status == 72 || status == 61 || status == 62 || status == 63 ){
	int isExcState = ((isExcitedState(event,prt,tmpId)) ? 1 : 0);
        sEvent->Build(event[prt].px(),event[prt].py(),event[prt].pz(),event[prt].e(), tmpId, event[prt].charge(), 0, 1, isExcState); 
      } 

    } 
    
    tree->Fill();  //fill the tree
    sEvent->Clear();
  }
  
  outFile = tree->GetCurrentFile(); //just in case we switched to a new file
  tree->AutoSave("Overwrite");
  //tree->Print();
  // We own the event (since we set the branch address explicitly), we need to delete it.
  delete sEvent;  sEvent = 0;
   
  outFile->Close();
  return 0;
} 