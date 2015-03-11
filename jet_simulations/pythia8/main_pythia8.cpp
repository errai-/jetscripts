//////////////////////////////////////////////////////////////
// This class sorts pythia8 jets with the fastjet algorithm.// 
// See READMEi_ScriptInfo for further details.              //
// Hannu Siikonen 11.03.2015                                //
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <ctime>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
#include "include/help_functions.h"
#include "pythia8/pythia8_functions.h"
#include "events/PrtclEvent.h"

using namespace Pythia8;

int main(int argc, char **argv)
{
   TApplication theApp("event_generation", &argc, argv);
   size_t nEvent = 400;    // by default create 400 events
   if (argc > 1)  nEvent = atoi(argv[1]);

   /* Create Pythia instance and set it up to generate hard QCD processes
    * above pTHat = 20 GeV for pp collisions at 14 TeV. */
   Pythia pythia;
   Event& event = pythia.event;
   //pythia.readFile("pythiaSettings.cmnd");
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
  
   /* Create file on which a particle data tree is saved 
    * (before sampling to jets) */
   TFile *outFile = new TFile("pythia8_particles.root", "RECREATE");
   outFile->SetCompressionLevel(1); // File is compressed
   
   TTree *tree = new TTree("Pythia8Tree","Tree filled with pythia8 data.");
   PrtclEvent *pEvent = new PrtclEvent();
   
   // autosave when 1 Gbyte written 
   tree->SetAutoSave(1000000000);
   // set a 10 MBytes cache (useless when writing local files) 
   tree->SetCacheSize(10000000);
   
   TTree::SetBranchStyle(1); // New branch style
   TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
   branch->SetAutoDelete(kFALSE);
   tree->BranchRef();

   // Simulation loop:
   Timer timer;
   timer.setParams(nEvent,100);
   timer.startTiming();  
   for (size_t ev = 0; ev != nEvent; ++ev) {
      if (!pythia.next()) continue;
      if (ev!=0&&ev%100==0) timer.printTime();
      
      for (size_t prt = 0; prt!=event.size(); ++prt){
         double status = abs( event[prt].status() );
         int tmpId = event[prt].id();
         if (event[prt].isFinal() && event[prt].isVisible() ) {  
            // Indicate pi0 photons
            int pi0Gamma = 0;
            if ((tmpId == 22) && gammaChecker( event, prt )) pi0Gamma = 1;
            // Add final particles into the tree
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(),tmpId, event[prt].charge(), pi0Gamma);
         } else if ( status == 71 || status == 72 ) {
            /* HadronAndPartonSelector.cc in cmssw contains information 
             * of interesting status codes. */
            int isExcState = ((isExcitedState(event,prt,tmpId)) ? 1 : 0);
            // Add non-permanent particles and partons for flavour indication
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(), tmpId, event[prt].charge(), 0, 1,isExcState); 
         } 
      } 
      tree->Fill();  //fill the tree
      pEvent->Clear();
   }

   // Close the processes:
   outFile = tree->GetCurrentFile();
   tree->AutoSave("Overwrite");

   delete pEvent;  pEvent = 0;
   
   outFile->Close();
   return 0;
}