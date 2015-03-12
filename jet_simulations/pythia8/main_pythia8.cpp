//////////////////////////////////////////////////////////////
// This class sorts pythia8 jets with the fastjet algorithm.// 
// See READMEi_ScriptInfo for further details.              //
// Hannu Siikonen 11.03.2015                                //
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <ctime>
// #define NDEBUG
#include <cassert>

/* Header file to access Pythia 8 program elements. */
#include "Pythia8/Pythia.h"

// ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

/* scripts */
#include "generic/help_functions.h"
#include "pythia8/pythia8_functions.h"
#include "events/PrtclEvent.h"

using namespace Pythia8;

int main(int argc, char **argv)
{
   TApplication theApp("event_generation", &argc, argv);
   size_t nEvent = 400;
   if (argc > 1) nEvent = atoi(argv[1]); assert( nEvent > 0 );
   
   /* Create Pythia instance and set it up to generate hard QCD processes
    * above pTHat = 20 GeV for pp collisions at 14 TeV. */
   Pythia pythia;
   Event& event = pythia.event;
   pythia.readFile("pythia8/pythiaSettings.cmnd");
   pythia.init();
   pythia.settings.listChanged();
  
   /* Create file on which a particle data tree is saved 
    * (before sampling to jets) */
   TFile *outFile = new TFile("pythia8_particles.root", "RECREATE");
   outFile->SetCompressionLevel(1); /* File is compressed */
   
   TTree *tree = new TTree("Pythia8Tree","Tree filled with pythia8 data.");
   PrtclEvent *pEvent = new PrtclEvent();
   
   /* Autosave when 1 Gbyte written */ 
   tree->SetAutoSave(1000000000);
   /* Set a 10 MBytes cache */
   tree->SetCacheSize(10000000);
   
   TTree::SetBranchStyle(1); /* New branch style */
   TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
   branch->SetAutoDelete(kFALSE);
   tree->BranchRef();

   /* Simulation loop: */
   Timer timer;
   timer.setParams(nEvent,100);
   timer.startTiming();  
   for (size_t ev = 0; ev != nEvent; ++ev) {
      if (!pythia.next()) continue;
      if (ev!=0&&ev%100==0) timer.printTime();
      
      for (size_t prt = 0; prt!=event.size(); ++prt){
         double status = abs( event[prt].status() );
         int tmpId = event[prt].id();
         
         /* Stable particles */
         if (event[prt].isFinal() && event[prt].isVisible()) {  
            int pi0Gamma = 0; 
            if ((tmpId == 22) && gammaChecker( event, prt )) pi0Gamma = 1;
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(),tmpId, event[prt].charge(), pi0Gamma ? 10 : 1);
         } 
         
         /* Ghost partons/particles. There might be a small overlap with stable hadrons. */
         int ghostStatus = 0;
         if ( status == 71 || status == 72 ) { 
            ghostStatus = 11; /* Partons */
         } else if ( (abs(tmpId) >= 100) && !isExcitedHadronState(event,prt,tmpId)) { 
            ghostStatus = 12; /* Hadrons */
         }
         
         if (ghostStatus) {
            pEvent->AddPrtcl(event[prt].px(),event[prt].py(),event[prt].pz(),
               event[prt].e(), tmpId, event[prt].charge(), ghostStatus);
         }
      } 
      tree->Fill();
      pEvent->Clear();
   }

   outFile = tree->GetCurrentFile();
   tree->AutoSave("Overwrite");

   delete pEvent;  pEvent = 0;
   outFile->Close();
   return 0;
}