#ifndef PYTHIA6_FUNCTIONS_H
#define PYTHIA6_FUNCTIONS_H

#include <cassert>

#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <cstdlib>
using namespace std;

#define FILENAME   "pythia.root"
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define HISTNAME   "ptSpectra"
#define PDGNUMBER  211

// nEvents is how many events we want.
int makeEventSample(Int_t nEvents)
{
   // Create an instance of the Pythia event generator ...
   TPythia6* pythia = new TPythia6;

   // ... and initialise it to run p+p at sqrt(200) GeV in CMS
   pythia->Initialize("cms", "p", "p", 200);

   // Open an output file
   TFile* file = TFile::Open(FILENAME, "RECREATE");
   if (!file || !file->IsOpen()) {
      Error("makeEventSample", "Couldn;t open file %s", FILENAME);
      return 1;
   }

   TTree* tree = new TTree(TREENAME, "Pythia 6 tree");

   // It's a TClonesArray of TMCParticle objects
   TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
   tree->Branch(BRANCHNAME, &particles);

   // Now we make some events
   for (Int_t i = 0; i < nEvents; i++) {
         
      if (i % 100 == 0) {
         cout << "Event # " << i << endl;
      }

      pythia->GenerateEvent();

      for (Int_t j = 0; j < pythia->GetNumberOfParticles(); ++j) {
         cout << pythia->GetK(j,1) << " " << pythia->GetK(j,2) << endl;
      }
      tree->Fill();
   }

   tree->Print();

   // and now we flush and close the file
   file->Write();
   file->Close();

   return 0;
}



# endif