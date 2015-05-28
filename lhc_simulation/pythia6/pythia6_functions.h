#ifndef PYTHIA6_FUNCTIONS_H
#define PYTHIA6_FUNCTIONS_H

#include <cassert>

#include "TROOT.h"
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
#include "TBranch.h"

#include <cstdlib>
#include <string>
using namespace std;

#include "../events/PrtclEvent.h"
#include "../generic/help_functions.h"


// nEvents is how many events we want.
int makeEventSample(Int_t nEvent, string fileName, Int_t nameId)
{
    // Create an instance of the Pythia event generator.
    TPythia6* pythia = new TPythia6;
    // Initialise it to run p+p at sqrt(200) GeV in CMS.
    pythia->Initialize("cms", "p", "p", 200);
    pythia->SetMRPY(1,10000*nameId);

    // Output file
    TFile* file = TFile::Open(fileName.c_str(), "RECREATE");
    if (!file || !file->IsOpen()) {
        Error("makeEventSample", "Couldn;t open file %s", fileName.c_str());
        return 1;
    }

    TTree* tree = new TTree("Pythia6Tree", "Tree filled with pythia6 data.");
    PrtclEvent *pEvent = new PrtclEvent;
    
    /* Autosave after 0.1 GByte */
    tree->SetAutoSave(100000000);
    /* Set a 10 MBytes cache */
    tree->SetCacheSize(10000000);
    
    TTree::SetBranchStyle(1); /* New branch style */
    TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
    branch->SetAutoDelete(kFALSE);
    tree->BranchRef();

    /* Event loop */
    Timer timer;
    timer.setParams(nEvent,100);
    timer.startTiming();
    for (Int_t ev = 0; ev != nEvent; ++ev) {
        pythia->GenerateEvent();
        if (ev!=0&&ev%100==0) timer.printTime();

        for (Int_t j = 0; j != pythia->GetNumberOfParticles(); ++j) {
            if (pythia->GetK(j,2) == 1) {
                pEvent->AddPrtcl(pythia->GetP(j,1),pythia->GetP(j,2),pythia->GetP(j,3),
                pythia->GetP(j,4),pythia->GetK(j,1),pythia->Pychge(pythia->GetK(j,1))/3.0,1);
            }
        }
        pEvent->Clear();
        tree->Fill();
    }

    tree->Print();

    file->Write();
    file->Close();

    return 0;
}



# endif