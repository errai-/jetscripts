#ifndef PYTHIA6_FUNCTIONS_H
#define PYTHIA6_FUNCTIONS_H

#include <cassert>

#include "TROOT.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
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


namespace
{
    bool Pythia6ParticleLoop(TPythia6*, PrtclEvent*, const int);
    int MuonFinder(TPythia6*, int);

    int Pythia6Settings(TPythia6* pythia, int mode) {
        if (mode == 1) {
            // Standard QCD
            pythia->SetMSEL(1);
            // Min and max pthat
            pythia->SetCKIN(3,30);
            pythia->SetCKIN(4,3000);
        } else if (mode == 2) {
            // photon+jets
            //pythia->SetMSEL(10);
            pythia->SetMSEL(0);
            pythia->SetMSUB(14,1);
            pythia->SetMSUB(29,1);
            pythia->SetMSUB(115,1);
            // Min and max pthat
            pythia->SetCKIN(3,10);
            pythia->SetCKIN(4,3000);
        } else if (mode == 3) {
            // Z+jets
            pythia->SetMSEL(13); // 11 would be the vanilla y*/Z
            // Leave only decay to muons on
            pythia->SetMDME( 174,1,0 ); // Z decay to d dbar
            pythia->SetMDME( 175,1,0 ); // Z decay to u ubar
            pythia->SetMDME( 176,1,0 ); // Z decay to s sbar
            pythia->SetMDME( 177,1,0 ); // Z decay to c cbar
            pythia->SetMDME( 178,1,0 ); // Z decay to b bbar
            pythia->SetMDME( 179,1,0 ); // Z decay to t tbar
            pythia->SetMDME( 182,1,0 ); // Zee
            pythia->SetMDME( 183,1,0 ); // Znuenue
            pythia->SetMDME( 184,1,1 ); // Zmumu
            pythia->SetMDME( 185,1,0 ); // Znumunumu
            pythia->SetMDME( 186,1,0 ); // Ztautau
            pythia->SetMDME( 187,1,0 ); // Znutaunutau
            // Min and max mhat
            pythia->SetCKIN(1,40);
            pythia->SetCKIN(2,-1);
            // Min and max pthat
            pythia->SetCKIN(3,15);
            pythia->SetCKIN(4,3000);
        } else {
            cout << "Select a mode! (dijet: 1, gammajet: 2, Zjet: 3)" << endl;
            return -1;
        }
        pythia->SetMSTU(21,1); // Check for errors
        pythia->SetMSTJ(22,2); // Unstable particle decay:
        pythia->SetPARJ(71,10); // ctau = 10 mm
        pythia->SetMSTP(33,0); // no K factors in hard cross sections
        pythia->SetMSTP(2,1); // which order running alphaS
        pythia->SetMSTP(51,10042); // Structure function (PDF CTEQ6L1)
        pythia->SetMSTP(52,2); // LHAPDF

        pythia->SetPARP(82,1.921); // pt cutoff, multiparton interactions
        pythia->SetPARP(89,1800.); // sqrts for which parp82 is set 
        pythia->SetPARP(90,0.227); // MPI: rescaling power

        pythia->SetMSTP(142,2); // Turn on Pt reweighting

        pythia->SetMSTP(95,6); // Color reconnection setParams
        pythia->SetPARP(77,1.016); // CR
        pythia->SetPARP(78,0.538); // CR

        pythia->SetPARP(80,0.1); // Prob. colored parton from BBR

        pythia->SetPARP(83,0.356); // MPI matter distribution
        pythia->SetPARP(84,0.651); // MPI matter distribution

        pythia->SetPARP(62,1.025); // ISR cutoff

        pythia->SetMSTP(91,1); // Gaussian primordial KT
        pythia->SetMSTP(93,10.0); // Primordial KT-max

        pythia->SetMSTP(81,21); // MPI
        pythia->SetMSTP(82,4); // MPI model

        //pythia->SetMSTP(61,0); // ISR off
        //pythia->SetMSTP(71,0); // FSR off
        //pythia->SetMSTP(81,0); // MPI off

        // pythia->SetMSTP(111,0); // Hadronization off

        pythia->Initialize("cms", "p", "p", 13000);
        //pythia->Initialize("cms", "p", "p", 8000);
        return 1;
    }

    // nEvents is how many events we want.
    int Pythia6EventLoop(Int_t nEvent, string fileName, Int_t nameId, const int mode)
    {
        /* Create an instance of the Pythia event generator: */
        TPythia6* pythia = new TPythia6;
        /* Set a seed value according to the run index and make sure it is used: */
        pythia->SetMRPY(1,10000*nameId);
        pythia->SetMRPY(2,0);
        /* Other settings: */
        if (Pythia6Settings(pythia,mode) == -1) return -1;

        /* Output file: */
        TFile* file = TFile::Open(fileName.c_str(), "RECREATE");
        if (!file || !file->IsOpen()) {
            Error("makeEventSample", "Couldn;t open file %s", fileName.c_str());
            return 1;
        }

        /* Output tree: */
        TTree* tree = new TTree("Pythia6Tree", "Tree filled with pythia6 data.");
        PrtclEvent *pEvent = new PrtclEvent;
        /* Autosave after 0.1 GByte */
        tree->SetAutoSave(100000000);
        /* Set a 10 MBytes cache */
        tree->SetCacheSize(10000000);
        /* New branch style */
        TTree::SetBranchStyle(1); 
        TBranch *branch = tree->Branch("event", &pEvent, 32000,4);
        branch->SetAutoDelete(kFALSE);
        tree->BranchRef();

        /* Simulation loop */
        int timerStep = 1000;
        Timer timer; timer.setParams(nEvent,timerStep); timer.startTiming();
        std::size_t ev = 0;
        while (ev != nEvent) { 
            pythia->GenerateEvent();

            if ( Pythia6ParticleLoop(pythia,pEvent,mode) ) { tree->Fill(); ++ev; }
            if (ev%timerStep==0) timer.printTime();
            pEvent->Clear();
        }

        pythia->Pylist(1);
        pythia->Pystat(1);

        tree->AutoSave("Overwrite");

        delete pEvent; pEvent = 0;
        file->Close();

        return 0;
    }

    bool Pythia6ParticleLoop(TPythia6* pythia, PrtclEvent* pEvent, const int mode) {
        int muon1 = 12, muon2 = 13, gamma = 9;
        pEvent->fWeight = 1./pythia->GetVINT(99);
        if (mode==3) while (abs(pythia->GetK(muon1,2))!=13) { ++muon1; ++muon2; }
        for (Int_t j = 1; j <= pythia->GetN(); ++j) {
            // j == 7,8: outgoing particles in the hardest subprocess
            if (j==7 || j==8) {
                if (mode==2 && pythia->GetK(j,2)==22) continue;
                if (mode==3 && pythia->GetK(j,2)==23) continue;

                if ( pythia->GetK(j,1) != 21 ) {
                    cout << "False functionality in hardest subprocess" << endl;
                    break;
                }
                pEvent->AddPrtcl(pythia->GetP(j,1),pythia->GetP(j,2),pythia->GetP(j,3),
                    pythia->GetP(j,4),pythia->GetK(j,2),3);
            }

            if (j==muon1 && mode==3) muon1 = MuonFinder( pythia, j );
            if (j==muon2 && mode==3) muon2 = MuonFinder( pythia, j );
            if (j==gamma && mode==2 && pythia->GetK(j,1)>10) gamma = pythia->GetK(j,4);

            // Stable particles
            if (pythia->GetK(j,1) <= 10) {
                if ( ((j==muon1 || j==muon2)&&mode==3) || (j==gamma&&mode==2) ) {
                    // Z+jets, Photon+jets
                    pEvent->AddPrtcl(pythia->GetP(j,1),pythia->GetP(j,2),pythia->GetP(j,3),
                        pythia->GetP(j,4),pythia->GetK(j,2),2);
                } else {
                    // Others
                    pEvent->AddPrtcl(pythia->GetP(j,1),pythia->GetP(j,2),pythia->GetP(j,3),
                        pythia->GetP(j,4),pythia->GetK(j,2),1);
                }
            }
        }
        if ( mode==3 && (muon1<0 || muon2<0 || pythia->GetK(muon1,1)>10 || pythia->GetK(muon2,1)>10) ) {
            cout << "Unexpected behaviour with Zmumu+jets" << endl;
        }
        if ( mode==2 && (gamma<9 || pythia->GetK(gamma,1)>10) ) {
            cout << "Unexpected behaviour with Photon+jets" << endl;
        }

        return true;
    }

    int MuonFinder(TPythia6* pythia, int j) {
        if ( abs(pythia->GetK(j,2))!=13 ) return -1;
        // Final state particle
        if ( pythia->GetK(j,1)<11 ) return j;
        
        for (int probe = pythia->GetK(j,4); probe <= pythia->GetK(j,5); ++probe) {
            if ( abs(pythia->GetK(probe,2)) == 13 ) return probe;
        }
        return -1;
    }

}

# endif
