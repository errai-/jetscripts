#include "HardProc.h"

void HardProc::EventProcessing(Long64_t)
{
    /* Save only the hard process */
    cerr << fPartonList.size() << endl;
    for ( auto i : fPartonList ) {
        fJetEvent->AddJet(fAuxInputs[i].px(),
                          fAuxInputs[i].py(),
                          fAuxInputs[i].pz(),
                          fAuxInputs[i].e(),
                          fJetVars,
                          fWeight,
                          fPDGCode[abs(fAuxInputs[i].user_index())]);
    }

    fOutTree->Fill();
}