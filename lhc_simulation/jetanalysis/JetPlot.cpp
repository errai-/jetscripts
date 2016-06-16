#include "JetPlot.h"


void JetPlot::PostProcessing()
{
    fJetVars.constituents = 1;
    int fl = -fFlavour-1;

    for ( auto part : fJetParts ) {
        if ( part.user_index() < 0 )
            continue;

        fJetEvent->AddJet(part.px(),
                          part.py(),
                          part.pz(),
                          part.e(),
                          fJetVars,
                          fWeight,
                          fl);
    }

    fJetVars.constituents = fJetParts.size();
}

void JetPlot::InitLoop()
{
    return;
}

void JetPlot::PostLoop()
{
    for (unsigned i = 0; i != fPrtcls_; ++i) {
        int status = fAnalysisStatus[i];
        int hFlav = fHistoryFlavor[i];
        int flavId = -fFlavour-1;
        if (status>4 && status!=8) continue;

        fastjet::PseudoJet part(fX[i],fY[i], fZ[i], fT[i]);

        if (status==1) {
            if (hFlav>=0)
                flavId = -1;
            else
                flavId = -2;
            fJetVars.constituents = 2;
        } else if (status==2) {
            fJetVars.constituents = 9;
        } else if (status==3) {
            if (fPDGCode[i]>0) {
                fJetVars.constituents = 3;
                flavId = -fPDGCode[i]-1;
            } else {
                fJetVars.constituents = 4;
                flavId = -abs(fPDGCode[i])-1;
            }
        } else if (status==8) {
            fJetVars.constituents = 5;
            flavId = -abs(fPDGCode[i])-1;
        } else if (status==4) {
            if (hFlav>=0)
                fJetVars.constituents = 6;
                fJetEvent->AddJet(part.px(),
                                  part.py(),
                                  part.pz(),
                                  part.e(),
                                  fJetVars,
                                  fWeight,
                                  flavId);

            if (fPDGCode[i]>0) {
                fJetVars.constituents = 7;
                flavId = -fPDGCode[i]-1;
            } else {
                fJetVars.constituents = 8;
                flavId = -abs(fPDGCode[i])-1;
            }
        } else {
            continue;
        }

        fJetEvent->AddJet(part.px(),
                          part.py(),
                          part.pz(),
                          part.e(),
                          fJetVars,
                          fWeight,
                          flavId);
        fJetVars.constituents = 0;
    }
}
