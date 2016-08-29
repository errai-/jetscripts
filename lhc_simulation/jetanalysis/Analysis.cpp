#include "Analysis.h"


void Analysis::PostProcessing(unsigned)
{
    Cuts();
    fJetVars.constituents = fCutJetParts.size();
    fJetVars.PTD = PTD();
    fJetVars.Sigma2 = Sigma2();

    /* Save jet properties*/
    TypeSort();
}


/* Throw the obtained values into temporary containers.
   For obtaining a good scaling summation is non-trivial.
   If Et is used instead of E, there will be error in the
   third or fourth decimal. */
void Analysis::TypeSort()
{
    PseudoJet tmpLorentz;

    tmpLorentz += fPiPlus;
    tmpLorentz += fPiMinus;
    tmpLorentz += fKaPlus;
    tmpLorentz += fKaMinus;
    tmpLorentz += fProton;
    tmpLorentz += fAproton;
    tmpLorentz += fSigma;
    tmpLorentz += fXiMinus;
    tmpLorentz += fOmMinus;
    fJetVars.chf = tmpLorentz.E();

    tmpLorentz = PseudoJet();

    tmpLorentz += fKSZero;
    tmpLorentz += fKLZero;
    tmpLorentz += fNeutron;
    tmpLorentz += fAneutron;
    tmpLorentz += fLambda0;
    tmpLorentz += fXiZero;
    fJetVars.nhf = tmpLorentz.E();

    tmpLorentz = PseudoJet();

    tmpLorentz += fPi0Gamma;
    tmpLorentz += fGamma;
    fJetVars.phf = tmpLorentz.E();

    fJetVars.elf = fElec.E();

    fJetVars.muf = fMuon.E();

    double scale = fEtSum.E();
    fJetVars.chf /= scale;
    fJetVars.nhf /= scale;
    fJetVars.phf /= scale;
    fJetVars.elf /= scale;
    fJetVars.muf /= scale;
}


///////////////
// qgl-studies:
///////////////


double Analysis::PTD()
{
    if (fMode==0) return 0;
    double square = 0, linear = 0;
    for(size_t q = 0; q != fCutJetParts.size(); ++q) {
        square += pow(fCutJetParts[q].pt(),2);
        linear += fCutJetParts[q].pt();
    }
    return sqrt(square)/linear;
}

double Analysis::Sigma2()
{
    if (fMode==0) return 0;
    double weightedDiffs[4] = {0,0,0,0};
    double phi = 0, eta = 0, pT2 = 0;

    for(size_t q = 0; q != fCutJetParts.size(); ++q) {
        pT2 += pow(fCutJetParts[q].pt(),2);
        eta += pow(fCutJetParts[q].pt(),2)*fCutJetParts[q].eta();
        phi += pow(fCutJetParts[q].pt(),2)*fCutJetParts[q].phi();
    }
    eta /= pT2; phi /= pT2;

    for(unsigned int q = 0; q != fCutJetParts.size(); ++q)
    {
        double deltaEta = eta-fCutJetParts[q].eta();
        double deltaPhi = TVector2::Phi_mpi_pi( phi-fCutJetParts[q].phi() );
        double pT2Tmp = pow(fCutJetParts[q].pt(),2);
        weightedDiffs[0] += pT2Tmp*deltaEta*deltaEta;
        weightedDiffs[3] += pT2Tmp*deltaPhi*deltaPhi;
        weightedDiffs[1] -= pT2Tmp*deltaEta*deltaPhi;
    }
    weightedDiffs[2] = weightedDiffs[1];

    TMatrixDSymEigen me( TMatrixDSym(2,weightedDiffs) );
    TVectorD eigenvals = me.GetEigenValues();

    return sqrt(eigenvals[1]/pT2);
}

void Analysis::Cuts()
{
    fCutJetParts.clear();
    bool cutMode = true;

    if (cutMode) {
        /* Explicit cuts (pt cut for photons and neutral hadrons) */
        vector<fastjet::PseudoJet> tmpParts;
        for ( auto q : fJetParts ) {
            if (q.user_index() < 0) continue;
            int id = abs(fPDGCode[ q.user_index() ]);
            if (!( q.pt()<1 && (id == 22 || (IsHadron(id) && !IsCharged(id)))) )
                tmpParts.push_back( q );
        }

        /* Implicit cuts (pt cut for hadrons) */
        for ( auto q : tmpParts ) {
            int id = abs(fPDGCode[ q.user_index() ]);
            if ( !IsHadron(id) || ( (IsCharged(id) && q.pt()>0.3) ||
                (!IsCharged(id) && q.pt()>3) ) )
            {
                fCutJetParts.push_back( q );
            }
        }
    } else {
        for ( auto q : fJetParts ) {
            if (q.user_index() < 0)
                continue;
            fCutJetParts.push_back(q);
        }
    }
}


bool Analysis::IsHadron(int pdg)
{
    if(pdg>99) return true;
    return false;
}

bool Analysis::IsCharged(int pdg)
{
    int charge = 0;
    /* photons and neutrinos */
    if (pdg==22 || pdg==12 || pdg==14 ||pdg==16 ) return false;
    /* charged leptons */
    if (pdg==11 || pdg==13 || pdg==15 ) return true;

    pdg = (pdg/10)%1000;
    if (pdg < 100) { /* Mesons */
        if ((pdg%10)%2 == 0) { charge += 2; }
        else { charge -= 1; }

        if ((pdg/10)%2 == 0) { charge -= 2; }
        else { charge += 1; }

    } else { /* Baryons */
        while (pdg != 0) {
            int digit = pdg%10;
            pdg = pdg/10;
            if (digit%2 == 0) { charge += 2; }
            else { charge -= 1; }
        }
    }
    if (charge == 0) return false;
    else return true;
}

