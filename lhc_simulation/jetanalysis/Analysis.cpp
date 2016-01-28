#include "Analysis.h"


void Analysis::PostProcessing()
{
    Cuts();
    fJetVars.constituents = fCutJetParts.size();
    fJetVars.PTD = PTD();
    fJetVars.Sigma2 = Sigma2();

    /* Save jet properties*/
    TypeSort();
}


/* Throw the obtained values into temporary containers.
   For obtaining a good scaling summation is non-trivial. */
void Analysis::TypeSort()
{
    PseudoJet tmpLorentz;
    double cumulator = 0;
    
    cumulator += fPiPlus.E();
    tmpLorentz += fPiPlus;
    cumulator += fPiMinus.E();
    tmpLorentz += fPiMinus;
    cumulator += fKaPlus.E();
    tmpLorentz += fKaPlus;
    cumulator += fKaMinus.E();
    tmpLorentz += fKaMinus;
    cumulator += fProton.E();
    tmpLorentz += fProton;
    cumulator += fAproton.E();
    tmpLorentz += fAproton;
    cumulator += fSigma.E();
    tmpLorentz += fSigma;
    cumulator += fXiMinus.E();
    tmpLorentz += fXiMinus;
    cumulator += fOmMinus.E();
    tmpLorentz += fOmMinus;
    fJetVars.chf = cumulator;
    fJetVars.chm = tmpLorentz.m();
    
    tmpLorentz = PseudoJet();
    cumulator = 0;
    
    cumulator += fKSZero.E();
    tmpLorentz += fKSZero;
    cumulator += fKLZero.E();
    tmpLorentz += fKLZero;
    cumulator += fNeutron.E();
    tmpLorentz += fNeutron;
    cumulator += fAneutron.E();
    tmpLorentz += fAneutron;
    cumulator += fLambda0.E();
    tmpLorentz += fLambda0;
    cumulator += fXiZero.E();
    tmpLorentz += fXiZero;
    fJetVars.nhf = cumulator;
    fJetVars.nhm = tmpLorentz.m();
    
    tmpLorentz = PseudoJet();
    cumulator = 0;
    
    cumulator += fPi0Gamma.E();
    tmpLorentz += fPi0Gamma;
    cumulator += fGamma.E();
    tmpLorentz += fGamma;
    fJetVars.phf = cumulator;
    fJetVars.phm = tmpLorentz.m();
    
    fJetVars.elf = fElec.E();
    fJetVars.elm = fElec.m();
    
    fJetVars.muf = fMuon.E();
    fJetVars.mum = fMuon.m();
    
    double scale = fJetVars.chf + fJetVars.nhf + fJetVars.phf + fJetVars.elf + fJetVars.muf;
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
    vector<fastjet::PseudoJet> tmpParts;
    bool cutMode = false;

    if (cutMode) {
        /* Explicit cuts (pt cut for photons and neutral hadrons) */
        for ( auto q : fJetParts ) {
            if ( q.user_index() < 0) continue;
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

