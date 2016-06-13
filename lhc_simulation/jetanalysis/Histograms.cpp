#include "Histograms.h"


void Histograms::InitFP(){
    for (int idx = 0; idx != 16; ++idx){
        std::stringstream tmpString("");
        tmpString << "g" << idx;
        fractionProfilesGluon.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
        tmpString.str("");
        tmpString << "q" << idx;
        fractionProfilesQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
        tmpString.str("");
        tmpString << "lq" << idx;
        fractionProfilesLQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
        tmpString.str("");
        tmpString << "hq" << idx;
        fractionProfilesHQuark.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
        tmpString.str("");
        tmpString << "a" << idx;
        fractionProfilesAll.push_back(new TProfile(tmpString.str().c_str(),"",ptBins,ptRange));
    }
}

void Histograms::PostProcessing(unsigned i)
{
    for (auto jet = 0u, N = TMath::Min(unsigned(fJetsPerEvent),unsigned(fSortedJets.size())); jet < N; ++jet)
        if (fabs(fSortedJets[jet].eta()) > 1.3)
            return;
    HistFill(i);
}

void Histograms::Finalize()
{
    WriteResults();
    JetBase::Finalize();
}

/* The charge sign of a quark jet */
int Histograms::ChargeSign( int id )
{
    if ( id == 1 ) return 1;
    if ( id == -2 ) return 1;
    if ( id == -3 ) return 1;
    if ( id == 4 ) return 1;
    if ( id == -5 ) return 1;
    if ( id == 6 ) return 1;
    if ( id == -1 ) return -1;
    if ( id == 2 ) return -1;
    if ( id == 3 ) return -1;
    if ( id == -4 ) return -1;
    if ( id == 5 ) return -1;
    if ( id == -6 ) return -1;
    return 1;
}

void Histograms::HistFill(int i)
{
    if ( fOpposite ) {
        if ( i == 0 )
            i = 1;
        else if ( 1 == 1 )
            i = 0;
        else
            return;
    } else if ( i > 1 ) {
        return;
    }
    
    FillerHandle( fractionProfilesAll, fSortedJets[i].pt(), fEtSum.E() );
    if (fFlavour==21) {
        gluonQuark->Fill( fSortedJets[i].pt(), 1);
        FillerHandle( fractionProfilesGluon, fSortedJets[i].pt(), fEtSum.E() );
    } else if (fFlavour==4 || fFlavour==5) {
        gluonQuark->Fill( fSortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesHQuark, fSortedJets[i].pt(), fEtSum.E() );
        FillerHandle( fractionProfilesQuark, fSortedJets[i].pt(), fEtSum.E() );
    } else if (fFlavour==1 || fFlavour==2 || fFlavour==3 ) {
        gluonQuark->Fill( fSortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesLQuark, fSortedJets[i].pt(), fEtSum.E() );
        FillerHandle( fractionProfilesQuark, fSortedJets[i].pt(), fEtSum.E() );
    }
}

void Histograms::FillerHandle( vector<TProfile*> &hists, double pt, double scale)
{
    hists[0 ]->Fill( pt, fPiPlus.E()  /scale, fWeight );
    hists[1 ]->Fill( pt, fPiMinus.E() /scale, fWeight );
    hists[2 ]->Fill( pt, fPi0Gamma.E()/scale, fWeight );
    hists[3 ]->Fill( pt, fKaPlus.E()  /scale, fWeight );
    hists[4 ]->Fill( pt, fKaMinus.E() /scale, fWeight );
    hists[5 ]->Fill( pt, fKSZero.E()  /scale, fWeight );
    hists[6 ]->Fill( pt, fKLZero.E()  /scale, fWeight );
    hists[7 ]->Fill( pt, fProton.E()  /scale, fWeight );
    hists[8 ]->Fill( pt, fAproton.E() /scale, fWeight );
    hists[9 ]->Fill( pt, fNeutron.E() /scale, fWeight );
    hists[10]->Fill( pt, fAneutron.E()/scale, fWeight );
    hists[11]->Fill( pt, fGamma.E()   /scale, fWeight );
    hists[12]->Fill( pt, fLambda0.E() /scale, fWeight );
    hists[13]->Fill( pt, fSigma.E()   /scale, fWeight );
    hists[14]->Fill( pt, (fElec+fMuon).E()/scale, fWeight );
    hists[15]->Fill( pt, fOthers.E()  /scale, fWeight );
}

void Histograms::WriteResults()
{
    fOutFile2 = new TFile(fOutFileName.c_str(),"RECREATE");
    fOutFile2->cd();

    TH1D *gq = gluonQuark->ProjectionX("gluonvsquark","");
    gq->Write();
    
    for (unsigned int i = 0; i != fractionProfilesGluon.size(); ++i){
        fractionProfilesGluon[i]->Write();
        fractionProfilesQuark[i]->Write();
        fractionProfilesLQuark[i]->Write();
        fractionProfilesHQuark[i]->Write();
        fractionProfilesAll[i]->Write();
    }
    
    fOutFile2->Close();
}
