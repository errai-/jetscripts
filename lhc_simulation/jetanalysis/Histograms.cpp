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

void Histograms::EventProcessing(Long64_t jentry)
{
    JetBase::EventProcessing(jentry);
    fQuarkJetCharge = ChargeSign(fFlavour);
    fOutTree->Fill();
    HistFill(jentry);
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
    FillerHandle( fractionProfilesAll, fSortedJets[i].pt(), fEtSum.Et() );
    if (fFlavour==21) {
        gluonQuark->Fill( fSortedJets[i].pt(), 1);
        FillerHandle( fractionProfilesGluon, fSortedJets[i].pt(), fEtSum.Et() );
    } else if (fFlavour==4 || fFlavour==5) {
        gluonQuark->Fill( fSortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesHQuark, fSortedJets[i].pt(), fEtSum.Et() );
        FillerHandle( fractionProfilesQuark, fSortedJets[i].pt(), fEtSum.Et() );
    } else if (fFlavour==1 || fFlavour==2 || fFlavour==3 ) {
        gluonQuark->Fill( fSortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesLQuark, fSortedJets[i].pt(), fEtSum.Et() );
        FillerHandle( fractionProfilesQuark, fSortedJets[i].pt(), fEtSum.Et() );
    }
}

void Histograms::FillerHandle( vector<TProfile*> &hists, double pt, double scale)
{
    hists[0 ]->Fill( pt, fPiPlus.Et()  /scale ); 
    hists[1 ]->Fill( pt, fPiMinus.Et() /scale );
    hists[2 ]->Fill( pt, fPi0Gamma.Et()/scale ); 
    hists[3 ]->Fill( pt, fKaPlus.Et()  /scale );
    hists[4 ]->Fill( pt, fKaMinus.Et() /scale ); 
    hists[5 ]->Fill( pt, fKSZero.Et()  /scale );
    hists[6 ]->Fill( pt, fKLZero.Et()  /scale ); 
    hists[7 ]->Fill( pt, fProton.Et()  /scale );
    hists[8 ]->Fill( pt, fAproton.Et() /scale ); 
    hists[9 ]->Fill( pt, fNeutron.Et() /scale );
    hists[10]->Fill( pt, fAneutron.Et()/scale ); 
    hists[11]->Fill( pt, fGamma.Et()   /scale );
    hists[12]->Fill( pt, fLambda0.Et() /scale ); 
    hists[13]->Fill( pt, fSigma.Et()   /scale );
    hists[14]->Fill( pt, (fElec.Et()+fMuon.Et())/scale ); 
    hists[15]->Fill( pt, fOthers.Et()  /scale );
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