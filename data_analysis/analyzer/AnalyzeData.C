#define AnalyzeData_cxx
#include "AnalyzeData.h"

using std::string;
using std::cin;
using std::endl;
using std::vector;
using std::cout;
using std::endl;

void AnalyzeData::Loop(string writeFile)
{
    if (fChain == 0) return;

    /* Scaling for DT, lumi or prescale */
    int lumiScale = 0;

    /* Using of an old pileup file format for DT */
    int oldPU = 0;
    
    /* Variables for testing */
    MemInfo_t memUsage;
    Long64_t initMem = 0;
    Long64_t goodJets = 0;
    Long64_t duplicEvents = 0;
    TStopwatch timer;
    if (testing) {
        gSystem->GetMemInfo(&memUsage);
        initMem = memUsage.fMemUsed;
        timer.Start();
    }

    /* A listing for avoiding duplicate events */
    set<EventID> usedEvents;

    /* Luminosities for each trigger for DT */
    const double triggerLumi[] = {0.000079299, 0.002120549, 0.055696629, 0.261372, 1.062908, 19.711789};

    /* Pileup handling */
    int puCheck = 1;
    vector<TH1D*> pileUpHists;
    for (int looper = 1; looper <= 6; ++looper) {
        std::stringstream tmpName;
        tmpName << "pileuphist" << looper;
        pileUpHists.push_back( new TH1D( tmpName.str().c_str(), "pileup", 50, 0, 50) );
    }

    TFile mcPileUp( (isMC == 1) ? "pileupcorr/PU_PU53RD_10M.root" : 
                                  "pileupcorr/PU_PY53X_10M.root","READ");
    TFile dtPileUp(oldPU ? "pileupcorr/pileup12px_Aug29_mb694.root" : 
                           "pileupcorr/pileup12px_Aug29_mb694_new.root","READ");
    TH1D *mcPiles = 0;
    vector<TH1D*> dtPiles;
    if (isMC == 1 || isMC == 2) {
        mcPiles = (TH1D*)mcPileUp.Get("pileup");
        mcPiles->Scale( 1./mcPiles->Integral() );
        int eTags[] = {40,80,140,200,260,320};
        for (int i=0; i<6; i++) {
            std::stringstream tmpName;
            tmpName << "pileup_jt" << eTags[i];
            dtPiles.push_back( (TH1D*)dtPileUp.Get(tmpName.str().c_str()) );
            dtPiles[i]->Scale(1./dtPiles[i]->Integral() );
        }
    }

    /* Setting up an environment which saves the fractions and the amount of samples */
    Long64_t ptBinAmount = 40;
    double minPt = 56;
    double maxPt = 2000;
    Long64_t etaBinAmount = 50;
    double minEta = -5;
    double maxEta = 5;

    /* 2D histograms */
    ProfileBuilder *analysisHelper = new ProfileBuilder(ptBinAmount,etaBinAmount,
                                                        minPt,maxPt,minEta,maxEta);
    /* Jet energy corrections */
    vector<JetCorrectorParameters> corParams;
    FactorizedJetCorrector *jetECor = EnergyCorrSetup( &corParams, isMC );

    double cor1 = 0, cor2 = 0;
    unsigned counter = 0;
    /* Event loop, until reaching the limit or tree size */
    for (Long64_t jentry=0; jentry<loopLimit; jentry++) {

        if ( jentry%100000 == 0 ) cout << jentry << " events analyzed\n";

        /* Load next entry if there are entries left in the chain */
        Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
        fChain->GetEntry(jentry);

        /* Check that the fraction of missing transverse energy is not too large */
        if (isMC != 3 && PFMet__et_ >= 0.3*PFMet__sumEt_) continue;
        
        /* Duplicate events are ignored */
        if (isMC != 3 && !usedEvents.emplace(mRun, mLumi, mEvent).second) {
            duplicEvents++;
            continue;
        }
        
        assert( PFJets_ < kMaxPFJets_ );
        
        /* Jet loop */
        for (Long64_t kentry=0; kentry<PFJets_; kentry++){
            if (isMC!=3 && !tightID_[kentry]) continue;

            TLorentzVector jetMomentum;
            jetMomentum.SetPxPyPzE(fX[kentry],fY[kentry],fZ[kentry],fT[kentry]);
        
            /* Applying energy corrections (for pythia, corr = 1) */
            double corr = 1;
            if (isMC!= 3) {
                jetECor->setRho( mPFRho );
                jetECor->setJetA( area_[kentry] );
                jetECor->setJetPt( jetMomentum.Pt() );
                jetECor->setJetEta( jetMomentum.Eta() );
                corr = jetECor->getCorrection();
            }
            jetMomentum *= corr;
            //cout << corr << " " << cor_[kentry] << " " << corr/cor_[kentry] << endl;
            cor1 += corr; cor2 += cor_[kentry];
            ++counter;

            double jetPt = jetMomentum.Pt();
            double jetEta = jetMomentum.Eta();

            int triggerType = TriggerType( jetPt );
            if (triggerType == -1) continue;

            if (isDT){
                /* For data, only the triggering jets are of interest */
                for (int triggerIdx = 7*triggerType; triggerIdx!=7*triggerType+7; ++triggerIdx) {
                    if ( TriggerDecision_[triggerIdx]!=1 ) continue;
                    
                    double dtweight = lumiScale ? triggerLumi[5]/triggerLumi[triggerType] : L1_[triggerIdx]*HLT_[triggerIdx];

                    analysisHelper->FillHelper( jetPt, jetEta,chf_[kentry]*betaStar_[kentry], 
                        chf_[kentry]*(1 - betaStar_[kentry]), phf_[kentry], nhf_[kentry] - hf_hf_[kentry],
                        elf_[kentry] + muf_[kentry], hf_hf_[kentry], hf_phf_[kentry], dtweight);

                    pileUpHists[triggerType]->Fill(mNVtxGood, dtweight);
                    goodJets++;
                    break;
                }
            } else if (isMC==1 || isMC==2) {
                if ( jetPt >= 1.5*mPthat ) continue;

                jetMomentum.SetPxPyPzE(gen_fX[kentry],gen_fY[kentry],gen_fZ[kentry],gen_fT[kentry]);
                if ( jetPt >= 1.5*jetMomentum.Pt() ) continue;

                /* Values used to scale the pileup */
                int puBinIdx = dtPiles[triggerType]->FindBin(mTrPu);
                double dtBinCont = dtPiles[triggerType]->GetBinContent(puBinIdx);
                double mcBinCont = mcPiles->GetBinContent(puBinIdx);
                double mcweight =(dtBinCont == 0 || mcBinCont == 0) ? 1. : dtBinCont/mcBinCont;
                mcweight *= pow(15.0f/mPthat , 4.5);

                analysisHelper->FillHelper( jetPt, jetEta, chf_[kentry]*betaStar_[kentry], 
                    chf_[kentry]*(1 - betaStar_[kentry]), phf_[kentry], nhf_[kentry] - hf_hf_[kentry],
                    elf_[kentry] + muf_[kentry], hf_hf_[kentry], hf_phf_[kentry], mcweight);
                pileUpHists[triggerType]->Fill(mNVtxGood, mcweight );
                goodJets++;
            } else {
                analysisHelper->FillHelper( jetPt, jetEta, 0, chf_[kentry], 
                    phf_[kentry], nhf_[kentry], elf_[kentry] + muf_[kentry], 
                    0, 0, fWeight );
                goodJets++;
            }
        }
    }

    if (puCheck) {
        analysisHelper->WriteToFile(&pileUpHists, writeFile);
    }

    delete jetECor;

    if (testing){
        gSystem->GetMemInfo(&memUsage);
        cout << "Memory usage (MB): " << memUsage.fMemUsed - initMem << endl;
        timer.Stop();
        cout << "CPU time used while processing: " << timer.CpuTime() << endl;
        cout << "Duplicate events observed: " << duplicEvents << endl;
        cout << "Triggering jets: " << goodJets << endl;
    }
    cout << cor1/counter << " " << cor2/counter << endl;
    
    if (isMC){
        mcPileUp.Close();
        dtPileUp.Close();
    }
}

