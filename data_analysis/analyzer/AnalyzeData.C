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

    /* Event loop, until reaching the limit or tree size */
    for (Long64_t jentry=0; jentry<loopLimit; jentry++) {

        if ( jentry%100000 == 0 ) cout << jentry << " events analyzed\n";

        /* Load next entry if there are entries left in the chain */
        Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
        fChain->GetEntry(jentry);

        /* Check that the fraction of missing transverse energy is not too large */
        if (isMC != 3 && PFMet__et_ >= 0.3*PFMet__sumEt_) continue;
        
        /* Duplicate events are ignored */
        if (!usedEvents.emplace(EvtHdr__mRun, EvtHdr__mLumi, EvtHdr__mEvent).second) {
            duplicEvents++;
            continue;
        }
        
        assert( PFJets__ < kMaxPFJets_ );
        
        /* Jet loop */
        for (Long64_t kentry=0; kentry<PFJets__; kentry++){
            if (!PFJets__tightID_[kentry]) continue;

            TLorentzVector jetMomentum;
            jetMomentum.SetPxPyPzE(PFJets__P4__fCoordinates_fX[kentry],
                                   PFJets__P4__fCoordinates_fY[kentry],
                                   PFJets__P4__fCoordinates_fZ[kentry],
                                   PFJets__P4__fCoordinates_fT[kentry]);
        
            /* Applying energy corrections (for pythia, corr = 1) */
            double corr = 1;     
            if (isMC != 3){
                jetECor->setRho( EvtHdr__mPFRho );
                jetECor->setJetA( PFJets__area_[kentry] );
                jetECor->setJetPt( jetMomentum.Pt() );
                jetECor->setJetEta( jetMomentum.Eta() );
                corr = jetECor->getCorrection();
            }
            jetMomentum *= corr;

            double jetPt = jetMomentum.Pt();
            double jetEta = jetMomentum.Eta();

            int triggerType = TriggerType( jetPt );
            if (triggerType == -1) continue;

            if (isDT){
                /* For data, only the triggering jets are of interest */
                for (int triggerIdx = 7*triggerType; triggerIdx!=7*triggerType+7; ++triggerIdx) {
                    if ( TriggerDecision_[triggerIdx]!=1 ) continue;
                    
                    analysisHelper->FillHelper( jetPt, jetEta,
                        PFJets__chf_[kentry]*PFJets__betaStar_[kentry], 
                        PFJets__chf_[kentry]*(1 - PFJets__betaStar_[kentry]), 
                        PFJets__phf_[kentry], PFJets__nhf_[kentry] - PFJets__hf_hf_[kentry],
                        PFJets__elf_[kentry] + PFJets__muf_[kentry],
                        PFJets__hf_hf_[kentry], PFJets__hf_phf_[kentry],
                        1.0, lumiScale ? triggerLumi[5]/triggerLumi[triggerType] : 
                        L1Prescale_[triggerIdx]*HLTPrescale_[triggerIdx]);

                    pileUpHists[triggerType]->Fill(EvtHdr__mNVtxGood);
                    goodJets++;
                    break;
                }
            } else if (isMC==1 || isMC==2) {
                if ( jetPt >= 1.5*EvtHdr__mPthat ) continue;
                
                jetMomentum.SetPxPyPzE(PFJets__genP4__fCoordinates_fX[kentry],
                                    PFJets__genP4__fCoordinates_fY[kentry],
                                    PFJets__genP4__fCoordinates_fZ[kentry],
                                    PFJets__genP4__fCoordinates_fT[kentry]);
                if ( jetPt >= 1.5*jetMomentum.Pt() ) continue;
        
                /* Values used to scale the pileup */
                int puBinIdx = dtPiles[triggerType]->FindBin(EvtHdr__mTrPu);
                double dtBinCont = dtPiles[triggerType]->GetBinContent(puBinIdx);
                double mcBinCont = mcPiles->GetBinContent(puBinIdx);
                
                analysisHelper->FillHelper( jetPt, jetEta,
                    PFJets__chf_[kentry]*PFJets__betaStar_[kentry], 
                    PFJets__chf_[kentry]*(1 - PFJets__betaStar_[kentry]), 
                    PFJets__phf_[kentry], PFJets__nhf_[kentry] - PFJets__hf_hf_[kentry],
                    PFJets__elf_[kentry] + PFJets__muf_[kentry], 
                    PFJets__hf_hf_[kentry], PFJets__hf_phf_[kentry],
                    (dtBinCont == 0 || mcBinCont == 0) ? 1. : dtBinCont/mcBinCont );
                pileUpHists[triggerType]->Fill(EvtHdr__mNVtxGood, 
                    (dtBinCont == 0 || mcBinCont == 0) ? 1. : dtBinCont/mcBinCont );
                goodJets++;
            } else {
                analysisHelper->FillHelper( jetPt, jetEta,
                    0, 
                    PFJets__chf_[kentry], 
                    PFJets__phf_[kentry], PFJets__nhf_[kentry] - PFJets__hf_hf_[kentry],
                    PFJets__elf_[kentry] + PFJets__muf_[kentry], 
                    PFJets__hf_hf_[kentry], PFJets__hf_phf_[kentry],
                    1 );
                pileUpHists[triggerType]->Fill(EvtHdr__mNVtxGood); 
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
    
    if (isMC){
        mcPileUp.Close();
        dtPileUp.Close();
    }
}

