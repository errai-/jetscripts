#include "JetBase.h"


JetBase::JetBase(TTree *tree, 
                 const char *outFile1, 
                 const char *outFile2, 
                 int mode, 
                 int definition ) :
                 
                 fChain         (0), 
                 fMode          (mode), 
                 fDefinition    (definition),
                 fJetCuts       (true),
                 fParamCuts     (true),
                 fInitialized   (true),
                 fAddNonJet     (false),
                 fParticleStudy (false),
                 fR             (0.5),
                 fMinPT         (10.),
                 fSuccessCount  (0)
{
    assert(tree);
    Init(tree);

    /* Tree: autosave every 1 Gb, 10 Mb cache */ 
    fOutFile = new TFile(outFile1, "RECREATE");
    fOutFile->SetCompressionLevel(1);
    fOutTree = new TTree("JetTree","Tree with jet data");
    fOutTree->SetAutoSave(1000000000); fOutTree->SetCacheSize(10000000);  
    TTree::SetBranchStyle(1);
    fOutFileName = outFile2;
    fJetEvent = new JetEvent;
    fJetBranch = fOutTree->Branch("event", &fJetEvent, 32000, 4);
    fJetBranch->SetAutoDelete(kFALSE);
    fOutTree->BranchRef();

    fJetsPerEvent = 2;
    if (mode==0) fJetsPerEvent = 100;
    if (mode==1) fJetsPerEvent = 2;
    if (mode==2||mode==3) fJetsPerEvent = 1;
    if (mode==4) fJetsPerEvent = 100;
    
    /* Fastjet algorithm (settings stated explicitly) */
    fJetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, fR, fastjet::E_scheme, fastjet::Best); 
}


void JetBase::Init(TTree *tree)
{
    /* Set branch addresses and branch pointers */
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("fWeight", &fWeight);
    fChain->SetBranchAddress("fPrtcls", &fPrtcls_);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fX", fX);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fY", fY);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fZ", fZ);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fT", fT);
    fChain->SetBranchAddress("fPrtcls.fPDGCode", fPDGCode);
    fChain->SetBranchAddress("fPrtcls.fAnalysisStatus", fAnalysisStatus);
    fChain->SetBranchAddress("fPrtcls.fHistoryFlavor",fHistoryFlavor);
}


//////////////////////////////////
// Loop over events and particles:
//////////////////////////////////


void JetBase::EventLoop()
{
    /* Sanity checks for the run to be continued */
    if (fChain == 0) return;
    if (!fInitialized)
        return;
    else
        fInitialized = false;

    Long64_t nentries = fChain->GetEntries();
    fTimer.setParams(nentries,2000); fTimer.startTiming();

    for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
        /* Logistics */
        if (jentry!=0&&jentry%2000==0) fTimer.printTime();

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);
        assert( fPrtcls_ < kMaxfPrtcls );
        
        /* Jet clusting and analysis cycle */
        ParticlesToJetsorterInput();
        EventProcessing();
        fJetEvent->Clear();
    }
    
    Finalize();
}

void JetBase::Finalize()
{
    fOutFile = fOutTree->GetCurrentFile();
    fOutTree->AutoSave("Overwrite");
    
    delete fJetEvent;  fJetEvent = 0;
    fOutFile->Close();
}


void JetBase::EventProcessing() {
    fastjet::ClusterSequence fClustSeq(fJetInputs, fJetDef);
    vector< fastjet::PseudoJet > unsorteds = fClustSeq.inclusive_jets( fMinPT );
    fSortedJets = sorted_by_pt( unsorteds );
    
    /* Abort if the event does not meet quality specifications */
    fJetVars.SetZero();
    if (!SelectionParams())
        return;
    
    InitLoop();
    
    if (!JetLoop())
        return;
    
    ++fSuccessCount;
    PostLoop();

    fOutTree->Fill();
}


bool JetBase::JetLoop()
{
    for (size_t i = 0; i < fSortedJets.size(); ++i) {
        if ( i == fJetsPerEvent ) break;

        fJetParts = fSortedJets[i].constituents();

        /* Flavour definition */
        if (fDefinition == 1)
            PhysicsFlavor(i);
        else if (fDefinition == 2)
            GhostPhysicsFlavor(i);
        else if (fDefinition == 3)
            GhostPhysicsFlavor(i);
        else if (fDefinition == 4)
            GhostPhysicsFlavor(i);
        else if (fDefinition == 5)
            PhysAlgoFlavor(i);
        else if (fDefinition == 6)
            PhysAlgoFlavor(i);
        else if (fDefinition == 7)
            HistoricPhysicsFlavor(i);
        else if (fDefinition == 8)
            HadronicFlavor(i);
        else if (fDefinition == 9)
            AlgorithmicFlavor(i);
        else if (fDefinition == 10)
            GhostAlgorithmicFlavor(i);

        if (fParticleStudy)
            ParticleLoop(); /* Operations on jet particles */

        PostProcessing(i);

        fJetEvent->AddJet(fSortedJets[i].px(),
                          fSortedJets[i].py(),
                          fSortedJets[i].pz(),
                          fSortedJets[i].e(),
                          fJetVars,
                          fWeight,
                          fFlavour);
    }
    return true;
}


void JetBase::ParticleLoop()
{  
    fPiPlus   = PseudoJet(); fPiMinus  = PseudoJet();
    fPi0Gamma = PseudoJet(); fGamma    = PseudoJet(); 
    fKaPlus   = PseudoJet(); fKaMinus  = PseudoJet(); 
    fKSZero   = PseudoJet(); fKLZero   = PseudoJet(); 
    fXiMinus  = PseudoJet(); fXiZero   = PseudoJet();
    fOmMinus  = PseudoJet();
    fProton   = PseudoJet(); fAproton  = PseudoJet(); 
    fNeutron  = PseudoJet(); fAneutron = PseudoJet();
    fLambda0  = PseudoJet(); fSigma    = PseudoJet(); 
    fElec     = PseudoJet(); fMuon     = PseudoJet();
    fOthers   = PseudoJet(); fEtSum    = PseudoJet();
    
    for (unsigned int j = 0; j != fJetParts.size(); ++j) {
        if ( fJetParts[j].user_index() < 0 ) continue;
        
        int id = fPDGCode[ fJetParts[j].user_index() ];
        int status = fAnalysisStatus[ fJetParts[j].user_index() ]; 
       
        fEtSum += fJetParts[j];
        if ( id == 211 ) { 
            fPiPlus += fJetParts[j];
        } else if ( id == -211 ) { 
            fPiMinus+= fJetParts[j];
        } else if ( id == 22 ) {
            if ( status == 2 ) {
                fPi0Gamma += fJetParts[j];
            } else {
                fGamma += fJetParts[j];
            }
        } else if ( id == 321 ) { 
            fKaPlus += fJetParts[j];
        } else if ( id == -321 ) { 
            fKaMinus += fJetParts[j];
        } else if ( abs( id ) == 310 ) { 
            fKSZero += fJetParts[j];
        } else if ( abs( id ) == 130 ) { 
            fKLZero += fJetParts[j];
        } else if ( id == 2212 ) { 
            fProton += fJetParts[j];
        } else if ( id == -2212 ) { 
            fAproton += fJetParts[j];
        } else if ( id == 2112 ) { 
            fNeutron += fJetParts[j];
        } else if ( id == -2112 ) { 
            fAneutron += fJetParts[j];
        } else if ( abs( id ) == 3122 ) {
            fLambda0 += fJetParts[j];
        } else if ( abs( id ) == 3112 || abs( id ) == 3222 ) {
            fSigma += fJetParts[j];
        } else if ( abs( id ) == 3312 ) {
            fXiMinus += fJetParts[j];
        } else if ( abs( id ) == 3322 ) {
            fXiZero += fJetParts[j];
        } else if ( abs( id ) == 3334 ) {
            fOmMinus += fJetParts[j];
        } else if ( abs( id ) == 11 ) {
            fElec += fJetParts[j];
        } else if ( abs( id ) == 13 ) {
            fMuon += fJetParts[j];
        } else {
            fOthers += fJetParts[j];   
        }
    }
    fOthers += fXiMinus + fXiZero + fOmMinus;
}


//////////////////////////////////////////////
// Particle addition and selection procedures:
//////////////////////////////////////////////


void JetBase::ParticlesToJetsorterInput()
{
    fJetInputs.clear();
    fLeptons.clear();
    fHardPartons.clear();
    fPartons.clear();
    fTheGamma = PseudoJet();
    fTheLepton = PseudoJet();
    fMET = PseudoJet();

    for (unsigned i = 0; i != fPrtcls_; ++i) {
        fastjet::PseudoJet particleTemp(fX[i],fY[i], fZ[i], fT[i]);
        particleTemp.set_user_index( i ); /* Save particle index */
        
        int stat = fAnalysisStatus[i];
        int pdgID = abs(fPDGCode[i]);

        if (stat==1) {
            /* General final-state particles: neutrinos are thrown into MET */
            
            if (pdgID == 12 || pdgID == 14 || pdgID == 16)
                fMET += particleTemp;
            else
                fJetInputs.push_back(particleTemp);
            
            /* Ttbar lepton+jet events: save leptons */
            if (fMode==4 && (pdgID == 11 || pdgID == 13))
                fLeptons.push_back(particleTemp);
            
        } else if (stat==2) {
            /* Special final-state particles, e.g. muons in Zmumu+jet cases
             * Depending on the event type these can be exlucded from jet
             * clustering and stored to AuxInputs. */
            
            if (fMode==0 || fMode==1) {
                /* pi0 photons in generic events */
                fJetInputs.push_back(particleTemp);
            } else if (fMode==2) {
                /* outgoing photons in gamma+jet events */
                fTheGamma = particleTemp;
            } else if (fMode==3) {
                /* outgoing muons in Zmumu events */
                fLeptons.push_back(particleTemp);
            } else if (fMode==4) {
                /* charged lepton from W stored to output */
                fTheLepton = particleTemp;
                
            if (fAddNonJet && (fMode==2||fMode==3))
                fJetVars.SetZero();
                fJetEvent->AddJet(particleTemp.px(),
                                    particleTemp.py(),
                                    particleTemp.pz(),
                                    particleTemp.e(),
                                    fJetVars,
                                    fWeight,
                                    pdgID);
            }
        } else if (stat==3) {
            /* Outgoing hard process particles - these are used with some of the
             * jet flavour definitions */
            
            /* Always set the hard process partons */
            fHardPartons.push_back(particleTemp);
                
            if (fDefinition==2 || fDefinition==4 || fDefinition==5) {
                /* Physics clustering definition: ghost partons from the hard process */
                particleTemp *= pow( 10, -18 );
                particleTemp.set_user_index( -i-1 );
                fJetInputs.push_back( particleTemp );
            }
        } else if (stat==8) {
            /* Outgoing momentum corrected hard process partons. */
            if (fDefinition==3 || fDefinition==4) {
                /* Final state parton physics definition. */
                particleTemp *= pow( 10, -18 );
                particleTemp.set_user_index( -i-1 );
                fJetInputs.push_back( particleTemp );
            }
        } else if (stat==4) {
            /* Algorithmic and hadronic definition: partons just before hadronization */
            
            if (fDefinition==5 || fDefinition==6 || fDefinition==8 || fDefinition==10) {
                particleTemp *= pow( 10, -18 );
                particleTemp.set_user_index( -i-1 );
                fJetInputs.push_back( particleTemp );
            } else if (fDefinition==9) {
                if (pdgID > 6 && pdgID!=21) continue;
                fPartons.push_back(particleTemp);
            }
        } else if (/*stat == 5 ||*/ stat == 6 || stat == 7 ) {
            /* Hadronic definition hadrons: charm and top omitted for now */
            
            if (fDefinition==8) {
                particleTemp *= pow( 10, -18 );
                particleTemp.set_user_index( -i-1 );
                fJetInputs.push_back( particleTemp );
            }
        }
        /* Unnecessary particles are implicitly discarded */
    }
    /* The MET-"jet" is given a status 10 */
    fJetVars.SetZero();
    if (fAddNonJet)
        fJetEvent->AddJet(fMET.px(),
                          fMET.py(),
                          fMET.pz(),
                          fMET.e(),
                          fJetVars,
                          fWeight,
                          10);
    
    assert( fJetInputs.size() ); /* The input should not be empty */
}


bool JetBase::SelectionParams()
{
    if ( fSortedJets.size() == 0 ) return false;

    if (fMode == 0) {
        
        fJetVars.Alpha = 0;
        fJetVars.DPhi = 0;
        fJetVars.matchPT = 0;
        return true;
        
    } else if (fMode == 1) {
        /* Dijet events: require always at least two jets 
         *
         * Example of a jet selection that should be done:
         *  -Minimum jet pT of 30 GeV.
         *  -Maximum jet eta of 2.5
         *  -Back-to-back angle of min 2.8 rad (2.5 rad)
         *  -A third jet has at most 30% of the average pt of the leading jets (alpha) */
        
        if ( fSortedJets.size() < 2 ) 
            return false;
        
        fJetVars.Alpha   = fSortedJets.size() < 3 ? 0 :
                           2*fSortedJets[2].pt()/fabs(fSortedJets[0].pt()+fSortedJets[1].pt());
        fJetVars.DPhi    = fabs(fSortedJets[0].delta_phi_to( fSortedJets[1] ));
        fJetVars.matchPT = 0;
        
    } else if (fMode == 2) {
        /* Gammajet events: require always sufficient resolution and cuts for gamma eta and pt
         *
         * Example of a selection that should be done
         *  -Minimum jet pT of 30 GeV
         *  -Max jet eta of 2.5
         *  -Back-to-back angle of min 2.8 rad
         *  -A cut for the subleading jet pT with respect to gamma pT (alpha) */
        
        bool gammaInJet = fTheGamma.delta_R( fSortedJets[0] ) < fR;
        bool gammaPT = fTheGamma.pt()<30;
        bool gammaEta = fabs(fTheGamma.eta())>2.5;
        
        if ( gammaInJet || gammaPT || gammaEta )
            return false;
        
        fJetVars.Alpha   = (fSortedJets.size()>1) ? fSortedJets[1].pt()/fTheGamma.pt() : 0;
        fJetVars.DPhi    = fabs(fTheGamma.delta_phi_to(fSortedJets[0]));
        fJetVars.matchPT = fTheGamma.pt();
        
    } else if (fMode == 3) {
        /* Zjet events: require always sufficient resolution and cuts for muon pt and eta
         *
         * Example of a selection that should be done:
         *  -Back-to-back angle of min 2.8 rad
         *  -The subleading jet has smaller than 30% pT compared to the muons. (alpha)
         *  -Min. jet pT of 30GeV
         *  -Max jet eta of 2.5 */
        
        if (fLeptons.size() != 2)
            return false;
        
        bool muonInJet =    fLeptons[0].delta_R(fSortedJets[0])<fR
                         || fLeptons[1].delta_R(fSortedJets[0])<fR; 
        bool muonPT =    (fLeptons[0].pt()<20 || fLeptons[1].pt()<10)
                      && (fLeptons[1].pt()<20 || fLeptons[0].pt()<10);
        bool muonEta =    fabs(fLeptons[0].eta())>2.5 
                       || fabs(fLeptons[1].eta())>2.5; 
                      
        if ( muonInJet || muonPT || muonEta )
            return false;

        /* Dimuon system: require always Z0 mass cut (70-110 GeV) */
        fastjet::PseudoJet tmpVec = fLeptons[0]+fLeptons[1]; 
        if ( fabs(tmpVec.m())<70 || fabs(tmpVec.m())>110 ) return false;
        
        fJetVars.Alpha = (fSortedJets.size() > 1 ) ? fSortedJets[1].pt()/tmpVec.pt() : 0;
        fJetVars.DPhi = tmpVec.delta_phi_to( fSortedJets[0] );
        fJetVars.matchPT = tmpVec.pt();
        
    } else if (fMode == 4) {
        /* ttbar events: require at least 4 jets with sufficient kinematics.
         * Additionally require single-charged lepton events */

        if ( fSortedJets.size() < 4 ) {
            return false;
        }

        if (fJetCuts) {
            for ( auto i = 0u; i < 4; ++i ) {
                if (fSortedJets[i].pt() < 30 || fabs(fSortedJets[i].eta()) > 2.4) {
                    return false;
                }
            }
        }

        /* Check that there is only one charged lepton.
         * Special measures if a "false lepton" passes through the filter. */
        bool unwanted_lepton = false;
        for ( auto lept : fLeptons ) {
            if (lept.pt() > 33 && fabs(lept.eta()) < 2.1) {
                if (!unwanted_lepton)
                    unwanted_lepton = true;
                else
                    return false;
            }
        }

        if (fTheLepton.pt() < 33 || fabs(fTheLepton.eta()) > 2.1) {
            if (unwanted_lepton)
                cerr << "Single-lepton signal" << endl;
            return false;
        } else if (unwanted_lepton) {
            return false;
        }
        
        return true;
    } else {
        throw std::runtime_error("Mode problematic");
    }

    // Cuts are done collectively except for ttbar and generic events
    if (fParamCuts) {
        if (fJetVars.Alpha > 0.3 || fJetVars.DPhi < 2.8)
            return false;
    }
    if (fJetCuts) {
        for (auto i = 0u, N = TMath::Min(unsigned(fJetsPerEvent),unsigned(fSortedJets.size())); i < N; ++i)
            if (fSortedJets[i].pt() < 30 || fabs(fSortedJets[i].eta()) > 2.5)
                return false;
    }

    return true;
}


///////////////////////
// Flavour definitions:
///////////////////////


/* dR_min: the smallest jet axis - parton distance
   dR_nextmin: the next smallest jet axis - parton distance */
void JetBase::PhysicsFlavor(unsigned i) 
{
    fFlavour = 0;
    fJetVars.partonPT = 0;
    double dR_min = 10, dR_nextmin = 10; 
    for ( auto part : fHardPartons ) {
        double dR = fSortedJets[i].delta_R( part );
        
        if ( dR < dR_min ) {
            if (fFlavour!=0) dR_nextmin = dR_min;
            dR_min = dR;
            fJetVars.partonPT = part.pt();
            fFlavour = abs(fPDGCode[part.user_index()]);
        } else if ( dR < dR_nextmin ) {
            dR_nextmin = dR;
        }
    }
    fJetVars.DR = dR_min;
    fJetVars.nextDR = dR_nextmin;
}

/* If there are more than one hard process parton within a jet, mark no flavour */
void JetBase::GhostPhysicsFlavor(unsigned i)
{
    fFlavour = 0;
    int flav_hist = -1;
    for ( auto part : fJetParts ) {
        if (part.user_index() > 0) continue; /* Not ghosts */
        int id = -part.user_index()-1;
        if (flav_hist != -1) { 
            if (fHistoryFlavor[id]==flav_hist) {
                fJetVars.partonPT += part.pt()*pow(10,18);
                fJetVars.partonPT /= 2.0;
            } else {
                flav_hist = -2;
            }
        } else {
            flav_hist = fHistoryFlavor[id];
            fJetVars.partonPT = part.pt()*pow(10,18);
        }
    }
    if (flav_hist >= 0)
        fFlavour = abs(fPDGCode[fHardPartons[flav_hist].user_index()]);
}

void JetBase::PhysAlgoFlavor(unsigned int)
{
    fFlavour = 0;
    int hardflav_hist = -1, algoflav = -1;
    double tmpPt = 0;
    for ( auto part : fJetParts ) {
        if (part.user_index() > 0) continue; /* Not ghosts */
        int id = -part.user_index()-1;
        int status = fAnalysisStatus[id];
        if (status == 3) {
            if (hardflav_hist != -1) { 
                if (fHistoryFlavor[id]==hardflav_hist) {
                    fJetVars.partonPT += part.pt()*pow(10,18);
                    fJetVars.partonPT /= 2.0;
                } else {
                    hardflav_hist = -2;
                }
            } else {
                hardflav_hist = fHistoryFlavor[id];
                fJetVars.partonPT = part.pt()*pow(10,18);
            }
        } else {
            if (tmpPt < part.pt()) {
                algoflav = fPDGCode[id];
                tmpPt = part.pt();
            }
        }
        
    }
    if (hardflav_hist >= 0) {
        fFlavour = abs(fPDGCode[fHardPartons[hardflav_hist].user_index()]);
    } else {
        fFlavour = abs(algoflav);
        if ( fFlavour > 6 && fFlavour != 21 )
            fFlavour = 0;
        fJetVars.partonPT = tmpPt;
    }
}

void JetBase::HistoricPhysicsFlavor(unsigned i) 
{
    map<int,double> et_sums;
    et_sums.emplace(-1,0);
    
    /* Calculate et sum for each outgoing hard process parton descendants */
    for ( auto part : fJetParts ){
        if (part.user_index() < 0) continue; /* Select non-ghosts */
        
        unsigned fl = fHistoryFlavor[part.user_index()];
        if ( !(et_sums.emplace(fl,part.E()).second) )
            et_sums[fl] += part.E();
    }
    
    double et_sum = 0;
    int max_id = 0;
    for ( auto fl : et_sums ) {
        et_sum += fl.second;
        
        if ( fl.second > et_sums[max_id] )
            max_id = fl.first;
    }
    
    /* A purity condition for the jets. Reduces the flavor and parton information
     * to a simple flavor information. */
    fFlavour = (max_id==-1) ? 0 : abs(fPDGCode[fHardPartons[max_id].user_index()]);
    
    fJetVars.DR = et_sums[max_id]/et_sum;
    if ( fJetVars.DR == 0.0 ) cout << "wot" << endl;
    fJetVars.nextDR = 0;
}


void JetBase::HadronicFlavor(unsigned)
{
    fFlavour = 0;
    int partonFlav = 0, hardestLightParton = 0;
    double lightPT = 0, partonPT = 0;
    
    for ( auto part : fJetParts ){
        if (part.user_index() > 0) continue; /* Select ghosts */
            
        int pdgCode = abs(fPDGCode[-part.user_index()-1]);
        int status = fAnalysisStatus[-part.user_index()-1]-2;
        double PT = part.pt();

        if (status == 4 || status == 5) {
            if (status > fFlavour || PT > fJetVars.partonPT) {
                fFlavour = status;
                fJetVars.partonPT = PT;
            }
        } else if (pdgCode == 4 || pdgCode == 5) {
            if (pdgCode > partonFlav || PT > partonPT) {
                partonFlav = pdgCode;
                partonPT = PT;
            }
        } else if (PT > lightPT) {
            hardestLightParton = pdgCode;
            lightPT = PT;
        }
    }
    if (!partonFlav) { 
        partonFlav = hardestLightParton;
        partonPT = lightPT;
    }

    /* flavour is determined with the domination of hadronFlav. If parton flavour
        * is used separately, partonFlav tells the complete parton flavour. */
    if (fFlavour == 0) {
        fFlavour = partonFlav;
        fJetVars.partonPT = partonPT;
    }
}


void JetBase::AlgorithmicFlavor(unsigned i)
{
    fFlavour = 0;
    int hardestLightParton = 0;
    double lightDR = 0, lightPT = 0; 
    
    for ( auto part : fPartons ) {
        double dR = fSortedJets[i].delta_R( part );
        int status = fAnalysisStatus[part.user_index()]-2;
        int pdgCode = abs(fPDGCode[part.user_index()]);

        if (dR > 0.3) continue;

        double PT = part.pt();
        if (status == 4 || status == 5) {
            if (status > fFlavour || PT > fJetVars.partonPT) {
                fFlavour = status;
                fJetVars.DR = dR;
                fJetVars.partonPT = PT;
            }
        } else if (PT > lightPT) {
            hardestLightParton = pdgCode;
            lightDR = dR;
            lightPT = PT;
        }
    }
    
    
    if (!fFlavour) {
        fFlavour = hardestLightParton;
        fJetVars.DR = lightDR;
        fJetVars.partonPT = lightPT;
    }
}


void JetBase::GhostAlgorithmicFlavor(unsigned i)
{
    int hardestLightParton = 0;
    double lightPT = 0;
    fFlavour = 0;
    
    for ( auto part : fJetParts ){
        if (part.user_index() > 0) continue; /* Select ghosts */
            
        int pdgCode = abs(fPDGCode[-part.user_index()-1]);
        double PT = part.pt();

        if (pdgCode == 4 || pdgCode == 5) {
            if (pdgCode > fFlavour || PT > fJetVars.partonPT) {
                fFlavour = pdgCode;
                fJetVars.partonPT = PT;
            }
        } else if (PT > lightPT) {
            hardestLightParton = pdgCode;
            lightPT = PT;
        }
    }
    
    if (!fFlavour) {
        fFlavour = hardestLightParton;
        fJetVars.partonPT = lightPT;
    }
}


///////////////////////////////////////////////////
// Generic functions, these should not be modified:
///////////////////////////////////////////////////



Int_t JetBase::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}


Long64_t JetBase::LoadTree(Long64_t entry)
{
    /* Set the environment to read one entry */
    if (!fChain) return 0;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}


void JetBase::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}
