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
                 fAddNonJet     (true),
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
    if (mode==1) fJetsPerEvent = 100;
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
        assert(fPrtcls_ < kMaxfPrtcls);

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
    vector<fastjet::PseudoJet> unsorted = fClustSeq.inclusive_jets(fMinPT);
    fSortedJets = sorted_by_pt(unsorted);

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
        if (i == fJetsPerEvent || fSortedJets[i].pt()<30) break;
        fJetVars.SetZero();

        fJetParts = fSortedJets[i].constituents();

        /* Flavour definition */
        fFlavour = 0;
        fJetVars.partonPT = 0;
        if (fDefinition == 1)
            PhysFlavor(i);
        else if (fDefinition == 2)
            GhostFlavor();
        else if (fDefinition == 3)
            GhostFlavor();
        else if (fDefinition == 4)
            GhostFlavor();
        else if (fDefinition == 5)
            EqualFlavor();
        else if (fDefinition == 6)
            LOFSFlavor();
        else if (fDefinition == 7)
            GhostFlavor();
        else if (fDefinition == 8)
            HadrFlavor();
        else if (fDefinition == 9)
            AlgoFlavor(i);
        else if (fDefinition == 10)
            FSFlavor();
        else if (fDefinition == 11)
            GhostFlavor();

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
        if (id == 211) {
            fPiPlus += fJetParts[j];
        } else if (id == -211) {
            fPiMinus+= fJetParts[j];
        } else if (id == 22) {
            if (status == 2) {
                fPi0Gamma += fJetParts[j];
            } else {
                fGamma += fJetParts[j];
            }
        } else if (id == 321) {
            fKaPlus += fJetParts[j];
        } else if (id == -321) {
            fKaMinus += fJetParts[j];
        } else if (abs(id) == 310) {
            fKSZero += fJetParts[j];
        } else if (abs(id) == 130) {
            fKLZero += fJetParts[j];
        } else if (id == 2212) {
            fProton += fJetParts[j];
        } else if (id == -2212) {
            fAproton += fJetParts[j];
        } else if (id == 2112) {
            fNeutron += fJetParts[j];
        } else if (id == -2112) {
            fAneutron += fJetParts[j];
        } else if (abs(id) == 3122) {
            fLambda0 += fJetParts[j];
        } else if (abs(id) == 3112 || abs(id) == 3222) {
            fSigma += fJetParts[j];
        } else if (abs(id) == 3312) {
            fXiMinus += fJetParts[j];
        } else if (abs(id) == 3322) {
            fXiZero += fJetParts[j];
        } else if (abs(id) == 3334) {
            fOmMinus += fJetParts[j];
        } else if (abs(id) == 11) {
            fElec += fJetParts[j];
        } else if (abs(id) == 13) {
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
    fGhosts.clear();
    fAux.clear();
    fTheGamma = PseudoJet();
    fTheLepton = PseudoJet();
    fMET = PseudoJet();

    for (unsigned i = 0; i != fPrtcls_; ++i) {
        fastjet::PseudoJet particleTemp(fX[i],fY[i], fZ[i], fT[i]);
        particleTemp.set_user_index(i); /* Save particle index */

        int stat = fAnalysisStatus[i];
        int pdgID = abs(fPDGCode[i]);

        if (stat==1) {
            /* General final-state particles: neutrinos are thrown into MET */
            if (pdgID==12 || pdgID==14 || pdgID==16)
                fMET += particleTemp;
            else
                fJetInputs.push_back(particleTemp);

            /* In gamma+jets events all photons are monitored */
            if (fMode==2 && pdgID==22)
                fAux.push_back(particleTemp);

            /* In Zmumu+jets events all muons are monitored */
            if (fMode==3 && pdgID==13)
                fAux.push_back(particleTemp);

            /* Ttbar lepton+jet events: save charged leptons */
            if (fMode==4 && (pdgID==11 || pdgID==13))
                fAux.push_back(particleTemp);

            if (fDefinition==11 && (pdgID==11 || pdgID==13 || pdgID==15 || pdgID==22)) {
                particleTemp.set_user_index(fPDGCode[i]);
                fPartons.push_back(particleTemp);
            }
        } else if (stat==2) {
            /* Pi0 photons */
            fJetInputs.push_back(particleTemp);

        } else if (stat==3 || stat==4 || stat==5) {
            /* Hadronic definition hadrons: (strange), charm and bottom */
            if (fDefinition==8) {
                particleTemp *= pow(10, -18);
                particleTemp.set_user_index(-i-1);
                fJetInputs.push_back(particleTemp);
            }

        } else if (stat==6) {
            /* Algorithmic and hadronic definition: partons just before hadronization */
            if (fDefinition==5 || fDefinition==6 || fDefinition==8 || fDefinition==10) {
                particleTemp *= pow(10, -18);
                particleTemp.set_user_index(-i-1);
                fJetInputs.push_back(particleTemp);
            } else if (fDefinition==9) {
                if (pdgID>6 && pdgID!=21) continue;
                fPartons.push_back(particleTemp);
            } else if (fDefinition==11) {
                if (pdgID>6 && pdgID!=21) continue;
                particleTemp.set_user_index(fPDGCode[i]);
                fPartons.push_back(particleTemp);
            }
        } else if (stat==7) {
            /* Outgoing momentum corrected hard process partons and charged leptons. */
            if (fDefinition==3 || fDefinition==4) {
                /* Final state parton physics definition. */
                particleTemp.set_user_index(pdgID);
                fGhosts.push_back(particleTemp);
                particleTemp *= pow(10,-18);
                particleTemp.set_user_index(-fGhosts.size());
                fJetInputs.push_back(particleTemp);
            }

            if (pdgID == 12 || pdgID == 14 || pdgID == 16) {
                fMET += particleTemp;
            }

            /* Add gamma/lepton in certain event types as a jet.
             * Status: lepton or gamma particle code. */
            if (fAddNonJet && (fMode==2||fMode==3||fMode==4)) {
                fJetVars.SetZero();
                fJetEvent->AddJet(particleTemp.px(),
                                  particleTemp.py(),
                                  particleTemp.pz(),
                                  particleTemp.e(),
                                  fJetVars,
                                  fWeight,
                                  fabs(pdgID) );
            } else if (fMode==2) {
                /* outgoing photons in gamma+jet events */
                fTheGamma = particleTemp;
            } else if (fMode==3) {
                /* outgoing muons in Zmumu events */
                fLeptons.push_back(particleTemp);
            } else if (fMode==4) {
                /* charged lepton from W stored to output */
                fTheLepton = particleTemp;
            }
        } else if (stat==8) {
            /* Special final-state particles, e.g. muons in Zmumu+jet cases
             * Depending on the event type these can be excluded from jet
             * clustering and stored to AuxInputs. */
            if (pdgID<=6 || pdgID==21) {
                /* Always set the hard process partons */
                fHardPartons.push_back(particleTemp);

                if (fDefinition==2 || fDefinition==4 || fDefinition==5) {
                    /* Physics clustering definition: ghost partons from the hard process */
                    if (fDefinition==5) {
                        particleTemp.set_user_index(-i-1);
                    } else {
                        particleTemp.set_user_index(-pdgID-1);
                        fGhosts.push_back(particleTemp);
                        particleTemp.set_user_index(-fGhosts.size());
                    }

                    particleTemp *= pow(10,-18);
                    fJetInputs.push_back(particleTemp);
                }
            }
            /* For now, ignore the uncorrected photons and muons */
        }
    }

    /* The MET-"jet" is given a status 10 */
    fJetVars.SetZero();
    if (fAddNonJet) {
        fJetEvent->AddJet(fMET.px(),
                          fMET.py(),
                          fMET.pz(),
                          fMET.e(),
                          fJetVars,
                          fWeight,
                          10);
    }

    /* QCD Aware operations */
    if (fDefinition==11) {
        fHardPartons.clear();

        AntiKtMeasure *aktm = new AntiKtMeasure(0.4);
        QCDAwarePlugin *qa_akt = new QCDAwarePlugin(aktm);
        fastjet::ClusterSequence qcdaw(fPartons, qa_akt);
        vector<PseudoJet> partonJets = sorted_by_pt(qcdaw.inclusive_jets());
        for (auto pJet : partonJets) {
            fGhosts.push_back(pJet);
            pJet *= pow(10,-18);
            pJet.set_user_index(fGhosts.size());
            fJetInputs.push_back(pJet);
        }
    }

    assert(fJetInputs.size()); /* The input should not be empty */
}


bool JetBase::SelectionParams()
{
    if ( fSortedJets.size() == 0 ) return false;

    unsigned studyJets = TMath::Min(unsigned(fJetsPerEvent),unsigned(fSortedJets.size()));
    if (fMode==0) {
        fJetVars.Alpha = 0;
        fJetVars.DPhi = 0;
        fJetVars.matchPT = 0;
        return true;

    } else if (fMode==1) {
        /* Dijet events: require always at least two jets
         *
         * Example of a jet selection that should be done:
         *  -Minimum jet pT of 30 GeV.
         *  -Maximum jet eta of 2.5
         *  -Back-to-back angle of min 2.8 rad (2.5 rad)
         *  -A third jet has at most 30% of the average pt of the leading jets (alpha) */
        studyJets = 2;

        if (fSortedJets.size()<2)
            return false;

        fJetVars.Alpha   = (fSortedJets.size() < 3) ? 0 :
                           2*fSortedJets[2].pt()/fabs(fSortedJets[0].pt()+fSortedJets[1].pt());
        fJetVars.DPhi    = fabs(fSortedJets[0].delta_phi_to( fSortedJets[1] ));
        fJetVars.matchPT = 0;

    } else if (fMode==2) {
        /* Gammajet events: require always sufficient resolution and cuts for gamma eta and pt
         *
         * Example of a selection that should be done
         *  -Minimum jet pT of 30 GeV
         *  -Max jet eta of 2.5
         *  -Back-to-back angle of min 2.8 rad
         *  -A cut for the subleading jet pT with respect to gamma pT (alpha) */
        studyJets = 1;

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
         *  -Min. jet pT of 30 GeV
         *  -Max jet eta of 2.5 */
        studyJets = 1;

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
        studyJets = 4;

        if (fSortedJets.size() < 4) {
            return false;
        }
        /* Quality check for primary lepton */
        if (fTheLepton.pt() < 33 || fabs(fTheLepton.eta()) > 2.1 || !Isolation(fTheLepton))
            return false;

        /* Quality check for secondary leptons (required to be low-quality) */
        for ( auto lept : fLeptons ) {
            if (lept.pt()>33 && fabs(lept.eta())<2.1 && Isolation(lept)) {
                cerr << "High-quality secondary lepton" << endl;
                return false;
            }
        }

        /* Treated as dummies here for now */
        fJetVars.Alpha = 0;
        fJetVars.DPhi = 3;
    } else {
        throw std::runtime_error("Mode problematic");
    }

    // Cuts are done collectively except for ttbar and generic events
    if (fParamCuts) {
        if (fJetVars.Alpha > 0.3 || fJetVars.DPhi < 2.8)
            return false;
    }

    if (fJetCuts) {
        for (auto i = 0u; i < studyJets; ++i)
            if (fSortedJets[i].pt() < 30 || fabs(fSortedJets[i].eta()) > 2.5)
                return false;
    }

    return true;
}


Bool_t JetBase::Isolation(PseudoJet prt, double R)
{
    unsigned id = prt.user_index();
    unsigned type = abs(fPDGCode[id]);
    double limit = 0.1;
    if (type==11) {
        R = 0.3;
        limit = 0.1;
    } else if (type==13) {
        R = 0.4;
        limit = 0.12;
    }

    double E_dR = 0;
    for (auto part : fJetInputs) {
        if (part.user_index()<0 || part.user_index()==id)
            continue;
        double dR = prt.delta_R(part);
        if (dR < R)
            E_dR += part.e();
    }

    if (E_dR/prt.e() > limit)
        return false;
    else
        return true;
}


///////////////////////
// Flavour definitions:
///////////////////////


/* dR_min: the smallest jet axis - parton distance
   dR_nextmin: the next smallest jet axis - parton distance */
void JetBase::PhysFlavor(unsigned jet)
{
    double dR_min = 10, dR_nextmin = 10;
    for (auto part : fHardPartons) {
        double dR = fSortedJets[jet].delta_R( part );

        if (dR < dR_min) {
            if (fFlavour!=0) dR_nextmin = dR_min;
            dR_min = dR;
            fJetVars.partonPT = part.pt();
            fFlavour = abs(fPDGCode[part.user_index()]);
        } else if ( dR < dR_nextmin ) {
            dR_nextmin = dR;
        }
    }
    fJetVars.DR = dR_min;
    //fJetVars.nextDR = dR_nextmin;
}


/* A generic method for obtaining the flavour of partons with the same priority. */
void JetBase::GhostFlavor()
{
    unsigned partons = 0;
    for (auto part : fJetParts) {
        if (part.user_index() > 0) continue; /* Not ghosts */
        unsigned pos = -part.user_index()-1;
        int flav = fGhosts[pos].user_index();
        if (fFlavour==0 || fFlavour==21) {
            fJetVars.partonPT += fGhosts[pos].pt();
            fFlavour = flav;
            ++partons;
        } else if (flav==21) {
            fJetVars.partonPT += fGhosts[pos].pt();
            ++partons;
        } else if (fFlavour+flav==0) {
            fJetVars.partonPT += fGhosts[pos].pt();
            fFlavour = 21;
            ++partons;
        } else {
            partons = 0;
            break;
        }
    }

    if (partons!=0) {
        fJetVars.partonPT /= partons;
        fFlavour = abs(fFlavour);
    } else {
        fJetVars.partonPT = 0.0;
        fFlavour = 0;
    }
}


/* Combine LO and CLO (or other similar purely partonic). UNTESTED */
void JetBase::EqualFlavor() {
    map<unsigned,int> flavs;
    map<unsigned,double> pt;
    map<unsigned,unsigned> partons;
    for (auto part : fJetParts) {
        if (part.user_index() > 0) continue; /* Not ghosts */
        unsigned pos = -part.user_index()-1;
        int flav = fPDGCode[pos];
        unsigned status = fAnalysisStatus[pos];
        if (flavs.find(status)==flavs.end()) {
            flavs[status] = flav;
            pt[status] = part.pt()*pow(10,18);
            partons[status] = 1;
        } else {
            if (flavs[status]==21) {
                flavs[status] = flav;
                pt[status] += part.pt()*pow(10,18);
                ++partons[status];
            } else if (flavs[status]==21) {
                pt[status] += part.pt()*pow(10,18);
                ++partons[status];
            } else if (fFlavour+flav==0) {
                flavs[status] = 21;
                pt[status] += part.pt()*pow(10,18);
                ++partons[status];
            } else {
                flavs[status] = 0;
                break;
            }
        }
    }

    /* Check that there is an agreement of the flavor */
    for (auto flav : flavs) {
        if (fFlavour!=0 && fFlavour!=flav.second) {
            fFlavour=0;
            return;
        } else if (fFlavour==0) {
            fFlavour = flav.second;
        }
    }

    if (flavs.find(8)!=flavs.end()) {
        fJetVars.partonPT = pt[8]/partons[8];
    } else {
        for (auto prts : partons) {
            if (prts.second > 0) {
                fJetVars.partonPT = pt[prts.first]/prts.second;
                return;
            }
        }
        fFlavour = 0;
    }
}

/* Priority for LO partons, then check FS. UNTESTED */
void JetBase::LOFSFlavor()
{
    int flav = -1, algoflav = -1;
    double tmpPt = 0;
    for (auto part : fJetParts) {
        if (part.user_index() > 0) continue; /* Disregard non-ghosts */
        unsigned pos = -part.user_index()-1;
        int pdgCode = fPDGCode[pos];
        int status = fAnalysisStatus[pos];
        if (status==8 || status==7) {
            if (flav != -1) {
                if (pdgCode==flav) {
                    fJetVars.partonPT += part.pt()*pow(10,18);
                    fJetVars.partonPT /= 2.0;
                } else {
                    flav = -2;
                }
            } else {
                flav = pdgCode;
                fJetVars.partonPT = part.pt()*pow(10,18);
            }
        } else {
            if (tmpPt < part.pt()) {
                algoflav = pdgCode;
                tmpPt = part.pt();
            }
        }

    }
    if (flav >= 0) {
        fFlavour = abs(flav);
    } else {
        fFlavour = abs(algoflav);
        if (fFlavour > 6 && fFlavour != 21)
            fFlavour = 0;
        fJetVars.partonPT = tmpPt;
    }
}


void JetBase::HadrFlavor()
{
    fFlavour = 0;
    int partonFlav = 0, hardestLightParton = 0;
    double lightPT = 0, partonPT = 0;

    for (auto part : fJetParts) {
        if (part.user_index() > 0) continue; /* Select ghosts */

        unsigned pos =-part.user_index()-1;
        int pdgCode = abs(fPDGCode[pos]);
        int status = fAnalysisStatus[pos]-2;
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


void JetBase::AlgoFlavor(unsigned jet)
{
    fFlavour = 0;
    int hardestLightParton = 0;
    double lightDR = 0, lightPT = 0;

    for (auto part : fPartons) {
        double dR = fSortedJets[jet].delta_R( part );
        int status = fAnalysisStatus[part.user_index()];
        int pdgCode = abs(fPDGCode[part.user_index()]);

        if (dR > 0.3) continue;

        double PT = part.pt();
        if (status==4 || status==5) {
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


void JetBase::FSFlavor()
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
