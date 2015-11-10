#include "JetTest.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

#define HARDSTUDY (0)

//////////////////////////////////
// Loop over events and particles:
//////////////////////////////////

/* Do jet clustering on an event level */
void JetAnalysis::EventLoop()
{
    if (fChain == 0) return;
    mAccepted = 0;

    Long64_t nentries = fChain->GetEntries();
    mTimer.setParams(nentries,2000); mTimer.startTiming();

    /* Fastjet algorithm */
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best); 
    //(fastjet::genkt_algorithm, R, power);

    for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
        if (jentry!=0&&jentry%2000==0) mTimer.printTime();

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);
        assert( fPrtcls_ < kMaxfPrtcls );

        ParticlesToJetsorterInput();

        /* Jet clustering */
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        vector< fastjet::PseudoJet > unsorteds = clustSeq.inclusive_jets( 10. );
        sortedJets = sorted_by_pt( unsorteds );
        
        /* Abort in some universal cases, calculate parameters for further selection */
        if (!SelectionParams()) continue;

        JetLoop(jentry);
    }

    cerr << "Accepted: " << mAccepted << endl;
    fOutFile->Close();
}

inline bool compatibility(double mass_sum, double mass_diff) {
    return (mass_sum < 400 && mass_sum > 300 && mass_diff < 50);
}

inline bool mass_study(double m1, double m2, double n1, double n2) {
    double sum_1 = m1+n2, sum_2 = m2+n1;
    double diff_1 = fabs(m1-n2), diff_2 = fabs(m2-n1);
    
    unsigned success_count = 0;
    if (compatibility(sum_1,diff_1)) {
        ++success_count;
        cout << "Lepton t " << m1 << "Jet t " << n2 << endl;
    }
    if (compatibility(sum_2,diff_2)) {
        ++success_count;
        cout << "Lepton t " << m2 << " Jet t " << n1 << endl;
    }
    
    if (success_count > 0)
        return true;
    return false;
}

/* Calculate variables for the newly clustered jets */
void JetAnalysis::JetLoop(int jentry)
{
    unsigned bcount = 0;
    int b1 = 0,b2 = 0, thirties = 0;
    vector<unsigned> flavours;
    for (unsigned i = 0; i < sortedJets.size(); ++i) {
        if ( i == jetsPerEvent ) break;

        jetParts = sorted_by_pt(sortedJets[i].constituents());

        /* Check the jet flavour if not a generic case */
        if (mMode!=0) {
            if (mDefinition == 1)
                PhysicsFlavor(i);
            else if (mDefinition == 2)
                HadronicFlavor(i);
            else if (mDefinition == 3)
                AlgorithmicFlavor(i);
            else if (mDefinition == 4)
                PhysClusterFlavor(i);
        }

        flavours.push_back(mFlavour);
        if (mFlavour == 5) {
            ++bcount;
            if (bcount == 1)
                b1 = i;
            else if (bcount == 2)
                b2 = i;
        }
        
        if (sortedJets[i].pt() > 10)
            ++thirties;

        ParticleLoop(i); /* Operations on jet particles */
    }
    
    if (bcount != 2) {
        return;
    }
    
    fastjet::PseudoJet t1, t2, t3, t4, t5, t6;
    t1 = mMET + auxInputs[mLeptonId];
    if (t1.m() < 60 || t1.m() > 110) {
        return;
    }
      
    vector<fastjet::PseudoJet> working;
    for (auto i = 0u; i < thirties-1; ++i) {
        if (i == b1 || i == b2) continue;
        for (auto j = i+1; j < thirties; ++j) {
            if (j == b1 || j == b2) continue;
            t2 = sortedJets[i] + sortedJets[j];
            if (t2.m() > 60 && t2.m() < 110)
                working.push_back(t2);
        }
    }
    if (working.size() == 0)
        return;
    if (working.size() > 1)
        cout << "     HOX!!!!" << endl;
    t3 = t1 + sortedJets[b1];
    t4 = t1 + sortedJets[b2];
    t5 = t2 + sortedJets[b1];
    t6 = t2 + sortedJets[b2];
    if (!mass_study(t3.m(),t4.m(),t5.m(),t6.m()))
        return;
    ++mAccepted;
    cout << "Lepton W:" << t1.m() << " Jet W:" << t2.m() << endl;
    for (auto i = 0u; i < thirties; ++i) {
        if (i == b1 || i == b2) continue;
        cout << "Light flavour: " << flavours[i] << endl;
    }
}

/* Study particle types in the clustered jets */
void JetAnalysis::ParticleLoop(unsigned i)
{
    
    for (unsigned int j = 0; j != jetParts.size(); ++j) {
        if ( jetParts[j].user_index() < 0 ) continue;
    }
}


/* fjInputs for jet clustering, other lists for additional event information */
void JetAnalysis::ParticlesToJetsorterInput()
{
    fjInputs.clear();
    auxInputs.clear();
    mPartonList.clear();
    mLeptonList.clear();
    int auxCount = 0;
    mMET = fastjet::PseudoJet();
    
    for (unsigned i = 0; i != fPrtcls_; ++i) {
        fastjet::PseudoJet particleTemp(fX[i],fY[i], fZ[i], fT[i]);
        int stat = fAnalysisStatus[i];
        int pdgID = abs(fPDGCode[i]);

        /* Ghost partons, hadrons and normal particles */
        if (stat==1) {
            particleTemp.set_user_index( i ); /* Save particle index */
            if (pdgID == 12 || pdgID == 14 || pdgID == 16)
                mMET += particleTemp;
            else
                fjInputs.push_back( particleTemp );
            if (pdgID == 11 || pdgID == 13) {
                mLeptonList.push_back(auxCount++);
                auxInputs.push_back(particleTemp);
            }
        } else if (stat==2) {
            /* Gammajet gammas and Zjet muons are excluded from jet clustering */
            if (mMode==2) {
                mGammaId = auxCount++;
                auxInputs.push_back( particleTemp );
            } else if (mMode==3) {
                mLeptonList.push_back(auxCount++);
                auxInputs.push_back( particleTemp );
            } else if (mMode==4) {
                if (pdgID < 20) {
                    particleTemp.set_user_index( -pdgID );
                    mLeptonId = auxCount++;
                    auxInputs.push_back( particleTemp );
                }
            } else if (mMode==0) {
                particleTemp.set_user_index( i );
                fjInputs.push_back( particleTemp );
            }
        } else if (mDefinition==1 && stat==3) {
            /* Physics definition: the hard process */
            mPartonList.push_back(auxCount++);
            particleTemp.set_user_index( pdgID );
            auxInputs.push_back(particleTemp);
        } else if (mDefinition==2 && (stat==4 ||/*stat == 5 ||*/ stat == 6 || stat == 7 /* || stat == 8 */)) {
            /* Hadronic definition: Ghost partons==4, hadrons: (strange==5), charm==6, bottom==7, (top==8) */
            particleTemp *= pow( 10, -18 );
            particleTemp.set_user_index( -i );
            fjInputs.push_back( particleTemp );
        } else if (mDefinition==3 && (stat==4 ) ) {
            /* Algorithmic definition: see hadronic definition */
            if (pdgID > 6 && pdgID!=21) continue;
            mPartonList.push_back(auxCount++);
            particleTemp.set_user_index( pdgID );
            auxInputs.push_back( particleTemp );
        } else if (mDefinition==4 && stat==3) {
            /* Physics clustering definition: ghost partons from the hard process */
            particleTemp *= pow( 10, -18 );
            particleTemp.set_user_index( -pdgID );
            fjInputs.push_back( particleTemp );
        } else {
            /* Discard unknown status codes */
            continue;
        }
    }
    
    assert( fjInputs.size() ); /* The input should not be empty */
}

/* Event type specific cuts */
bool JetAnalysis::SelectionParams()
{
    if ( sortedJets.size() == 0 ) return false;

    if (mMode == 0) {
        return true;
    } else if (mMode == 1) {
        
        /* Dijet events: require always at least two jets */
        if ( sortedJets.size() < 2 ) return false;
        
        /** Example of a selection that should be done:
          *  -Back-to-back angle of min 2.8 rad (2.5 rad)
          *  -Minimum pT of 30 GeV.
          *  -Max eta of 2.5
          *  -A third jet has at most 30% of the average pt of the leading jets (alpha) */
        mJetVars.Alpha = (sortedJets.size() > 2) ? 2*sortedJets[2].pt()/fabs(sortedJets[0].pt()+sortedJets[1].pt()) : 0;
        mJetVars.DPhi = fabs(sortedJets[0].delta_phi_to( sortedJets[1] ));
        mJetVars.matchPT = 0;

    } else if (mMode == 2) {
        
        /* Gammajet events: require always sufficient resolution and cuts for gamma eta and pt */
        if (   ( auxInputs[mGammaId].delta_R( sortedJets[0] ) < R )
            || ( auxInputs[mGammaId].pt()<30 || fabs(auxInputs[mGammaId].eta())>2.5 ) )
        {
            return false;
        }
        
        /** Example of a selection that should be done
          *  -Back-to-back angle of min 2.8 rad
          *  -Minimum jet pT of 30 GeV
          *  -A cut for the subleading jet pT with respect to gamma pT (alpha)
          *  -Max jet eta of 2.5 */
        mJetVars.Alpha = (sortedJets.size()>1) ? sortedJets[1].pt()/auxInputs[mGammaId].pt() : 0;
        mJetVars.DPhi = fabs(auxInputs[mGammaId].delta_phi_to(sortedJets[0]));
        mJetVars.matchPT = auxInputs[mGammaId].pt();
        
    } else if (mMode == 3) {
        
        /* Zjet events: require always sufficient resolution and cuts for muon pt and eta */
        if (   ( mLeptonList.size()!=2 )
            || ( auxInputs[mLeptonList[0]].delta_R(sortedJets[0])<R
            ||   auxInputs[mLeptonList[1]].delta_R(sortedJets[0])<R )
            || ((auxInputs[mLeptonList[0]].pt()<20 || auxInputs[mLeptonList[1]].pt()<10)
            &&  (auxInputs[mLeptonList[1]].pt()<20 || auxInputs[mLeptonList[0]].pt()<10))
            || ( fabs(auxInputs[mLeptonList[0]].eta())>2.5 
            ||   fabs(auxInputs[mLeptonList[1]].eta())>2.5 ) )
        {
            return false;
        }

        /* Dimuon system: require always Z0 mass cut (70-110 GeV) */
        fastjet::PseudoJet tmpVec = auxInputs[mLeptonList[0]]; 
        tmpVec += auxInputs[mLeptonList[1]];
        if ( fabs(tmpVec.m())<70 || fabs(tmpVec.m())>110 ) return false;
        
        /** Example of a selection that should be done:
          *  -Back-to-back angle of min 2.8 rad
          *  -The subleading jet has smaller than 30% pT compared to the muons. (alpha)
          *  -Min. jet pT of 30GeV
          *  -Max jet eta of 2.5 */
        mJetVars.Alpha = (sortedJets.size() > 1 ) ? sortedJets[1].pt()/tmpVec.pt() : 0;
        mJetVars.DPhi = tmpVec.delta_phi_to( sortedJets[0] );
        mJetVars.matchPT = tmpVec.pt();
        
    } else if (mMode == 4) {
        
        if ( sortedJets.size() < 4 ) {
            return false;
        }
        
        for ( auto i = 0u; i < 4; ++i ) {
            if (sortedJets[i].pt() < 30 || fabs(sortedJets[i].eta()) > 2.4) {
                return false;
            }
        }
        
        bool trouble = false;
        for ( auto i : mLeptonList ) {
            unsigned idx = auxInputs[i].user_index();
            if (auxInputs[i].pt() > 33 && fabs(auxInputs[i].eta() < 2.1)) {
                trouble = true;
                break;
            }
        }

        /* Requirements for the relevant charged lepton */
        if (auxInputs[mLeptonId].pt() < 33 || fabs(auxInputs[mLeptonId].eta() > 2.1)) {
            return false;
        }
        if (trouble) {
            return false;
        }
        
    } else {
        throw std::runtime_error("Mode problematic");
    }

    return true;
}


/* If there are conflicts, save the preferred flavour with a minus sign. 
   Assumes a 2 -> 2 hard process () 
   dR_min: the smallest jet axis - parton distance
   dR_nextmin: the next smallest jet axis - parton distance */
void JetAnalysis::PhysicsFlavor(unsigned i) 
{
    mFlavour = 0;
    mJetVars.partonPT = 0;
    double dR_min = 10, dR_nextmin = 10; 
    for ( auto k : mPartonList ) {
        double dR = sortedJets[i].delta_R( auxInputs[k] );
        
        if ( dR < dR_min ) {
            if (mFlavour!=0) dR_nextmin = dR_min;
            dR_min = dR;
            mJetVars.partonPT = auxInputs[k].pt();
            mFlavour = abs(auxInputs[k].user_index());
        } else if ( dR < dR_nextmin ) {
            dR_nextmin = dR;
        }
    }
    mJetVars.DR = dR_min;
    mJetVars.nextDR = dR_nextmin;
}

/* Particle identification the modern "Hadronic definition".
 * Determine whether a jet is dominated by quarks or by gluons.
 * Looping stops when a corresponding jet is found.
 * Hadron flavour is used as a dominating feature.
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
 * for further information. */
void JetAnalysis::HadronicFlavor(unsigned i)
{
    int hadronFlav = 0, partonFlav = 0, hardestLightParton = 0;

    // Jetparts is sorted from the greatest to the smallest pt
    for ( size_t k = 0; k != jetParts.size(); ++k ){
        if (jetParts[k].user_index() > 0) continue; /* Not ghosts */
        int idx = -jetParts[k].user_index();
        int id = abs(fPDGCode[idx]);
        int stat = fAnalysisStatus[idx];

        // Hadrons, set the id's to correspond to the hadron flavour
        if (stat == 7) {
            hadronFlav = 5;
        } else if (stat == 6 && hadronFlav!=5) {
            hadronFlav = 4;
        } else if (stat == 4) {
            if (!hardestLightParton && (id==1 || id==2 || id==3 || id==21))
                hardestLightParton = abs(id);
            partonFlav = (id==4 || id==5) && id > partonFlav ? id : partonFlav;
        }
    }
    if (!partonFlav) partonFlav = hardestLightParton;

    /* mFlavour is determined with the domination of hadronFlav. If parton flavour
        * is used separately, partonFlav tells the complete parton flavour. */
    if (hadronFlav != 0) {
        mFlavour = hadronFlav;
    } else {
        mFlavour = partonFlav;
    }
}

/* Algorithmic flavor tagging is somewhat similar to hadronic tagging.
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
 * for further information. */
void JetAnalysis::AlgorithmicFlavor(unsigned i)
{
    int hardestLightParton = 0;
    double lightDR = 0, lightPT = 0;
    mFlavour = 0; mJetVars.partonPT = 0;
    
    for ( auto k : mPartonList ) {
        double dR = sortedJets[i].delta_R( auxInputs[k] );
        int status = auxInputs[k].user_index();

        if (dR > 0.3) continue;
            
        if (status == 4 || status == 5) {
            if (status > mFlavour) {
                mFlavour = status;
                mJetVars.DR = dR;
                mJetVars.partonPT = auxInputs[k].pt();
            }
        } else if (auxInputs[k].pt() > mJetVars.partonPT) {
            hardestLightParton = status;
            lightDR = dR;
            lightPT = auxInputs[k].pt();
        }
    }
    
    
    if (!mFlavour) {
        mFlavour = hardestLightParton;
        mJetVars.DR = lightDR;
        mJetVars.partonPT = lightPT;
    }
}

/* A variant of the physics definition:
 * use ghost particle clustering for the hard process partons. */
void JetAnalysis::PhysClusterFlavor(unsigned i)
{
    mFlavour = 0;
    
    for ( auto part : jetParts ){
        if (part.user_index() > 0) continue; /* Not ghosts */
        
        /* If there are more than one hard process parton within a jet, mark no flavour */
        if (mFlavour != 0) {
            mFlavour = 0;
            mJetVars.partonPT = 0;
            break;
        }
        mFlavour = -part.user_index();
        mJetVars.partonPT = part.pt()*pow(10,18);
    }
    
    // DR is used mainly for cuts, so set this to zero to avoid unwanted pruning
    mJetVars.DR = 0;
}



/////////
// Setup:
/////////

JetAnalysis::JetAnalysis(TTree *tree, const char *outFile1, const char *outFile2, 
                         int mode, int definition ) 
                         : fChain(0), mMode(mode), mDefinition(definition)
{
    assert(tree);
    Init(tree);

    /* Tree: autosave every 1 Gb, 10 Mb cache */ 
    fOutFile = new TFile(outFile1, "RECREATE");
    fOutFile->SetCompressionLevel(1);

    jetsPerEvent = 2;
    if (mode==0) jetsPerEvent = 100;
    if (mode==1) jetsPerEvent = 2;
    if (mode==2||mode==3) jetsPerEvent = 1;
    if (mode==4) jetsPerEvent = 4;
}

JetAnalysis::~JetAnalysis()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}


/* Initializes the tree that is read */
void JetAnalysis::Init(TTree *tree)
{
    /* Set branch addresses and branch pointers */
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("fWeight", &fWeight, &b_fWeight);
    fChain->SetBranchAddress("fPrtcls", &fPrtcls_, &b_fPrtcls_);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fX", fX, &b_fX);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fY", fY, &b_fY);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fZ", fZ, &b_fZ);
    fChain->SetBranchAddress("fPrtcls.fP4.fCoordinates.fT", fT, &b_fT);
    fChain->SetBranchAddress("fPrtcls.fPDGCode", fPDGCode, &b_fPDGCode);
    fChain->SetBranchAddress("fPrtcls.fAnalysisStatus", fAnalysisStatus, &b_fAnalysisStatus);
}


///////////////////////////////////////////////////
// Generic functions, these should not be modified:
///////////////////////////////////////////////////

Int_t JetAnalysis::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}


Long64_t JetAnalysis::LoadTree(Long64_t entry)
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


void JetAnalysis::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}
