#include "JetAnalysis.h"
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

        if ( HARDSTUDY && mDefinition == 1 ) {
            /* Save only the hard process */
            for ( auto i : mPartonList ) {
                fjEvent->AddJet(auxInputs[i].px(),auxInputs[i].py(),auxInputs[i].pz(),
                                auxInputs[i].e(),mJetVars,fWeight,
                                abs(auxInputs[i].user_index()));
            }
        } else {
            /* Jet clustering */
            fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
            vector< fastjet::PseudoJet > unsorteds = clustSeq.inclusive_jets( 10. );
            sortedJets = sorted_by_pt( unsorteds );
            
            /* Abort in some universal cases, calculate parameters for further selection */
            if (!SelectionParams()) continue;

            JetLoop(jentry);
        }

        fOutTree->Fill();
        fjEvent->Clear();
    }

    fOutFile = fOutTree->GetCurrentFile();
    fOutTree->AutoSave("Overwrite");
    delete fjEvent;  fjEvent = 0;
    fOutFile->Close();
}

/* Calculate variables for the newly clustered jets */
void JetAnalysis::JetLoop(int jentry)
{
    for (size_t i = 0; i < sortedJets.size(); ++i) {
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

            Cuts();
            mJetVars.constituents = cutJetParts.size();
            mJetVars.PTD = PTD();
            mJetVars.Sigma2 = Sigma2();
        }

        ParticleLoop(i); /* Operations on jet particles */

        TypeSort(); /* Get ready for adding the jet */

        //HistFill(i); /* Information histograms (redundant) */

        fjEvent->AddJet(sortedJets[i].px(),
                        sortedJets[i].py(),
                        sortedJets[i].pz(),
                        sortedJets[i].e(),
                        mJetVars,fWeight,mFlavour);
    }
}

/* Study particle types in the clustered jets */
void JetAnalysis::ParticleLoop(unsigned i)
{  
    TLorentzVector zero(0,0,0,0);

    mPiPlus = zero; mPiMinus = zero;  mPi0Gamma = zero; mGamma = zero; 
    mKaPlus = zero; mKaMinus = zero; mKSZero = zero; mKLZero = zero; 
    mProton = zero; mAproton = zero; mNeutron = zero; mAneutron = zero;
    mLambda0 = zero; mSigma = zero; mElec = zero, mMuon = zero;
    mOthers = zero; mEtSum = zero;
    
    for (unsigned int j = 0; j != jetParts.size(); ++j) {
        if ( jetParts[j].user_index() < 0 ) continue;
        TLorentzVector tmpP( jetParts[j].px(), jetParts[j].py(), jetParts[j].pz(), 
            jetParts[j].e() );

        mEtSum += tmpP;
        int id = fPDGCode[ jetParts[j].user_index() ];
        if ( id == 211 ) { 
            mPiPlus += tmpP;
        } else if ( id == -211 ) { 
            mPiMinus+= tmpP;
        } else if ( id == 22 ) {
            if ( fAnalysisStatus[ jetParts[j].user_index() ] == 2 ) {
                mPi0Gamma += tmpP;
            } else {
                mGamma += tmpP;
            }
        } else if ( id == 20 ) { // pi0 gamma
            mPi0Gamma += tmpP;
        } else if ( id == 321 ) { 
            mKaPlus += tmpP;
        } else if ( id == -321 ) { 
            mKaMinus += tmpP;
        } else if ( abs( id ) == 310 ) { 
            mKSZero += tmpP;
        } else if ( abs( id ) == 130 ) { 
            mKLZero += tmpP;
        } else if ( id == 2212 ) { 
            mProton += tmpP;
        } else if ( id == -2212 ) { 
            mAproton += tmpP;
        } else if ( id == 2112 ) { 
            mNeutron += tmpP;
        } else if ( id == -2112 ) { 
            mAneutron += tmpP;
        } else if ( abs( id ) == 3122 ) {
            mLambda0 += tmpP;
        } else if ( abs( id ) == 3112 ||
            abs( id ) == 3222 ) {
            mSigma += tmpP;
        } else if ( abs( id ) == 11 ) {
            mElec += tmpP;
        } else if ( abs( id ) == 13 ) {
            mMuon += tmpP;
        } else {
            mOthers += tmpP;        
        }
    }
}

/* Throw the obtained values in temporary containers */
void JetAnalysis::TypeSort()
{
    TLorentzVector zero(0,0,0,0);
    TLorentzVector tmpLorentz = zero;
    
    tmpLorentz += mPiPlus;
    tmpLorentz += mPiMinus;
    tmpLorentz += mKaPlus;
    tmpLorentz += mKaMinus;
    tmpLorentz += mProton;
    tmpLorentz += mAproton;
    tmpLorentz += mSigma;
    mJetVars.chf = tmpLorentz.Et()/mEtSum.Et();
    mJetVars.chm = tmpLorentz.M();
    tmpLorentz = zero;
    
    tmpLorentz += mKSZero;
    tmpLorentz += mKLZero;
    tmpLorentz += mNeutron;
    tmpLorentz += mAneutron;
    tmpLorentz += mLambda0;
    mJetVars.nhf = tmpLorentz.Et()/mEtSum.Et();
    mJetVars.nhm = tmpLorentz.M();
    tmpLorentz = zero;
    
    tmpLorentz += mPi0Gamma;
    tmpLorentz += mGamma;
    mJetVars.phf = tmpLorentz.Et()/mEtSum.Et();
    mJetVars.phm = tmpLorentz.M();
    
    mJetVars.elf = mElec.Et()/mEtSum.Et();
    mJetVars.elm = mElec.M();
    
    mJetVars.muf = mMuon.Et()/mEtSum.Et();
    mJetVars.mum = mMuon.M();
}


//////////////////////////////////////////////
// Particle addition and selection procedures:
//////////////////////////////////////////////


/* fjInputs for jet clustering, other lists for additional event information */
void JetAnalysis::ParticlesToJetsorterInput()
{
    fjInputs.clear();
    auxInputs.clear();
    mPartonList.clear();
    mMuonList.clear();
    int auxCount = 0;

    for (unsigned i = 0; i != fPrtcls_; ++i) {
        fastjet::PseudoJet particleTemp(fX[i],fY[i], fZ[i], fT[i]);
        int stat = fAnalysisStatus[i];
        int pdgID = abs(fPDGCode[i]);

        /* Ghost partons, hadrons and normal particles */
        if (stat==1) {
            particleTemp.set_user_index( i ); /* Save particle index */
            fjInputs.push_back( particleTemp );
        } else if (stat==2) {
            /* Gammajet gammas and Zjet muons are excluded from jet clustering */
            if (mMode==2) {
                mGammaId = auxCount++;
                auxInputs.push_back( particleTemp );
            } else if (mMode==3) {
                mMuonList.push_back(auxCount++);
                auxInputs.push_back( particleTemp );
            } else if (mMode==4) {
                mJetVars.SetZero();
                fjEvent->AddJet(particleTemp.px(),particleTemp.py(),particleTemp.pz(),
                                particleTemp.e(),mJetVars,fWeight,pdgID);
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
            if (mMode==4) {
                mJetVars.SetZero();
                fjEvent->AddJet(particleTemp.px(),particleTemp.py(),particleTemp.pz(),
                                particleTemp.e(),mJetVars,fWeight,pdgID+30);
            }
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
        if (   ( mMuonList.size()!=2 )
            || ( auxInputs[mMuonList[0]].delta_R(sortedJets[0])<R
            ||   auxInputs[mMuonList[1]].delta_R(sortedJets[0])<R )
            || ((auxInputs[mMuonList[0]].pt()<20 || auxInputs[mMuonList[1]].pt()<10)
            &&  (auxInputs[mMuonList[1]].pt()<20 || auxInputs[mMuonList[0]].pt()<10))
            || ( fabs(auxInputs[mMuonList[0]].eta())>2.5 
            ||   fabs(auxInputs[mMuonList[1]].eta())>2.5 ) )
        {
            return false;
        }

        /* Dimuon system: require always Z0 mass cut (70-110 GeV) */
        fastjet::PseudoJet tmpVec = auxInputs[mMuonList[0]]; 
        tmpVec += auxInputs[mMuonList[1]];
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
        cerr << "Nothing to see here" << endl;
    } else {
        throw std::runtime_error("Mode problematic");
    }

    return true;
}


///////////////////////
// Flavour definitions:
///////////////////////


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
    mJetVars.maxDR = dR_nextmin;
    mQuarkJetCharge = ChargeSign(mFlavour);
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
    mQuarkJetCharge = ChargeSign(mFlavour);
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
    
    mQuarkJetCharge = ChargeSign(mFlavour);
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


///////////////
// qgl-studies:
///////////////


double JetAnalysis::PTD()
{
    if (mMode==0) return 0;
    double square = 0, linear = 0;
    for(size_t q = 0; q != cutJetParts.size(); ++q) {
        square += pow(cutJetParts[q].pt(),2);
        linear += cutJetParts[q].pt();
    }
    return sqrt(square)/linear;
}

double JetAnalysis::Sigma2()
{
    if (mMode==0) return 0;
    double weightedDiffs[4] = {0,0,0,0};
    double phi = 0, eta = 0, pT2 = 0;
    
    for(size_t q = 0; q != cutJetParts.size(); ++q) {
        pT2 += pow(cutJetParts[q].pt(),2);
        eta += pow(cutJetParts[q].pt(),2)*cutJetParts[q].eta();
        phi += pow(cutJetParts[q].pt(),2)*cutJetParts[q].phi();
    }
    eta /= pT2; phi /= pT2;

    for(unsigned int q = 0; q != cutJetParts.size(); ++q) 
    {
        double deltaEta = eta-cutJetParts[q].eta();
        double deltaPhi = TVector2::Phi_mpi_pi( phi-cutJetParts[q].phi() );
        double pT2Tmp = pow(cutJetParts[q].pt(),2);
        weightedDiffs[0] += pT2Tmp*deltaEta*deltaEta;
        weightedDiffs[3] += pT2Tmp*deltaPhi*deltaPhi;
        weightedDiffs[1] -= pT2Tmp*deltaEta*deltaPhi;    
    }
    weightedDiffs[2] = weightedDiffs[1];

    TMatrixDSymEigen me( TMatrixDSym(2,weightedDiffs) );
    TVectorD eigenvals = me.GetEigenValues();

    return sqrt(eigenvals[1]/pT2);
}

void JetAnalysis::Cuts()
{
    cutJetParts.clear();
    vector<fastjet::PseudoJet> tmpParts;
    bool cutMode = false;
    
    if (cutMode) {
        /* Explicit cuts */
        for ( auto q : jetParts ) {
            if ( q.user_index() < 0) continue;
            int id = abs(fPDGCode[ q.user_index() ]);
            if (!( q.pt()<1 && (id == 22 || (IsHadron(id) && !IsCharged(id)))) ) {
                tmpParts.push_back( q );
            }
        }
        
        /* Implicit cuts */
        for ( auto q : tmpParts ) {
            int id = abs(fPDGCode[ q.user_index() ]);
            if ( !IsHadron(id) || ( (IsCharged(id) && q.pt()>0.3) || 
                (!IsCharged(id) && q.pt()>3) ) )
            {
                cutJetParts.push_back( q );
            }
        }
    } else {
        for ( auto q : jetParts ) {
            cutJetParts.push_back(q);
        }
    }
}


bool JetAnalysis::IsHadron(int pdg)
{
    if(pdg>99) return true;
    return false;
}

bool JetAnalysis::IsCharged(int pdg)
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


////////////////
// OLD ROUTINES:
////////////////


/* The charge sign of a quark jet */
int JetAnalysis::ChargeSign( int id )
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

void JetAnalysis::HistFill(int i)
{
    FillerHandle( fractionProfilesAll, sortedJets[i].pt(), mEtSum.Et() );
    if (mFlavour==21) {
        gluonQuark->Fill( sortedJets[i].pt(), 1);
        FillerHandle( fractionProfilesGluon, sortedJets[i].pt(), mEtSum.Et() );
    } else if (mFlavour==4 || mFlavour==5) {
        gluonQuark->Fill( sortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesHQuark, sortedJets[i].pt(), mEtSum.Et() );
        FillerHandle( fractionProfilesQuark, sortedJets[i].pt(), mEtSum.Et() );
    } else if (mFlavour==1 || mFlavour==2 || mFlavour==3 ) {
        gluonQuark->Fill( sortedJets[i].pt(), 0);
        FillerHandle( fractionProfilesLQuark, sortedJets[i].pt(), mEtSum.Et() );
        FillerHandle( fractionProfilesQuark, sortedJets[i].pt(), mEtSum.Et() );
    }
}

void JetAnalysis::FillerHandle( vector<TProfile*> &hists, double pt, double scale)
{
    hists[0]->Fill( pt, mPiPlus.Et()   /scale ); 
    hists[1]->Fill( pt, mPiMinus.Et()  /scale );
    hists[2]->Fill( pt, mPi0Gamma.Et() /scale ); 
    hists[3]->Fill( pt, mKaPlus.Et()   /scale );
    hists[4]->Fill( pt, mKaMinus.Et()  /scale ); 
    hists[5]->Fill( pt, mKSZero.Et()   /scale );
    hists[6]->Fill( pt, mKLZero.Et()   /scale ); 
    hists[7]->Fill( pt, mProton.Et()   /scale );
    hists[8]->Fill( pt, mAproton.Et()  /scale ); 
    hists[9]->Fill( pt, mNeutron.Et()  /scale );
    hists[10]->Fill( pt, mAneutron.Et()/scale ); 
    hists[11]->Fill( pt, mGamma.Et()   /scale );
    hists[12]->Fill( pt, mLambda0.Et() /scale ); 
    hists[13]->Fill( pt, mSigma.Et()   /scale );
    hists[14]->Fill( pt, (mElec.Et()+mMuon.Et())/scale ); 
    hists[15]->Fill( pt, mOthers.Et()  /scale );
}

void JetAnalysis::WriteResults()
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
    fOutTree = new TTree("JetTree","Tree with jet data");
    fOutTree->SetAutoSave(1000000000); fOutTree->SetCacheSize(10000000);  
    TTree::SetBranchStyle(1);
    fOutFileName = outFile2;
    fjEvent = new JetEvent;
    fJetBranch = fOutTree->Branch("event", &fjEvent, 32000, 4);
    fJetBranch->SetAutoDelete(kFALSE);
    fOutTree->BranchRef();

    gluonQuark = new TProfile("gq","gq",ptBins,ptRange);

    InitFP();

    jetsPerEvent = 2;
    if (mode==0) jetsPerEvent = 100;
    if (mode==1) jetsPerEvent = 2;
    if (mode==2||mode==3) jetsPerEvent = 1;
    if (mode==4) jetsPerEvent = 4;
}

JetAnalysis::~JetAnalysis()
{
    delete fjEvent; fjEvent = 0;
    if (!fChain) return;
    delete fChain->GetCurrentFile();
    if (!fOutTree) return;
    delete fOutTree->GetCurrentFile();
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

void JetAnalysis::InitFP(){
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
