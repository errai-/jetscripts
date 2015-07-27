#include "JetAnalysis.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

/////////
// Setup:
/////////

JetAnalysis::JetAnalysis(TTree *tree, const char *outFile1, const char *outFile2, 
                         int mode, int definition ) 
                         : fChain(0), mMode(mode), mDefinition(definition)
{
    assert(tree);
    Init(tree);
    
    jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
    //fastjet::JetDefinition(fastjet::genkt_algorithm, R, power); 
    
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
    
    InitCI();
    InitFP();  
    
    jetsPerEvent = 2;
    if (mode==1) jetsPerEvent = 2;
    if (mode==2||mode==3) jetsPerEvent = 1;
}

JetAnalysis::~JetAnalysis() 
{
    delete fjEvent; fjEvent = 0;
    if (!fChain) return;
    delete fChain->GetCurrentFile();
    if (!fOutTree) return;
    delete fOutTree->GetCurrentFile();
    
    delete jetDef;
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

void JetAnalysis::InitCI(){
  // Create file on which histogram(s) can be saved.
  chargeIndicator.push_back(new TH1D("gluonjet amount","",150,0,150) );
  chargeIndicator.push_back(new TH1D("quarkjet amount","",150,0,150) );
  chargeIndicator.push_back(new TH1D("gluonjet charge","",30,-15,15) );
  chargeIndicator.push_back(new TH1D("quarkjet charge","",30,-15,15) );
  chargeIndicator.push_back(new TH1D("gluonjet wcharge","",200,-1,1) );
  chargeIndicator.push_back(new TH1D("quarkjet wcharge","",200,-1,1) );
  chargeIndicator.push_back(new TH1D("gluonjet w2charge","",200,-0.8,0.8) );
  chargeIndicator.push_back(new TH1D("quarkjet w2charge","",200,-0.8,0.8) );
  chargeIndicator.push_back(new TH1D("gluonjet w","",250,0,1) );
  chargeIndicator.push_back(new TH1D("quarkjet w","",250,0,1) );        
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
    if (!fChain) return -5;
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


//////////////////////////////////
// Loop over events and particles:
//////////////////////////////////

void JetAnalysis::EventLoop() 
{
    bool hardStudy = true;
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    mTimer.setParams(nentries,2000); mTimer.startTiming();  
    
    cout << mMode << endl;
    int good=0,bad=0;
    for (Long64_t jentry=0; jentry!=nentries; ++jentry) {
        mFlavorIndices.clear();
        if (jentry!=0&&jentry%2000==0) mTimer.printTime();
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);
        assert( fPrtcls_ < kMaxfPrtcls );

        ParticlesToJetsorterInput();

        if ( hardStudy ) {
            for ( auto i : mPartonList ) {
                fjEvent->AddJet(hiddenInputs[i].px(),hiddenInputs[i].py(),hiddenInputs[i].pz(),
                                hiddenInputs[i].e(),0,0,0,0,0,0,0,0,0,0,fWeight,
                                abs(hiddenInputs[i].user_index()),0,0,0);
            }
        } else {
            /* Fastjet algorithm */
            fastjet::JetDefinition jotDof(fastjet::genkt_algorithm, R, power); //(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
            fastjet::ClusterSequence clustSeq(fjInputs, jotDof);
            vector< fastjet::PseudoJet > unsorteds = clustSeq.inclusive_jets( 10. );
            sortedJets = sorted_by_pt( unsorteds );

            if (!GoodEvent()) continue;

            JetLoop(jentry);
        }

        fOutTree->Fill();

        fjEvent->Clear();
    }

    fOutFile = fOutTree->GetCurrentFile();
    fOutTree->AutoSave("Overwrite");
    delete fjEvent;  fjEvent = 0;

    WriteResults();

    fOutFile->Close();
}

bool JetAnalysis::GoodEvent()
{
    if (sortedJets.size()==0) return false;

    if (mMode == 1) {
        /** dijet events: for the 2 leading jets: 
          *  -Back-to-back angle of min 2.8 rad (2.5 rad)
          *  -Minimum pT of 30 GeV.
          *  -Max eta of 2.5
          *  -A third jet has at most 30% of the average pt of the leading jets */ 
        if (   (sortedJets.size()<2) )
        {
            return false;
        }
        //if (   ( sortedJets.size()<2 ) 
        //    || ( sortedJets.size()>2 
        //    &&   0.15*fabs(sortedJets[0].pt()+sortedJets[1].pt())<sortedJets[2].pt())
        //    || ( fabs(sortedJets[0].eta())>2.5 || fabs(sortedJets[1].eta())>2.5 )
        //    || ( sortedJets[1].pt()<30 )
        //    || ( fabs(sortedJets[0].delta_phi_to( sortedJets[1] ))<2.8 ))
        //{
        //    return false;
        //}
        return true;
    } else if (mMode == 2) {
        /** gamma-jet events:
          *  -Check for a sufficient resolution.
          *  -Back-to-back angle of min 2.8 rad
          *  -Minimum pT of 30 GeV
          *  -A cut for the subleading jet pT with respect to gamma pT 
          *  -Max eta of 2.5 */
        if (   ( hiddenInputs[mGammaId].delta_R( sortedJets[0] )<R )
            || ( hiddenInputs[mGammaId].pt()<30 || sortedJets[0].pt()<30 )
            || ( fabs(hiddenInputs[mGammaId].delta_phi_to(sortedJets[0])) < 2.8 )
            || ( sortedJets.size()>1 && sortedJets[1].pt()>0.3*hiddenInputs[mGammaId].pt() )
            || ( fabs(sortedJets[0].eta()) > 2.5 || fabs(hiddenInputs[mGammaId].eta()) > 2.5 ) )
        {
            return false;
        }
        return true;
    } else if (mMode == 3) {
        /** Z-jet events:
          *  -Require sufficient resolution between the leading jet and the muon system.
          *  -Back-to-back angle of min 2.8 rad
          *  -The subleading jet has to have a small-enough pT compared to the muons.
          *  -Min. muon pT 20GeV and 10GeV 
          *  -Min. leading jet pT of 30GeV
          *  -Subleading jet with a pT smaller than 30% of the dimuon system 
          *  -The dimuon invariant mass is required to fall in the 70-110 GeV range 
          *  -Max eta of 2.5 */
        if (   ( mMuonList.size()!=2 )
            || ( hiddenInputs[mMuonList[0]].delta_R(sortedJets[0])<R 
            ||   hiddenInputs[mMuonList[1]].delta_R(sortedJets[0])<R ) 
            || ((hiddenInputs[mMuonList[0]].pt()<20 || hiddenInputs[mMuonList[1]].pt()<10) 
            && ( hiddenInputs[mMuonList[1]].pt()<20 || hiddenInputs[mMuonList[0]].pt()<10))
            || ( sortedJets[0].pt()<30 ) )
        {
            return false;
        }

        /* Dimuon system */
        fastjet::PseudoJet tmpVec = hiddenInputs[mMuonList[0]]; 
        tmpVec += hiddenInputs[mMuonList[1]];
        if (   ( sortedJets.size()>1 && sortedJets[1].pt()>0.3*tmpVec.pt() ) 
            || ( fabs(tmpVec.m())<70 || fabs(tmpVec.m())>110 )
            || ( fabs(tmpVec.delta_phi_to( sortedJets[0] )) < 2.5 )
            || ( fabs(sortedJets[0].eta())>2.5 || fabs(tmpVec.eta())>2.5 ) )
        {
            return false;
        }

        return true;
    }
    return false;
}

void JetAnalysis::ParticlesToJetsorterInput()
{
    fjInputs.clear();
    hiddenInputs.clear();
    mPartonList.clear();
    mMuonList.clear();
    int hiddenCount = 0;
    
    for (size_t i = 0; i != fPrtcls_; ++i) {
        fastjet::PseudoJet particleTemp(fX[i],fY[i], fZ[i], fT[i]);
        int stat = fAnalysisStatus[i];
        
        /* Ghost partons, hadrons and normal particles */
        if (stat==1) {
            particleTemp.set_user_index( i ); /* Save particle index */
            fjInputs.push_back( particleTemp );
        } else if (stat==2) {
            /* Gammajet gammas and Zjet muons are excluded from jet clustering */
            if (mMode==2) {
                mGammaId = hiddenCount++;
                hiddenInputs.push_back( particleTemp );
            } else if (mMode==3) {
                mMuonList.push_back(hiddenCount++);
                hiddenInputs.push_back( particleTemp );
            } else if (mMode==0) {
                particleTemp.set_user_index( i );
                fjInputs.push_back( particleTemp );
            }
        } else if (stat==3 && mDefinition==1) {
            mPartonList.push_back(hiddenCount++);
            particleTemp.set_user_index( fPDGCode[i] );
            hiddenInputs.push_back(particleTemp);
        } else if (mDefinition==2 && (stat==4 ||/*stat == 5 ||*/ stat == 6 || stat == 7 /* || stat == 8 */)) {
            /* Ghost partons==4, hadrons: (strange==5), charm==6, bottom==7, (top==8) */
            particleTemp *= pow( 10, -18 );
            particleTemp.set_user_index( -i );
            fjInputs.push_back( particleTemp );
        } else {
            /* Discard unknown status codes */
            continue;
        }
    }
    assert( fjInputs.size() ); /* The input should not be empty */
}

void JetAnalysis::JetLoop(int jentry)
{
    for (size_t i = 0; i < sortedJets.size(); i++) {
        if ( i == jetsPerEvent ) break;

        jetParts = sorted_by_pt(sortedJets[i].constituents());

        if (mDefinition == 1) {
            PhysicsFlavor(i);
        } else if (mDefinition == 2) {
            FlavorLoop(i);
        }

        Cuts();
        int multiplicity = cutJetParts.size();

        ParticleLoop(i); /* Operations on jet particles */

        TypeSort(); /* Get ready for adding the jet */

        HistFill(i);

        fjEvent->AddJet(sortedJets[i].px(),sortedJets[i].py(),sortedJets[i].pz(),
            sortedJets[i].e(),mChf,mNhf,mPhf,mElf,mMuf,mChm,mNhm,mPhm,mElm,mMum,
            fWeight,mFlavour,multiplicity,PTD(),Sigma2());
    }
}

void JetAnalysis::PhysicsFlavor(size_t i) 
{
    mFlavour = 0;
    for ( auto k : mPartonList ) {
        double dR = sortedJets[i].delta_R( hiddenInputs[k] );
        if ( dR < 0.3 ) {
            if (mFlavour!=0) {
                mFlavour = 0;
                break;
            } else {
                mFlavour = abs(hiddenInputs[k].user_index());
            }
        }
    }
    mIsHadron = 0;
    mQuarkJetCharge = ChargeSign(mFlavour);
}

void JetAnalysis::FlavorLoop(size_t i)
{
   /* Particle identification.
    * Determine whether a jet is dominated by quarks or by gluons.
    * Looping stops when a corresponding jet is found.
    * Hadron flavour is used as a dominating feature.
    * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    * for further information. */
    int hadronFlav = 0, partonFlav = 0, hardestLightParton = 0;
    
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
        } else if (stat = 4) {
            if (!hardestLightParton && (id==1 || id==2 || id==3 || id==21)) { 
                hardestLightParton = abs(id);
            }
            if (!partonFlav && (id==4 || id==5)) {
                partonFlav = id;
            } else if ( id==5 ) {
                partonFlav = id;
            }
        }
    }
    if (!partonFlav) partonFlav = hardestLightParton;
    
    /* mFlavour is determined with the domination of hadronFlav. If parton flavour
        * is used separately, partonFlav tells the complete parton flavour. */
    if (hadronFlav != 0) {
        mFlavour = hadronFlav;
        mIsHadron = 1;
    } else {
        mFlavour = hardestLightParton;
        mIsHadron = 0;
    }
    mQuarkJetCharge = ChargeSign(mFlavour);
}

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


void JetAnalysis::ParticleLoop(size_t i){
  
    TLorentzVector zero(0,0,0,0);

    mPiPlus = zero; mPiMinus = zero;  mPi0Gamma = zero; mGamma = zero; 
    mKaPlus = zero; mKaMinus = zero; mKSZero = zero; mKLZero = zero; 
    mProton = zero; mAproton = zero; mNeutron = zero; mAneutron = zero;
    mLambda0 = zero; mSigma = zero; mElec = zero, mMuon = zero;
    mOthers = zero; mEtSum = zero;
    
    mChargSum=0; mChargWSum=0; mChargW2Sum=0; mW2=0;
    for (unsigned int j = 0; j != jetParts.size(); ++j) {
        if ( jetParts[j].user_index() < 0 ) continue;
        TLorentzVector tmpP( jetParts[j].px(), jetParts[j].py(), jetParts[j].pz(), 
            jetParts[j].e() );

        mChargSum +=  1;
        mChargWSum += (1)
            *jetParts[j].perp()/sortedJets[i].perp();
        mChargW2Sum += (1)
            *pow(jetParts[j].perp()/sortedJets[i].perp(),2);
        mW2 += pow(jetParts[j].perp()/sortedJets[i].perp(),2);
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
    mChf = tmpLorentz.Et()/mEtSum.Et();
    mChm = tmpLorentz.M();
    tmpLorentz = zero;
    
    tmpLorentz += mKSZero;
    tmpLorentz += mKLZero;
    tmpLorentz += mNeutron;
    tmpLorentz += mAneutron;
    tmpLorentz += mLambda0;
    mNhf = tmpLorentz.Et()/mEtSum.Et();
    mNhm = tmpLorentz.M();
    tmpLorentz = zero;
    
    tmpLorentz += mPi0Gamma;
    tmpLorentz += mGamma;
    mPhf = tmpLorentz.Et()/mEtSum.Et();
    mPhm = tmpLorentz.M();
    
    mElf = mElec.Et()/mEtSum.Et();
    mElm = mElec.M();
    
    mMuf = mMuon.Et()/mEtSum.Et();
    mMum = mMuon.M();
}

void JetAnalysis::HistFill(int i){
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
    
    if ( abs(sortedJets[i].eta()) < 1.3 
        && ( sortedJets[i].pt() > 80 && sortedJets[i].pt() < 120 ) )
    {
        if ( mFlavour==21 ){
            chargeIndicator[0]->Fill(jetParts.size());
            chargeIndicator[2]->Fill(mChargSum);
            chargeIndicator[4]->Fill(mChargWSum);
            chargeIndicator[6]->Fill(mChargW2Sum);
            chargeIndicator[8]->Fill(mW2);
        } else if ( mFlavour < 9 ){
            chargeIndicator[1]->Fill(jetParts.size());
            chargeIndicator[3]->Fill(mIsHadron ? mChargSum : mQuarkJetCharge*mChargSum);
            chargeIndicator[5]->Fill(mIsHadron ? mChargWSum : mQuarkJetCharge*mChargWSum);
            chargeIndicator[7]->Fill(mIsHadron ? mChargW2Sum : mQuarkJetCharge*mChargW2Sum);
            chargeIndicator[9]->Fill(mW2);
        }
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

void JetAnalysis::Cuts()
{
    cutJetParts.clear();
    vector<fastjet::PseudoJet> tmpParts;
    
    /* Explicit cuts */
    for (size_t q = 0; q != jetParts.size(); ++q) {
        if (jetParts[q].user_index() < 0) continue;
        int id = abs(fPDGCode[ jetParts[q].user_index() ]);
        if (!(jetParts[q].pt()<1 && (id == 22 || (IsHadron(id) && !IsCharged(id)))) ) {
            tmpParts.push_back(jetParts[q]);
        }
    }
    
    /* Implicit cuts */
    for (size_t q = 0; q != tmpParts.size(); ++q) {
        int id = abs(fPDGCode[ abs(tmpParts[q].user_index()) ]);
        if ( !IsHadron(id) || ( (IsCharged(id) && tmpParts[q].pt()>0.3) || 
            (!IsCharged(id) && tmpParts[q].pt()>3) ) )
        {
            cutJetParts.push_back(tmpParts[q]);
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

double JetAnalysis::PTD()
{  
    double square = 0, linear = 0;
    for(size_t q = 0; q != cutJetParts.size(); ++q) {
        square += pow(cutJetParts[q].pt(),2);
        linear += cutJetParts[q].pt();
    }
    return sqrt(square)/linear;
}

double JetAnalysis::Sigma2()
{
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

    for (int i = 0; i != 10; ++i){
        chargeIndicator[i]->Write();
    }
    
    fOutFile2->Close();
}
