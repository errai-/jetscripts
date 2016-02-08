#include "Pythia8Tree.h"

bool TTBarSelector::doVetoProcessLevel(Event& process)
{
    unsigned leptons = 0;
    for (unsigned prt = 0; prt!=process.size(); ++prt) {
        if (process[prt].statusAbs() == 23 && process[prt].isLepton())
            ++leptons;
    }
    if (leptons != 2)
        return true;
    return false;
}


Pythia8Tree::Pythia8Tree(string settings, string fileName, int mode) :
    mEvent(mPythia.event),
    mProcess(mPythia.process),
    mInitialized(true)
{
    mMode = mode;
    
    /* Initialization of the Pythia8 run */
    if (mMode == 4) mPythia.setUserHooksPtr(&mTTBarSelector);
    if (!mPythia.readFile(settings.c_str())) throw std::invalid_argument("Error while reading settings"); 
    if (!mPythia.init()) throw runtime_error("Pythia8 initialization failed");
    mPythia.settings.listChanged();

    mNumEvents = mPythia.mode("Main:numberOfEvents");
    mCounter = 0;
    mNumErrs = 0;

    /* Try to create a file to write */
    mFile = new TFile(fileName.c_str(), "RECREATE");
    if(!mFile->IsOpen()) throw runtime_error("Creating an output file failed");
    mFile->SetCompressionLevel(1);

    /* Create a tree. Autosave every 100 Mb, cache of 10 Mb */
    mTree = new TTree("Pythia8Tree","Pythia8 particle data.");
    if(!mTree) throw runtime_error("Creating a tree failed");
    mTree->SetAutoSave(100000000); /* 0.1 GBytes */
    mTree->SetCacheSize(10000000); /* 100 MBytes */
    TTree::SetBranchStyle(1); /* New branch style */

    /* Connect an event to the tree */
    mPrtclEvent = new PrtclEvent();
    if (!mPrtclEvent) throw runtime_error("Creating an event handle failed");
    mPrtclEvent->SetBit(kCanDelete);
    mPrtclEvent->SetBit(kMustCleanup);
    mBranch = mTree->Branch("event", &mPrtclEvent, 32000,4);
    if (!mBranch) throw runtime_error("Associating the event handle with the tree failed");
    mBranch->SetAutoDelete(kFALSE);
    mTree->BranchRef();
    
    /* Setup a custom event timer */
    mTimerStep = 1000;
    mTimer.setParams(mNumEvents,mTimerStep);       
    mTimer.startTiming();
} // Pythia8Tree


void Pythia8Tree::EventLoop()
{
    if (!mInitialized) {
        cerr << "Event loop can be performed only after proper initialization" << endl;
        return;
    }
    
    std::size_t numEvent = 0;
    while (numEvent != mNumEvents) {
        if (!mPythia.next()) continue;
        if (!ParticleLoop()) continue;
        
        /* Print event listing */
        //mEvent.list();

        mTree->Fill();

        /* Update progress */
        ++numEvent;
        if (numEvent%mTimerStep==0) mTimer.printTime();
    }
    /* Print event statistics */
    //mPythia.stat();

    /* Cleaning up: */
    mFile = mTree->GetCurrentFile();
    mTree->AutoSave("Overwrite");
    mFile->Close();
    
    if (mCounter != 0)
        cerr << "Non-zero counter value: " << mCounter << endl;
    
    delete mPrtclEvent; mPrtclEvent = 0;
    mInitialized = false;
} // EventLoop


bool Pythia8Tree::ParticleLoop()
{
    mPrtclEvent->Clear();
    /* Special particle indices are saved to eliminate saving overlap. */
    mSpecialIndices.clear();
    mHistory.clear();
    mPartonHistory.clear();
    
    mPrtclEvent->fWeight = mPythia.info.weight();
    
    /* Particle loop */
    mHardProcCount = 0, mPartonCount = 0;
    for (unsigned prt = 0; prt!=mEvent.size(); ++prt) {
        if (!ProcessParticle(prt))
            return false;
    }
    /* Sanity checks */
    if (   (mMode<4 && mHardProcCount !=2)
        || (mMode==2 && mSpecialIndices.size()!=1)
        || (mMode==3 && mSpecialIndices.size()!=2) 
        || (mMode==4 && mHardProcCount != 4) ) 
    {
        throw std::logic_error("Unexpected hard process structure");
    }
    
    int history_count = 0;
    for (auto prt : mPartonHistory) {
        mPrtclEvent->AddPrtcl(  prt.second.Px(),
                                prt.second.Py(),
                                prt.second.Pz(),
                                prt.second.E(),
                                mEvent[prt.first].id(),
                                8,
                                history_count++);
    }
    
    return true;
} // ParticleLoop


void Pythia8Tree::ParticleAdd(unsigned prt, int saveStatus)
{
    int history = -1;
    auto it = mHistory.find(prt);
    if ( it != mHistory.end() && it->second!=-1 ) {
        history = std::distance(mPartonHistory.begin(),mPartonHistory.find(it->second));
    }
    
    mPrtclEvent->AddPrtcl(  mEvent[prt].px(),
                            mEvent[prt].py(),
                            mEvent[prt].pz(),
                            mEvent[prt].e(),
                            mEvent[prt].id(),
                            saveStatus,
                            history);
} // ParticleAdd


void Pythia8Tree::GhostHadronAdd(unsigned prt, bool useStrange)
{
    int id = mEvent[prt].idAbs();
    unsigned ghostStatus = 0;

    if (HadrFuncs::HasBottom(id) && !IsExcitedHadronState(prt,5)) {
        ghostStatus = 7; /* b Hadrons */
    } else if (HadrFuncs::HasCharm(id) && !IsExcitedHadronState(prt,4)) {
        ghostStatus = 6; /* c Hadrons */
    } else if (useStrange && (HadrFuncs::HasStrange(id) && !IsExcitedHadronState(prt,3))) {
        ghostStatus = 5; /* s Hadrons */
    }

    if (ghostStatus) { ParticleAdd(prt,ghostStatus); }
} // GhostHadronAdd


bool Pythia8Tree::IsExcitedHadronState(unsigned prt, int quarkId) 
{
    vector<int> daughters = mEvent[prt].daughterList();
    for (int& dtr : daughters) {
        if ( HadrFuncs::StatusCheck(quarkId, mEvent[dtr].id()) ) return true;
    }
    return false;
} // IsExcitedHadronState

void Pythia8Tree::PropagateHistory(unsigned int prt, int hard_count)
{
    if ( mHistory.emplace(prt,hard_count).second ) {
        /* Propagate history info recursively */
        for (auto d : mEvent[prt].daughterList())
            PropagateHistory(d,hard_count);
    } else {
        /* Mark repeated cases with a zero flavor */
        if ( mHistory[prt] != hard_count )
            mHistory[prt] = -1;
    }
}

///////////////////////////////////////////////////
// Event type specific ProcessParticle functions //
///////////////////////////////////////////////////


bool Pythia8Tree::ProcessParticle(unsigned prt)
{
    if ( mEvent[prt].isFinalPartonLevel() ) {
        auto it = mHistory.find(prt);
        if (it != mHistory.end() && it->second!=-1)
            mPartonHistory[mHistory[prt]] += TLorentzVector(mEvent[prt].px(),
                                                            mEvent[prt].py(),
                                                            mEvent[prt].pz(),
                                                            mEvent[prt].e());
    }
    
    if ( mEvent[prt].statusAbs() == 71 || mEvent[prt].statusAbs() == 72 ) {
        /* Save final parton level */
        ParticleAdd( prt, 4 );
        return true;
    } else if (mEvent[prt].isParton()) {
        if ( mEvent[prt].statusAbs()==23 ) {
            /* Hard process activities */
            
            /* Initiate corrected hard process parton momentum */
            mPartonHistory.emplace( prt, TLorentzVector() );
            /* Propagate history information */
            PropagateHistory(prt, prt);
            ++mHardProcCount;

            /* Save hard process outgoing partons */
            ParticleAdd( prt, 3 );
        }
        return true;
    } else if (mEvent[prt].isHadron()) {
        /* Ghost hadrons can be final state particles */
        GhostHadronAdd(prt);
    }

    /* Return true only if a parton is added */
    return false;
} // ProcessParticle : base


bool P8GenericTree::ProcessParticle(unsigned int prt)
{
    if (Pythia8Tree::ProcessParticle(prt))
        return true;
    
    /* pi0 photons in a generic event have the status 2 */
    if ( mEvent[prt].isFinal() ) {
        int saveStatus = 1;
        if ( mEvent[prt].id()==22 && GammaChecker(prt) ) saveStatus = 2; 
        ParticleAdd( prt, saveStatus );
    }
    
    return true;
} // ProcessParticle : Generic


bool P8DijetTree::ProcessParticle(unsigned prt)
{
    if (Pythia8Tree::ProcessParticle(prt))
        return true;

    /* Final particles */
    if ( mEvent[prt].isFinal() ) {
        int saveStatus = 1;
        if ( mEvent[prt].id()==22 && GammaChecker(prt) ) saveStatus = 2;
        ParticleAdd( prt, saveStatus );
    }

    return true;
} // ProcessParticle : Dijet


bool P8GammajetTree::ProcessParticle(unsigned prt)
{
    if (Pythia8Tree::ProcessParticle(prt))
        return true;

    /* The first status 62 hits correspond to the hard process gamma */
    if ( mSpecialIndices.size()==0 && mEvent[prt].statusAbs()==62 )
        return GammaAdd(prt);

    /* Final particles */
    if ( mEvent[prt].isFinal() )
        ParticleAdd( prt, 1 );

    return true;
} // ProcessParticle : Gammajet


bool P8GammajetTree::GammaAdd(unsigned prt)
{
    if (mEvent[prt].idAbs()==22) {
        if (!mEvent[prt].isFinal())
            throw std::runtime_error("Unstable gamma!");
        ParticleAdd(prt,2);
        mSpecialIndices.push_back(prt);
        ++mHardProcCount;
        return true;
    } else if (mPartonCount++ < 1) {
        return true;
    }
    
    cerr << "Pair production, no final-state gamma! " << ++mNumErrs << endl;
    return false;
} // GammaAdd


bool P8ZmumujetTree::ProcessParticle(unsigned prt)
{
    if (Pythia8Tree::ProcessParticle(prt))
        return true;

    /* The first status 62 hits correspond to the hard process Z0 */
    if ( mSpecialIndices.size()==0 && mEvent[prt].statusAbs()==62 )
        return MuonAdd(prt);

    /* Skip the muons that have been already added */
    if ( std::count( mSpecialIndices.begin(), mSpecialIndices.end(), prt)>0 )
        return true;

    /* Final particles */
    if ( mEvent[prt].isFinal() )
        ParticleAdd( prt, 1 );
    
    return true;
} // ProcessParticle : Zmumujet


bool P8ZmumujetTree::MuonAdd(unsigned prt)
{
    if (mEvent[prt].idAbs()!=23) {
        std::cerr << "Expected Z, found " << mEvent[prt].name() << endl;
        throw std::logic_error("Z identification malfunction.");
    }
    
    for ( int daughter : mEvent[prt].daughterList() ) {
        if (mEvent[daughter].idAbs()==13) {
            mSpecialIndices.push_back(daughter);
        }
    }
    
    /* Descend to the final muon forms */
    for ( unsigned i = 0; i < mSpecialIndices.size(); ++i ) {
        unsigned counter = 0;
        while (!mEvent[mSpecialIndices[i]].isFinal() && counter++ < 100) {
            vector<int> mus = mEvent[mSpecialIndices[i]].daughterList();
            for (int daughter : mus) {
                if (mEvent[daughter].idAbs()==13) {
                    mSpecialIndices[i] = daughter; break;
                }
            }
        }
        ParticleAdd( mSpecialIndices[i], 2 );
    }
    
    if ( mSpecialIndices.size() == 2 ) { 
        ++mHardProcCount;
        return true;
    }
    
    cerr << "Failed to locate muon pair! " << ++mNumErrs << endl;
    return true;
} // MuonAdd


bool P8ttbarjetTree::ProcessParticle(unsigned prt)
{
    if (Pythia8Tree::ProcessParticle(prt))
        return true;

    /* Add the outgoing hard process lepton */
    if ( mEvent[prt].statusAbs()==23 && mEvent[prt].idAbs() < 20 )
        return LeptonAdd( prt );

    /* Special final-state particles have already been added */
    if ( std::count( mSpecialIndices.begin(), mSpecialIndices.end(), prt)>0 )
        return true;

    /* pi0 photons in a generic event have the status 2 */
    if ( mEvent[prt].isFinal() )
        ParticleAdd( prt, 1 );

    return true;
} // ProcessParticle : ttbarjet


bool P8ttbarjetTree::LeptonAdd(unsigned int prt)
{
    /* Charged lepton input: find a final-state charged lepton
     * neutrino input: add the parent W */
    int type = mEvent[prt].idAbs()%2;
    if (type) {
        /* Charged leptons */
        while (!mEvent[prt].isFinal()) {
            vector<int> leptons = mEvent[prt].daughterList();
            prt = 0;
            for (int daughter : leptons) {
                int dType = mEvent[daughter].idAbs()%2;
                if (mEvent[daughter].isLepton() && dType==1)
                    prt = daughter;
            }
            
            /* This occurs around 25-30% of the time originating from tau decay */
            if (prt == 0)
                return false; /* Charged lepton decays to partons and a neutrino */
        }
        
        if (mEvent[prt].idAbs()==15)
            cerr << "No tau decay, check settings." << endl;
    } else {
        /* Neutrinos - also saved. (Secondary neutrinos ignored.) */
        while (!mEvent[prt].isFinal()) {
            vector<int> leptons = mEvent[prt].daughterList();
            if (leptons.size()>1)
                cerr << "Neutrino decay, check settings." << endl;
            prt = leptons[0];
        }
    }
    mSpecialIndices.push_back(prt);
    ParticleAdd( prt, 2 );
    return true;
} // LeptonAdd


/////////////////////////////////////////////////////////
// Not in current production but kept in for reference //
/////////////////////////////////////////////////////////


/* Used in the generic events - these are not used to anything. */
bool Pythia8Tree::GammaChecker(unsigned prt)
{
    assert( mEvent.size() > prt );

    /* One mother, which is pi0 */
    vector<int> mothers = mEvent[prt].motherList();
    if ( mothers.size()!=1 || abs(mEvent[mothers[0]].id())!=111 ) return false;

    vector<int> daughters = mEvent[mothers[0]].daughterList();
    
    double eDifference = mEvent[mothers[0]].e();
    for ( auto daugh : daughters )
        eDifference -= mEvent[daugh].e();

    if ( fabs(eDifference) > 0.001 )
        return false;

    return true;
}

/* Done within the history propagation */
TLorentzVector Pythia8Tree::LastParton(unsigned prt)
{
    if (mEvent[prt].statusAbs() == 71 || mEvent[prt].statusAbs() == 72) {
        TLorentzVector handle(mEvent[prt].px(),mEvent[prt].py(),mEvent[prt].pz(),mEvent[prt].e());
        return handle;
    }
    
    TLorentzVector cumulator(0,0,0,0);
    for (auto &daughter : mEvent.daughterList(prt) ) {
        cumulator += LastParton(daughter);
    }
    return cumulator;
}


