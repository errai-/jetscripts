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
    if (!mTree) throw runtime_error("Creating a tree failed");
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
} // Pythia8Tree initializer


void Pythia8Tree::EventLoop()
{
    if (!mInitialized) {
        cerr << "Event loop can be performed only after proper initialization" << endl;
        return;
    }

    /* The event loop */
    unsigned evtNo = 0;
    while (evtNo != mNumEvents) {
        if (!mPythia.next()) continue;
        if (!ParticleLoop()) continue;
        mTree->Fill();

        ++evtNo;
        if (evtNo%mTimerStep==0) mTimer.printTime();
    }

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
    mHardProcCount = 0;
    for (unsigned prt = 0; prt!=mEvent.size(); ++prt)
        if (!ProcessParticle(prt))
            return false;

    /* Sanity checks (for some ROOT-related reason the error produces a segfault) */
    if (   (mMode<4 && mHardProcCount !=2)
        || (mMode==2 && mSpecialIndices.size()!=1)
        || (mMode==3 && mSpecialIndices.size()!=2)
        || (mMode==4 && mHardProcCount != 4) )
    {
        mEvent.list();
        throw std::runtime_error("Unexpected hard process structure");
    }

    /* Adding corrected parton momenta */
    for (auto prt : mPartonHistory)
        mPrtclEvent->AddPrtcl(prt.second.Px(),
                              prt.second.Py(),
                              prt.second.Pz(),
                              prt.second.E(),
                              mEvent[prt.first].id(),
                              8);

    return true;
} // ParticleLoop


void Pythia8Tree::ParticleAdd(unsigned prt, int saveStatus)
{
    mPrtclEvent->AddPrtcl(mEvent[prt].px(),
                          mEvent[prt].py(),
                          mEvent[prt].pz(),
                          mEvent[prt].e(),
                          mEvent[prt].id(),
                          saveStatus);
} // ParticleAdd


void Pythia8Tree::GhostHadronAdd(unsigned prt, bool useStrange)
{
    int id = mEvent[prt].idAbs();

    if (HadrFuncs::HasBottom(id) && !IsExcitedHadronState(prt,5))
        ParticleAdd(prt,7);
    else if (HadrFuncs::HasCharm(id) && !IsExcitedHadronState(prt,4))
        ParticleAdd(prt,6);
    else if (useStrange && (HadrFuncs::HasStrange(id) && !IsExcitedHadronState(prt,3)))
        ParticleAdd(prt,5);
} // GhostHadronAdd


bool Pythia8Tree::IsExcitedHadronState(unsigned prt, int quarkId)
{
    vector<int> daughters = mEvent[prt].daughterList();
    for (int& dtr : daughters)
        if (HadrFuncs::StatusCheck(quarkId, mEvent[dtr].idAbs()))
            return true;
    return false;
} // IsExcitedHadronState


void Pythia8Tree::PropagateHistory(unsigned prt, int hard_count)
{
    if (!mEvent[prt].isParton()) {
        return;
    } else if (mHistory.emplace(prt,hard_count).second) {
        /* Propagate history info recursively */
        for (auto d : mEvent[prt].daughterList())
            PropagateHistory(d,hard_count);
    } else {
        /* Mark repeated cases with a zero flavor */
        if (mHistory[prt]!=hard_count)
            mHistory[prt] = -1;
    }
} // PropagateHistory


inline bool Pythia8Tree::Absent(unsigned int prt)
{
    return std::find(mSpecialIndices.begin(),mSpecialIndices.end(),prt)==mSpecialIndices.end();
}


///////////////////////////////////////////////////
// Event type specific ProcessParticle functions //
///////////////////////////////////////////////////


bool Pythia8Tree::ProcessParticle(unsigned prt)
{
    int id = mEvent[prt].idAbs();

    if (mEvent[prt].isFinalPartonLevel()) {
        /* Accumulate the corrected momentum of FS partons */
        auto it = mHistory.find(prt);
        if (it != mHistory.end() && it->second!=-1)
            mPartonHistory[mHistory[prt]] += TLorentzVector(mEvent[prt].px(),
                                                            mEvent[prt].py(),
                                                            mEvent[prt].pz(),
                                                            mEvent[prt].e());
    }

    if (mEvent[prt].statusAbs()==71 || mEvent[prt].statusAbs()==72) { /* Alt. 61-63 */
        /* Save final parton level */
        ParticleAdd(prt, 4);
        return true;
    } else if (mEvent[prt].isParton()) {
        if (mEvent[prt].statusAbs()==23) { /* Hard process activities */
            /* Initiate corrected hard process parton momentum */
            mPartonHistory.emplace(prt, TLorentzVector());
            /* Propagate history information */
            PropagateHistory(prt, prt);
            ++mHardProcCount;

            /* Save hard process outgoing partons */
            ParticleAdd(prt, 3);
        }
        return true;
    } else if (mEvent[prt].isHadron()) {
        /* Ghost hadrons can be final state particles */
        GhostHadronAdd(prt);
    }

    int cp = CustomProcess(prt);
    if (cp==0)
        return false;
    else if (cp==1)
        return true;

    /* Final-state particles */
    if (mEvent[prt].isFinal() && Absent(prt)) {
        int saveStatus = 1;
        if (mEvent[prt].id()==22 && GammaChecker(prt)) saveStatus = 9;
        ParticleAdd( prt, saveStatus );
    }

    return true;
} // ProcessParticle : Base


inline int P8GenericTree::CustomProcess(unsigned prt)
{
    return 2;
} // CustomProcess : Generic


inline int P8DijetTree::CustomProcess(unsigned prt)
{
    return 2;
} // CustomProcess : Dijet


inline int P8GammajetTree::CustomProcess(unsigned prt)
{
    if (mSpecialIndices.size()==0 && mEvent[prt].statusAbs()==23)
        return (GammaAdd(prt) ? 1 : 0);

    return 2;
} // CustomProcess : Gammajet


bool P8GammajetTree::GammaAdd(unsigned prt)
{
    while (!mEvent[prt].isFinal()) {
        if (mEvent[prt].daughterList().size()>1)
            return false; // Detect pair production
        prt = mEvent[prt].daughter1();
    }
    mSpecialIndices.push_back(prt);
    ParticleAdd(prt, 2);
    ++mHardProcCount;
    return true;
} // GammaAdd


/* This is based on an assumption made about status codes - there are safeguards in muonadd. */
inline int P8ZmumujetTree::CustomProcess(unsigned prt)
{
    /* The first status 62 hits correspond to the hard process Z0 */
    if (mSpecialIndices.size()==0 && mEvent[prt].statusAbs()==62)
        return (MuonAdd(prt) ? 1 : 0);

    return 2;
} // CustomProcess : Zmumujet


bool P8ZmumujetTree::MuonAdd(unsigned prt)
{
    if (mEvent[prt].idAbs()!=23) {
        cerr << "Expected Z, found " << mEvent[prt].name() << endl;
        throw std::logic_error("Z identification malfunction.");
    }

    for (int daughter : mEvent[prt].daughterList()) {
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
        ParticleAdd(mSpecialIndices[i], 2);
    }

    if ( mSpecialIndices.size() == 2 ) {
        ++mHardProcCount;
        return true;
    }

    cerr << "Failed to locate muon pair! " << ++mNumErrs << endl;
    return true;
} // MuonAdd


inline int P8ttbarjetTree::CustomProcess(unsigned prt)
{
    /* Add the outgoing hard process lepton */
    if (mEvent[prt].statusAbs()==23 && mEvent[prt].idAbs() < 20)
        return (LeptonAdd(prt) ? 1 : 0);

    return 2;
} // ProcessParticle : ttbarjet


bool P8ttbarjetTree::LeptonAdd(unsigned prt)
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
    ParticleAdd(prt, 2);
    return true;
} // LeptonAdd


/* Does a photon originate from pions? */
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
