#include "Pythia8Tree.h"

void Pythia8Tree::EventLoop()
{
    std::size_t numEvent = 0;
    while (numEvent != mNumEvents) {
        if (!mPythia.next()) continue;
        //mEvent.list();

        if ( !ParticleLoop() ) continue; 
            
        /* Add the candidates to the tree */
        for (auto add = mCandidates.begin(), end = mCandidates.begin()+mNextCand; add != end; ++add) {
            mPrtclEvent->AddPrtcl( mEvent[add->first].px(),mEvent[add->first].py(),
                                   mEvent[add->first].pz(),mEvent[add->first].e(),
                                   mEvent[add->first].id(),add->second);
        }
        mTree->Fill();

        /* Update progress */
        ++numEvent;
        if (numEvent%mTimerStep==0) mTimer.printTime();
    }
    //mPythia.stat();

    /* Cleaning up: */
    mFile = mTree->GetCurrentFile();
    mTree->AutoSave("Overwrite");
    mFile->Close();
}

/* A handle for adding particle information */
void Pythia8Tree::ParticleAdd(std::size_t prt, int saveStatus)
{
    while ( mNextCand >= mCandidates.size() ) {
        mCandidates.push_back( std::make_pair(0,0) );
    }
    mCandidates[mNextCand].first = prt;
    mCandidates[mNextCand].second = saveStatus;
    ++mNextCand;
}

TLorentzVector Pythia8Tree::Vogel(unsigned prt)
{
    if (mEvent[prt].statusAbs() == 71 || mEvent[prt].statusAbs() == 72) {
        TLorentzVector handle(mEvent[prt].px(),mEvent[prt].py(),mEvent[prt].pz(),mEvent[prt].e());
        return handle;
    }
    
    TLorentzVector cumulator(0,0,0,0);
    for (auto &daughter : mEvent.daughterList(prt) ) {
        cumulator += Vogel(daughter);
    }
    return cumulator;
}

/* Returns true if event is to be saved */
bool Pythia8Tree::ParticleLoop()
{
    mPrtclEvent->Clear();
    mNextCand = 0;
    /* Special particle indices are saved to eliminate saving overlap. */
    mSpecialIndices.clear();
    
    mPrtclEvent->fWeight = mPythia.info.weight();
    
    /* Particle loop */
    unsigned hardProcCount = 0, partonCount = 0;
    for (unsigned prt = 0; prt!=mEvent.size(); ++prt) {
        
        // TODO: Chase down the hard process partons on horseback, like men once did

        /* Save partons and such */
        if ( mEvent[prt].statusAbs()==23 ) {
            /* Add the outgoing hard process parton (and lepton in ttbar events) */
            if ( mEvent[prt].isParton() ) {
                ParticleAdd( prt, 3 );
                ++hardProcCount;
                //TLorentzVector juu = Vogel(prt);
                //mPrtclEvent->AddPrtcl( juu.X(), juu.Y(), juu.Z(), juu.T(), mEvent[prt].idAbs(), 3 );
                //cout << "1:" << juu.X() << " " << juu.Y() << " " << juu.Z() << " " << juu.T() << endl;
                //cout << "2:" << mEvent[prt].px() << " " << mEvent[prt].py() << " " << mEvent[prt].pz() << " " << mEvent[prt].e() << endl;
            } else if ( mMode==4 && mEvent[prt].idAbs() < 20 ) {
                if (!LeptonAdd( prt )) return false;
            }
        } else if ( mMode>1 && mSpecialIndices.size()==0 && mEvent[prt].statusAbs()==62 ) {
            /* The first status 62 hits correspond to the hard process, use for special particles */
            if (mMode == 2) {
                if (mEvent[prt].idAbs()==22) {
                    if (!GammaAdd( prt ) ) return false;
                    ++hardProcCount;
                } else if (mEvent[prt].isParton() && partonCount++<1) {
                    continue;
                } else {
                    std::cerr << "Pair production, " << mEvent[prt].name() << endl;
                    return false;
                }
            } else if (mMode == 3) {
                if (mEvent[prt].idAbs()==23) {
                    if (!MuonAdd( prt ) ) return false;
                    ++hardProcCount;
                } else {
                    std::cerr << "Expected Z, found " << mEvent[prt].name() << endl;
                    return false;
                }
            }
        } else if (mEvent[prt].isParton() && (mEvent[prt].statusAbs()==71 || mEvent[prt].statusAbs()==72)) {
            /* Adding ghost partons, used by the algorithmic and the modern "ghost" flavour definition */
            ParticleAdd( prt, 4 );
        }
        
        if (mEvent[prt].isHadron())
            GhostHadronAdd(prt);

        /* Special final-state particles have already been added */
        if ( std::count( mSpecialIndices.begin(), mSpecialIndices.end(), prt)>0 ) {
            continue;
        }

        /* pi0 photons in a generic event have the status 2 */
        if ( mEvent[prt].isFinal() ) {
            int saveStatus = 1;
            if ( mMode==0 && mEvent[prt].id()==22 && GammaChecker(prt) ) saveStatus = 2; 
            ParticleAdd( prt, saveStatus );
        }
    }
    
    /* ttbar events: seek lepton+jets (one w to leptons, one to quarks) */
    if ( mMode==4 && mSpecialIndices.size()!=2 ) return false;
    
    /* Sanity checks */
    if (   (mMode<4 && hardProcCount !=2)
        || (mMode==2 && mSpecialIndices.size()!=1)
        || (mMode==3 && mSpecialIndices.size()!=2) 
        || (mMode==4 && hardProcCount != 4) ) 
    {
        throw std::logic_error("Unexpected hard process structure");
    }

    return true;
}


bool Pythia8Tree::GammaAdd(std::size_t prt)
{
    if (mEvent[prt].isFinal()) {
        mSpecialIndices.push_back(prt);
        ParticleAdd(prt,2);
        return true;
    }
    std::cerr << "Unstable gamma found" << endl;
    return false;
}
    
bool Pythia8Tree::MuonAdd(std::size_t prt)
{
    for ( int daughter : mEvent[prt].daughterList() ) {
        if (mEvent[daughter].idAbs()==13) {
            mSpecialIndices.push_back(daughter);
        }
    }
    
    /* Descend to the final muon forms */
    for ( std::size_t i = 0; i < mSpecialIndices.size(); ++i ) {
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
    
    if ( mSpecialIndices.size() == 2 ) return true;
    std::cerr << "Failed to locate muon pair" << endl;
    return true;
}

bool Pythia8Tree::LeptonAdd(std::size_t prt )
{
    /* For a given neutrino expect neutrino and for a charged lepton
     * expect a charged lepton (indicated by type) */
    int type = mEvent[prt].idAbs()%2;
    while (!mEvent[prt].isFinal()) {
        vector<int> leptons = mEvent[prt].daughterList();
        bool stuck = true;
        for (int daughter : leptons) {
            int absId = mEvent[daughter].idAbs();
            if (absId<20 && absId>10 && absId%2==type) { 
                prt = daughter; stuck = false; break;
            }
        }
        /* Check if stuck in a loop (for instance if lepton goes to hadrons) */
        if ( stuck ) return false;
    }
    mSpecialIndices.push_back(prt);
    ParticleAdd( prt, 2 );
    return true;
}

/* See: HadronAndPartonSelector.cc in CMSSW. Indicates whether a ghost hadron 
 * is in an excited state or not. Checks whether a hadron has a daughter of 
 * the same flavour. Parameter quarkId is a PDG quark flavour. */
bool Pythia8Tree::IsExcitedHadronState(std::size_t idx, int quarkId) 
{
    assert( mEvent.size() > idx );
    assert( quarkId>=0 && quarkId<=6 );

    vector<int> daughters = mEvent[idx].daughterList();
    for (int& dtr : daughters) {
        if ( HadrFuncs::StatusCheck(quarkId, mEvent[dtr].id()) ) return true;
    }
    return false;
}

/* Particles needed by the hadronic flavor definition */
void Pythia8Tree::GhostHadronAdd(std::size_t prt, bool useStrange)
{
    int id = mEvent[prt].idAbs();
    unsigned ghostStatus = 0;

    if (HadrFuncs::HasBottom(id) && !IsExcitedHadronState(prt,5)) {
        ghostStatus = 7; /* b Hadrons */
    } else if (HadrFuncs::HasStrange(id) && !IsExcitedHadronState(prt,3)) {
        ghostStatus = 6; /* c Hadrons */
    } else if (useStrange && (HadrFuncs::HasCharm(id) && !IsExcitedHadronState(prt,4))) {
        ghostStatus = 5; /* s Hadrons */
    }

    if (ghostStatus) { ParticleAdd(prt,ghostStatus); }
}

/////////////////////////////////////////////////////////
// Not in current production but kept in for reference //
/////////////////////////////////////////////////////////

/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
bool Pythia8Tree::GammaChecker( std::size_t prt )
{
    assert( mEvent.size() > prt );

    /* One mother, which is pi0 */
    vector<int> mothers = mEvent[prt].motherList();
    if ( mothers.size()!=1 || mEvent[mothers[0]].id()!=111 ) return false;

    vector<int> daughters = mEvent[mothers[0]].daughterList();
    if ( daughters.size()!=2 ) return false;

    double eDifference = fabs( mEvent[mothers[0]].e() - 
        mEvent[daughters[0]].e() - mEvent[daughters[1]].e() );
    if ( eDifference > 0.001 ) return false;

    return true;
}


