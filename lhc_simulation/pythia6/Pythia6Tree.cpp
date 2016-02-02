#include "Pythia6Tree.h"


Pythia6Tree::Pythia6Tree(Int_t nEvent, string fileName, Int_t nameId, const int mode) :
    mMode(mode),
    mNumEvents(nEvent),
    mInitialized(true)
{
    /* Create an instance of the Pythia event generator: */
    mPythia = new TPythia6;
    /* Set a seed value according to the run index and make sure it is used: */
    if ( nameId < 1 || nameId > 10 )
        throw std::runtime_error("Incompatible seed ID");
    mPythia->SetMRPY(1,mSeeds[nameId-1]);
    mPythia->SetMRPY(2,0);
    
    /* Event type: */
    ModeSettings();
    /* Other settings: */
    GeneralSettings();

    /* Try to create a file to write */
    mFile = TFile::Open(fileName.c_str(), "RECREATE");
    if(!mFile->IsOpen()) throw runtime_error("Creating an output file failed");
    mFile->SetCompressionLevel(1);

    /* Output tree: */
    mTree = new TTree("Pythia6Tree", "Pythia6 particle data.");
    if(!mTree) throw runtime_error("Creating a tree failed");
    mTree->SetAutoSave(100000000); /* 0.1 GBytes */
    mTree->SetCacheSize(10000000); /* 100 MBytes */
    TTree::SetBranchStyle(1); /* New branch style */

    /* Connect an event to the tree */
    mPrtclEvent = new PrtclEvent();
    if (!mPrtclEvent) throw runtime_error("Creating an event handle failed");
    mBranch = mTree->Branch("event", &mPrtclEvent, 32000,4);
    if (!mBranch) throw runtime_error("Associating the event handle with the tree failed");
    mBranch->SetAutoDelete(kFALSE);
    mTree->BranchRef();
    
    /* Setup a custom event timer */
    mTimerStep = 1000;
    mTimer.setParams(mNumEvents,mTimerStep);       
    mTimer.startTiming();
} // Pythia6Tree

void Pythia6Tree::GhostHadronAdd(unsigned int prt, bool useStrange)
{
    return;
}

bool Pythia6Tree::IsExcitedHadronState(unsigned int prt, int quarkID)
{
    return false;
}

TLorentzVector Pythia6Tree::LastParton(unsigned int prt)
{
    return TLorentzVector();
}

bool Pythia6Tree::ProcessParticle(unsigned int prt)
{
    return true;
}

void Pythia6Tree::ModeSettings() {
    if (mMode == 1 || mMode == 0) {
        // Standard QCD
        mPythia->SetMSEL(1);
        // Min and max pthat
        mPythia->SetCKIN(3,25);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 2) {
        // photon+jets
        //pythia->SetMSEL(10);
        mPythia->SetMSEL(0);
        mPythia->SetMSUB(14,1);
        mPythia->SetMSUB(29,1);
        mPythia->SetMSUB(115,1);
        // Min and max pthat
        mPythia->SetCKIN(3,10);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 3) {
        // Z+jets
        mPythia->SetMSEL(13); // 11 would be the vanilla y*/Z
        // Leave only decay to muons on
        mPythia->SetMDME( 174,1,0 ); // Z decay to d dbar
        mPythia->SetMDME( 175,1,0 ); // Z decay to u ubar
        mPythia->SetMDME( 176,1,0 ); // Z decay to s sbar
        mPythia->SetMDME( 177,1,0 ); // Z decay to c cbar
        mPythia->SetMDME( 178,1,0 ); // Z decay to b bbar
        mPythia->SetMDME( 179,1,0 ); // Z decay to t tbar
        mPythia->SetMDME( 182,1,0 ); // Zee
        mPythia->SetMDME( 183,1,0 ); // Znuenue
        mPythia->SetMDME( 184,1,1 ); // Zmumu
        mPythia->SetMDME( 185,1,0 ); // Znumunumu
        mPythia->SetMDME( 186,1,0 ); // Ztautau
        mPythia->SetMDME( 187,1,0 ); // Znutaunutau
        // Min and max mhat
        mPythia->SetCKIN(1,40);
        mPythia->SetCKIN(2,-1);
        // Min and max pthat
        mPythia->SetCKIN(3,15);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 4) {
        // ttbar events
        mPythia->SetMSEL(6); // choose top quark
        mPythia->SetMSUB(81,1); // qqbar -> qqbar
        mPythia->SetMSUB(82,1); // gg->qqbar
        mPythia->SetPMAS(6,1,172);
        mPythia->SetCKIN(3,25);
        mPythia->SetCKIN(4,3000);
    } else {
        throw std::runtime_error("The selected mode is nonsense");
    }
} // ModeSettings


void Pythia6Tree::GeneralSettings() {
    int tune = 0;
    
    mPythia->SetMSTU(21,1); // Check for errors
    mPythia->SetMSTJ(22,2); // Unstable particle decay:
    mPythia->SetPARJ(71,10); // ctau = 10 mm
    
    mPythia->SetMSTP(33,0); // no K factors in hard cross sections
    mPythia->SetMSTP(2,1); // which order running alphaS
    mPythia->SetMSTP(51,10042); // Structure function (PDF CTEQ6L1)
    mPythia->SetMSTP(52,2); // LHAPDF

    mPythia->SetMSTP(142,2); // Turn on Pt reweighting

    if (tune==0) {
        /* Z2*, a classic CMS tune */
        mPythia->SetPARP(82,1.921); // pt cutoff, multiparton interactions
        mPythia->SetPARP(89,1800.); // sqrts for which parp82 is set 
        mPythia->SetPARP(90,0.227); // MPI: rescaling power

        mPythia->SetMSTP(95,6); // Color reconnection setParams
        mPythia->SetPARP(77,1.016); // CR
        mPythia->SetPARP(78,0.538); // CR

        mPythia->SetPARP(80,0.1); // Prob. colored parton from BBR

        mPythia->SetPARP(83,0.356); // MPI matter distribution
        mPythia->SetPARP(84,0.651); // MPI matter distribution

        mPythia->SetPARP(62,1.025); // ISR cutoff

        mPythia->SetMSTP(91,1); // Gaussian primordial KT
        mPythia->SetMSTP(93,10.0); // Primordial KT-max

        mPythia->SetMSTP(81,21); // MPI
        mPythia->SetMSTP(82,4); // MPI model
    } else if (tune==1) {
        /* cuep6s1 tune by CMS, the newest but not used that much */
        
        mPythia->SetPARP(82,1.9096); // pt cutoff, multiparton interactions
        mPythia->SetPARP(89,1800.); // sqrts for which parp82 is set 
        mPythia->SetPARP(90,0.2479); // MPI: rescaling power

        mPythia->SetMSTP(95,6); // Color reconnection setParams
        mPythia->SetPARP(77,0.6646); // CR
        mPythia->SetPARP(78,0.5454); // CR

        mPythia->SetPARP(80,0.1); // Prob. colored parton from BBR

        mPythia->SetPARP(83,0.8217); // MPI matter distribution
        mPythia->SetPARP(84,0.651); // MPI matter distribution

        mPythia->SetPARP(62,1.025); // ISR cutoff

        mPythia->SetMSTP(91,1); // Gaussian primordial KT
        mPythia->SetMSTP(93,10.0); // Primordial KT-max

        mPythia->SetMSTP(81,21); // MPI
        mPythia->SetMSTP(82,4); // MPI model
        
        mPythia->SetPARJ(1,0.08); // HAD diquark suppression
        mPythia->SetPARJ(2,0.21); // HAD strangeness suppression
        mPythia->SetPARJ(3,0.94); // HAD strange diquark suppression
        mPythia->SetPARJ(4,0.04); // HAD vectior diquark suppression
        mPythia->SetPARJ(11,0.35); // HAD P(vector meson), u and d only
        mPythia->SetPARJ(12,0.35); // HAD P(vector meson) contains
        mPythia->SetPARJ(13,0.54); // HAD P(vector meson), heavy quarks
        mPythia->SetPARJ(21,0.34); // HAD fragmentation pt
        mPythia->SetPARJ(25,0.63); // HAD eta0 suppression
        mPythia->SetPARJ(26,0.12); // HAD eta0 suppression
    }
    
    //mPythia->SetMSTP(61,0); // ISR off
    //mPythia->SetMSTP(71,0); // FSR off
    //mPythia->SetMSTP(81,0); // MPI off

    //mPythia->SetMSTP(111,0); // Hadronization off

    //mPythia->Initialize("cms", "p", "p", 13000);
    mPythia->Initialize("cms", "p", "p", 8000);
} // GeneralSettings


void Pythia6Tree::EventLoop()
{
    /* Simulation loop */
    std::size_t ev = 0;
    while (ev != mNumEvents) { 
        mPythia->GenerateEvent();

        if ( !ParticleLoop() ) continue;
        /* Print event listing */
        //mPythia->Pylist(2);
        
        mTree->Fill();
        
        /* Update progress */
        ++ev;
        if (ev%mTimerStep==0) mTimer.printTime();
    }
    /* Print event statistics */
    //mPythia->Pystat(1);

    mFile = mTree->GetCurrentFile();
    mTree->AutoSave("Overwrite");
    mFile->Close();
    
    if (mCounter != 0)
        cerr << "Non-zero counter value: " << mCounter << endl;
    
    delete mPythia;
    mPythia = 0;
    delete mPrtclEvent;
    mPrtclEvent = 0;  
} // EventLoop


bool Pythia6Tree::ParticleLoop()
{
    mPrtclEvent->Clear();
    
    /* Special particle indices are saved to eliminate saving overlap */
    mSpecialIndices.clear();
    if (mMode==2) {
        GammaAdd();
    } else if (mMode==3) {
        MuonAdd();
    }
    
    mPrtclEvent->fWeight = 1./mPythia->GetVINT(99);
    
    for (Int_t prt = 1; prt <= mPythia->GetN(); ++prt) {
        int status = mPythia->GetK(prt,1);
        int id = mPythia->GetK(prt,2);
        
        // prt == 7,8: outgoing particles in the hardest subprocess
        if (mMode > 0) {
            if (prt==7 || prt==8) {
                if (mMode==2 && id==22) continue;
                if (mMode==3 && id==23) continue;

                if ( status != 21 ) throw std::runtime_error("False functionality in hardest subprocess");
                
                ParticleAdd(prt,3);
                continue;
            }
        }

        /* Special final-state particles have already been added */
        if ( std::count( mSpecialIndices.begin(), mSpecialIndices.end(), prt)>0 ) {
            continue;
        }
        
        /* Hadronic and algorithmic definition */
        if (status >= 11 && status <= 20) {
            if (abs(id) <= 6 || id == 21)
                ParticleAdd(prt,4);
            else if (HadrFuncs::HasCharm(id))
                ParticleAdd(prt,6);
            else if (HadrFuncs::HasBottom(id))
                ParticleAdd(prt,7);
        }
        
        /* Stable particles */
        if (status <= 10) {
            int saveStatus = 1;
            if ( (mMode==0||mMode==1) && id==22 && GammaChecker(prt))
                saveStatus = 2;
            ParticleAdd(prt,saveStatus);
        }
    }

    return true;
}


/* A handle for adding particle information */
void Pythia6Tree::ParticleAdd(unsigned prt, int saveStatus)
{
    // TODO: history flavor?
    int history = 0;
    
    mPrtclEvent->AddPrtcl( mPythia->GetP(prt,1),
                           mPythia->GetP(prt,2),
                           mPythia->GetP(prt,3),
                           mPythia->GetP(prt,4),
                           mPythia->GetK(prt,2),
                           saveStatus,
                           history);
} // ParticleAdd


bool Pythia6Tree::GammaAdd()
{
    mSpecialIndices.push_back(9);
    
    while ( mPythia->GetK(mSpecialIndices[0],1)>10 ) {
        mSpecialIndices[0] = mPythia->GetK(mSpecialIndices[0],4);
    }
    
    if ( mPythia->GetK(mSpecialIndices[0],2) != 22 ) throw std::logic_error("Expected photon is not photon");
    ParticleAdd(mSpecialIndices[0],2);
}

bool Pythia6Tree::MuonAdd() 
{
    mSpecialIndices.push_back(12); mSpecialIndices.push_back(13);
    while (abs(mPythia->GetK(mSpecialIndices[0],2))!=13) {
        ++mSpecialIndices[0]; ++mSpecialIndices[1];
    }
    
    for ( std::size_t i = 0; i < mSpecialIndices.size(); ++i ) {
        if ( abs(mPythia->GetK(mSpecialIndices[i],2))!=13 ) throw std::logic_error("Expected muon is not muon");
        
        while ( mPythia->GetK(mSpecialIndices[i],1)>10 ) {
            for (int probe = mPythia->GetK(mSpecialIndices[i],4); probe <= mPythia->GetK(mSpecialIndices[i],5); ++probe) {
                if ( abs(mPythia->GetK(probe,2)) == 13 ) {
                    mSpecialIndices[i] = probe;
                    break;
                }
            }
        }
        ParticleAdd(mSpecialIndices[i],2);
    }
    return true;
}


bool Pythia6Tree::LeptonAdd(unsigned int prt)
{
    return true;
}


bool Pythia6Tree::GammaChecker(unsigned prt)
{
    /* One mother, which is pi0? */
    int mother = mPythia->GetK(prt,3);
    if ( !mother || abs(mPythia->GetK(mother,2))!=111 ) return false;

    int d1 = mPythia->GetK(mother,4);
    int d2 = mPythia->GetK(mother,5);
    double eDifference = mPythia->GetP(mother,4);
    for (int daugh = mPythia->GetK(mother,4), N = mPythia->GetK(mother,5); daugh <= N; ++daugh)
        eDifference -= mPythia->GetP(daugh,4);

    if ( fabs(eDifference) > 0.001 )
        return false;

    return true;
}

