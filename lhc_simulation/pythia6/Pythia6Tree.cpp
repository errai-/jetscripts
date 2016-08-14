#include "Pythia6Tree.h"


Pythia6Tree::Pythia6Tree(Int_t nEvent, string fileName, Int_t nameId, const int mode) :
    mMode(mode),
    mNumEvents(nEvent),
    mInitialized(true),
    mEnergy(13000),
    mTune(0),
    mMPI(true),
    mISR(true),
    mFSR(true),
    mUseStrange(false)
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


void Pythia6Tree::ModeSettings() {
    cerr << "Initializing Pythia6 ..." << endl;
    if (mMode == 1 || mMode == 0) {
        cerr << "  Generating standard QCD events." << endl;
        mPythia->SetMSEL(1);
        // Min and max pthat
        mPythia->SetCKIN(3,20);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 2) {
        cerr << "  Generating photon+jet events." << endl;
        //pythia->SetMSEL(10); //(includes gamma-gamma bkg)
        mPythia->SetMSEL(0);
        mPythia->SetMSUB(14,1);
        mPythia->SetMSUB(29,1);
        mPythia->SetMSUB(115,1);
        // Min and max pthat
        mPythia->SetCKIN(3,20);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 3) {
        cerr << "  Generating Zmumu+jet events." << endl;
        mPythia->SetMSEL(13); // 11 would be the vanilla y*/Z
        // Leave only decay to muons on
        mPythia->SetMDME(174,1,0); // Z decay to d dbar
        mPythia->SetMDME(175,1,0); // Z decay to u ubar
        mPythia->SetMDME(176,1,0); // Z decay to s sbar
        mPythia->SetMDME(177,1,0); // Z decay to c cbar
        mPythia->SetMDME(178,1,0); // Z decay to b bbar
        mPythia->SetMDME(179,1,0); // Z decay to t tbar
        mPythia->SetMDME(182,1,0); // Zee
        mPythia->SetMDME(183,1,0); // Znuenue
        mPythia->SetMDME(184,1,1); // Zmumu
        mPythia->SetMDME(185,1,0); // Znumunumu
        mPythia->SetMDME(186,1,0); // Ztautau
        mPythia->SetMDME(187,1,0); // Znutaunutau
        // Min and max mhat
        mPythia->SetCKIN(1,40);
        mPythia->SetCKIN(2,-1);
        // Min and max pthat
        mPythia->SetCKIN(3,20);
        mPythia->SetCKIN(4,3000);
    } else if (mMode == 4) {
        cerr << "  Generating ttbar lepton+jets events." << endl;
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
    mPythia->SetMSTU(21,1); // Check for errors
    mPythia->SetMSTJ(22,2); // Unstable particle decay:
    mPythia->SetPARJ(71,10); // ctau = 10 mm

    mPythia->SetMSTP(33,0); // no K factors in hard cross sections
    mPythia->SetMSTP(2,1); // which order running alphaS
    mPythia->SetMSTP(51,10042); // Structure function (PDF CTEQ6L1)
    mPythia->SetMSTP(52,2); // LHAPDF

    mPythia->SetMSTP(142,2); // Turn on Pt reweighting

    if (mTune==0) {
        /* Z2*, a classic CMS tune */
        cerr << "  Z2* tune has been selected." << endl;
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
    } else if (mTune==1) {
        /* cuep6s1 tune by CMS, the newest but not used that much */
        cerr << "  CUEP6S1 tune has been selected." << endl;
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

    if (!mISR) {
        mPythia->SetMSTP(61,0);
        cerr << "  WARNING: ISR has been turned off!" << endl;
    }
    if (!mFSR) {
        mPythia->SetMSTP(71,0); // FSR off
        cerr << "  WARNING: FSR has been turned off!" << endl;
    }
    if (!mMPI) {
        mPythia->SetMSTP(81,0); // MPI off
        cerr << "  WARNING: MPI has been turned off!" << endl;
    }

    cerr << "  CMS energy of " << mEnergy/1000 << " TeV has been chosen." << endl;
    mPythia->Initialize("cms", "p", "p", mEnergy);
} // GeneralSettings


// mPythia->Pylist(2);
// mPythia->Pystat(1);
void Pythia6Tree::EventLoop()
{
    /* Simulation loop */
    unsigned evtNo = 0;
    while (evtNo != mNumEvents) {
        mPythia->GenerateEvent();

        //PrintEvent();
//         if (evtNo==63)
//             mPythia->Pylist(2);

        if (!ParticleLoop()) continue;

        mTree->Fill();

        /* Update progress */
        ++evtNo;
        if (evtNo%mTimerStep==0) mTimer.printTime();
    }

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
    mPartonHistory.clear();
    /* Special particle indices are saved to eliminate saving overlap. */
    mSpecialIndices.clear();

    /* Pythia 6 speciality: do the gamma/muon analysis before the particle loop. */
    if (mMode==2)
        GammaAdd();
    else if (mMode==3)
        MuonAdd();
    else if (mMode==4) {
        unsigned partno = 0, leptno = 0;
        for (unsigned i = 13; i <= 16; ++i) {
            unsigned id = abs(mPythia->GetK(i,2));
            ParticleAdd(i,8);

            if (IsParton(i)) {
                ++partno;
                mPartonHistory.emplace(i,TLorentzVector());
            } else if (id < 20) {
                ++leptno;
            }
        }

        if (partno!=2 || leptno!=2 || partno%2!=0 || leptno%2!=0 || (partno+leptno)!=4)
            return false;
    }

    mPrtclEvent->fWeight = 1./mPythia->GetVINT(99);

    /* Particle loop */
    for (unsigned prt = 1; prt <= mPythia->GetN(); ++prt)
        if (!ProcessParticle(prt))
            return false;

    /* Adding corrected parton momenta */
    for (auto prt : mPartonHistory) {
        mPrtclEvent->AddPrtcl(prt.second.Px(),
                              prt.second.Py(),
                              prt.second.Pz(),
                              prt.second.E(),
                              mPythia->GetK(prt.first,2),
                              7);
    }

    return true;
} // ParticleLoop


bool Pythia6Tree::ProcessParticle(unsigned prt)
{
    int status = mPythia->GetK(prt,1);
    int id = abs(mPythia->GetK(prt,2));

    // prt == 7,8: outgoing LO particles (hard process)
    if (mMode > 0) {
        unsigned parent = mPythia->GetK(prt,3);
        if (mMode < 4) {
            if (prt==7 || prt==8) {
                if (status!=21)
                    throw std::runtime_error("False functionality in hard process.");

                mPartonHistory.emplace(prt,TLorentzVector());

                ParticleAdd(prt,8);
                return true;
            } else if (parent==7 || parent==8) {
                mPartonHistory[parent] += LastParton(prt);
            }
        } else if (mMode==4) {
            if (parent >= 13 && parent <= 16) {
                if (IsParton(prt))
                    mPartonHistory[parent] += LastParton(prt);
                else
                    LeptonAdd(prt);
            }
        }
    }

    /* Hadronic and algorithmic definitions: add FS partons or interesting hadrons. */
    if (status >= 11 && status <= 20) {
        if (IsParton(prt)) {
            if (IsFSParton(prt))
                ParticleAdd(prt,6);
            return true;
        } else if (id>=100 && id%100>9) {
            GhostHadronAdd(prt);
        }
    }

    /* Stable particles, check that this is no special case. */
    if (status <= 10 && Absent(prt)) {
        int saveStatus = 1;
        if (id==22 && GammaChecker(prt)) /* Gamma from pi0 */
            saveStatus = 2;
        ParticleAdd(prt,saveStatus);
    }

    return true;
}


inline bool Pythia6Tree::Absent(unsigned int prt)
{
    return std::find(mSpecialIndices.begin(),mSpecialIndices.end(),prt)==mSpecialIndices.end();
}


inline bool Pythia6Tree::IsParton(unsigned int prt)
{
    unsigned id = abs(mPythia->GetK(prt,2));
    return (id <= 6 || id == 21);
}


/* This assumes implicitly that the particle prt is a parton */
inline bool Pythia6Tree::IsFSParton(unsigned int prt)
{
    /* No decays in sight */
    if (mPythia->GetK(prt,4)>mPythia->GetK(prt,5))
        return false;

    /* The parton has a parton daughter */
    for (unsigned dtr = mPythia->GetK(prt,4), N = mPythia->GetK(prt,5); dtr <= N; ++dtr) {
        unsigned id = abs(mPythia->GetK(prt,2));
        if (id > 90 && id < 94)
            return false;
    }

    return true;
}


/* A handle for adding particle information */
void Pythia6Tree::ParticleAdd(unsigned prt, int saveStatus)
{
    mPrtclEvent->AddPrtcl(mPythia->GetP(prt,1),
                          mPythia->GetP(prt,2),
                          mPythia->GetP(prt,3),
                          mPythia->GetP(prt,4),
                          mPythia->GetK(prt,2),
                          saveStatus);
} // ParticleAdd


void Pythia6Tree::GhostHadronAdd(unsigned int prt, bool useStrange)
{
    int id = abs(mPythia->GetK(prt,2));

    if (HadrFuncs::HasBottom(id) && !IsExcitedHadronState(prt,5))
        ParticleAdd(prt,5);
    else if (HadrFuncs::HasCharm(id) && !IsExcitedHadronState(prt,4))
        ParticleAdd(prt,4);
    else if (useStrange && HadrFuncs::HasStrange(id) && !IsExcitedHadronState(prt,3))
        ParticleAdd(prt,3);
}


bool Pythia6Tree::IsExcitedHadronState(unsigned int prt, int quarkId)
{
    for (unsigned dtr = mPythia->GetK(prt,4), N = mPythia->GetK(prt,5); dtr <= N; ++dtr)
        if (HadrFuncs::StatusCheck(quarkId, mPythia->GetK(dtr,2)))
            return true;
    return false;
}


bool Pythia6Tree::GammaAdd()
{
    mSpecialIndices.push_back(9);

    while (mPythia->GetK(mSpecialIndices[0],1)>10)
        mSpecialIndices[0] = mPythia->GetK(mSpecialIndices[0],4);

    if (mPythia->GetK(mSpecialIndices[0],2) != 22)
        throw std::logic_error("The photon finder did not find a photon.");

    ParticleAdd(mSpecialIndices[0],7);
}


bool Pythia6Tree::MuonAdd()
{
    mSpecialIndices.push_back(12); mSpecialIndices.push_back(13);
    while (abs(mPythia->GetK(mSpecialIndices[0],2))!=13) {
        ++mSpecialIndices[0];
        ++mSpecialIndices[1];
    }

    for (unsigned i = 0; i < mSpecialIndices.size(); ++i) {
        if (abs(mPythia->GetK(mSpecialIndices[i],2))!=13)
            throw std::logic_error("The muon finder did not find a muon.");

        while (mPythia->GetK(mSpecialIndices[i],1)>10) {
            for (unsigned probe = mPythia->GetK(mSpecialIndices[i],4);
                 probe <= mPythia->GetK(mSpecialIndices[i],5); ++probe) {
                if (abs(mPythia->GetK(probe,2))==13) {
                    mSpecialIndices[i] = probe;
                    break;
                }
            }
        }
        ParticleAdd(mSpecialIndices[i],7);
    }
    return true;
}


bool Pythia6Tree::LeptonAdd(unsigned prt)
{
    /* Charged lepton input: find a final-state charged lepton. */
    int type = abs(mPythia->GetK(prt,2))%2;
    if (type) {
        /* Charged leptons require processing. */
        while (mPythia->GetK(prt,1)>9) {
            unsigned tmpPrt = 0;
            for (unsigned i = mPythia->GetK(prt,4); i <= mPythia->GetK(prt,5); ++i) {
                unsigned dType = abs(mPythia->GetK(i,2))%2;
                unsigned dId = abs(mPythia->GetK(i,2));
                if ((dId>10 && dId < 20) && dType==1)
                    tmpPrt = i;
            }

            /* This occurs around 25-30% of the time originating from tau decay */
            if (tmpPrt == 0)
                return false; /* Charged lepton decays to partons and a neutrino */

            prt = tmpPrt;
        }

        if (abs(mPythia->GetK(prt,2))==15)
            throw std::logic_error("Final-state tau spotted.");
    }
    mSpecialIndices.push_back(prt);
    ParticleAdd(prt, 7);
    return true;
}


bool Pythia6Tree::GammaChecker(unsigned prt)
{
    /* One mother, which is pi0? */
    int mother = mPythia->GetK(prt,3);
    if (!mother || abs(mPythia->GetK(mother,2))!=111) return false;

    int d1 = mPythia->GetK(mother,4);
    int d2 = mPythia->GetK(mother,5);
    double eDifference = mPythia->GetP(mother,4);
    for (int daugh = mPythia->GetK(mother,4), N = mPythia->GetK(mother,5); daugh <= N; ++daugh)
        eDifference -= mPythia->GetP(daugh,4);

    if (fabs(eDifference) > 0.001)
        return false;

    return true;
}


/* Corrected parton momentum */
TLorentzVector Pythia6Tree::LastParton(unsigned prt)
{
    if (!IsParton(prt) || IsFSParton(prt)) {
        TLorentzVector handle(mPythia->GetP(prt,1),mPythia->GetP(prt,2),mPythia->GetP(prt,3),mPythia->GetP(prt,4));
        return handle;
    }

    TLorentzVector cumulator;
    for (unsigned i = mPythia->GetK(prt,4); i <= mPythia->GetK(prt,5); ++i)
        cumulator += LastParton(i);

    return cumulator;
}


/* Graphviz printing for a single particle */
void Pythia6Tree::PrintParticle(unsigned prt)
{
    unsigned stat = mPythia->GetK(prt,2);
    unsigned mother = mPythia->GetK(prt,3);
    if (mother == 0 && prt > 2)
        mother = prt-2;
    cout << mother << " -> " << prt << " [label=\"" << prt << " " << mPythia->GetK(prt,1) << " " << mPythia->GetK(prt,2) << "\\n";
    cout << mPythia->GetK(prt,3) << " " << mPythia->GetK(prt,4) << " " << mPythia->GetK(prt,5) << "\\n";
    cout << mPythia->GetP(prt,4) << " " << mPythia->GetP(prt,5) << "\\n";
    cout << '"';
    cout << ",penwidth=2,color=\"";
    if (prt>2 && prt < 9)
        cout << "blue";
    else
        cout << "violet";
    cout << "\"];" << endl;

    if (prt == 7 || prt == 8) {
        if (prt==7)
            mother += 1;
        if (prt==8)
            mother -= 1;
        cout << mother << " -> " << prt << " [label=\"" << prt << " " << mPythia->GetK(prt,1) << " " << mPythia->GetK(prt,2) << "\\n";
        cout << mPythia->GetK(prt,3) << " " << mPythia->GetK(prt,4) << " " << mPythia->GetK(prt,5) << "\\n";
        cout << mPythia->GetP(prt,4) << " " << mPythia->GetP(prt,5) << "\\n";
        cout << '"';
        cout << ",penwidth=2,color=\"";
        cout << "blue";
        cout << "\"];" << endl;
    }

    mother = prt;
    prt = mPythia->GetK(prt,4);
    if (prt==mPythia->GetK(mother,5) && mPythia->GetK(prt,2)==92) {
        cout << mother << " -> " << prt << " [label=\"" << prt << " " << mPythia->GetK(prt,1) << " " << mPythia->GetK(prt,2) << "\\n";
        cout << mPythia->GetK(prt,3) << " " << mPythia->GetK(prt,4) << " " << mPythia->GetK(prt,5) << "\\n";
        cout << mPythia->GetP(prt,4) << " " << mPythia->GetP(prt,5) << "\\n";
        cout << '"';
        cout << ",penwidth=2,color=\"";
        cout << "yellow";
        cout << "\"];" << endl;
    }
}


/* Graphviz printing for a whole event */
void Pythia6Tree::PrintEvent()
{
    cout << "digraph test {" << endl;
    cout << "randir=LR;" << endl;
    cout << "ranksep=1.5;" << endl;
    cout << "node [width=0.03,height=0.03,shape=point,label=\"\"];" << endl;
    for (unsigned prt = 1; prt <= mPythia->GetN(); ++prt)
        PrintParticle(prt);
    cout << "}" << endl;
}
