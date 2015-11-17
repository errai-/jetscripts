// This is a monstrous script for comparing jet data from pythia6, pythia8 and Herwig++

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "THStack.h"
#include "TH1D.h"
#include "TPad.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <utility>

#include "tdrstyle_mod14.C"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;

const int ptBins = 45.;//29.;//61.;
const double ptRange[]=
    //{18, 21, 24, 
    {28, 32, 37, 43, 49,
     56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
     3637, 3832, 4037};

inline bool compatibility(double mass_sum, double mass_diff) {
    return (mass_sum < 400 && mass_sum > 300 && mass_diff < 50);
}

inline bool mass_study(double m1, double m2, double n1, double n2, bool noisy, int& id) {
    double sum_1 = m1+n2, sum_2 = m2+n1;
    double diff_1 = fabs(m1-n2), diff_2 = fabs(m2-n1);
    
    unsigned success_count = 0;
    if (compatibility(sum_1,diff_1)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m1 << "Jet t " << n2 << endl;
        id = 0;
    }
    if (compatibility(sum_2,diff_2)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m2 << " Jet t " << n1 << endl;
        id = 1;
    }
    
    if (success_count > 0)
        return true;
    return false;
}


void Fracs(string file) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    TH1D *Wlepton = new TH1D("","m_W lepton",50,60.0,110.0);
    TH1D *Wjet = new TH1D("","m_W jet",50,60.0,110.0);
    TH1D *Wjetc = new TH1D("","m_W jet",50,60.0,110.0);
    TH1D *tlepton = new TH1D("","m_t lepton",100,140.0,240.0);
    TH1D *tjet = new TH1D("","m_t jet",100,140.0,240.0);
    TH1D *bboth = new TH1D("","m_b jet",100,3,23);
    TH1D *rbq = new TH1D("","pTb/pTW",100,0,3);
    
    static const Int_t kMaxfJets = 100;

    Int_t           mJets;
    Double_t        mX[kMaxfJets];   //[mJets]
    Double_t        mY[kMaxfJets];   //[mJets]
    Double_t        mZ[kMaxfJets];   //[mJets]
    Double_t        mT[kMaxfJets];   //[mJets]

    Double_t        mWeight;
    Int_t           mFlav[kMaxfJets];   //[mJets]

    /* Tree setup */
    TTree* jetTree;

    TChain* jetChain = new TChain("JetTree","");
    jetChain->Add(file.c_str()); jetTree = jetChain;

    jetTree->SetMakeClass(1);

    jetTree->SetBranchAddress("fJets", &mJets);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fX", mX);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fY", mY);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", mZ);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fT", mT);
    jetTree->SetBranchAddress("fWeight", &mWeight);
    jetTree->SetBranchAddress("fJets.fFlav", mFlav);

    /* event loop */
    std::size_t mCount = 0;
    std::size_t mN = jetTree->GetEntries();
    int success = 0, nonb1 = 0, nonb0 = 0, fluu = 0;
    for(size_t x=0; x != mN; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>mJets);

        vector<unsigned> flavours;
        vector<TLorentzVector> bjets, ljets;
        TLorentzVector MET, lepton;
        unsigned flav_count = 0;
        for (int i = 0; i < mJets; ++i) {
            TLorentzVector tmpVec(mX[i],mY[i],mZ[i],mT[i]);
            if (mFlav[i]==10)
                MET = tmpVec;
            else if (mFlav[i]==11||mFlav[i]==13||mFlav[i]==15)
                lepton = tmpVec;
            else {
                if (tmpVec.Pt() < 30)
                    continue;

                if (mFlav[i]==5)
                    bjets.push_back(tmpVec);
                else {
                    if (mFlav[i]!=0)
                        ++flav_count;
                    ljets.push_back(tmpVec);
                    flavours.push_back(mFlav[i]);
                }
            }
        }
        if (bjets.size()!=2) {
            if (bjets.size()==1)
                ++nonb1;
            else
                ++nonb0;
            continue;
        }
        if (flav_count!=2) {
            ++fluu;
            continue;
        }

        TLorentzVector t1, t2, t3, t4, t5, t6;
        t1 = MET + lepton;
        if (t1.M() < 60 || t1.M() > 110)
            continue;

        vector<TLorentzVector> working;
        vector< pair<unsigned,unsigned> > working_idx;
        for (auto i = 0u; i < ljets.size()-1; ++i) {
            if (flavours[i]==0) continue;
            for (auto j = i+1; j < ljets.size(); ++j) {
                if (flavours[j]==0) continue;
                t2 = ljets[i] + ljets[j];
                if (t2.M() > 60 && t2.M() < 110) {
                    working.push_back(t2);
                    working_idx.push_back( std::make_pair(i,j) );
                }
            }
        }
        if (working.size() == 0)
            continue;
        
        t3 = t1 + bjets[0];
        t4 = t1 + bjets[1];
        unsigned saccess = 0;
        unsigned best;
        int id;
        for ( unsigned i = 0; i < working.size(); ++i ) {
            t5 = working[i] + bjets[0];
            t6 = working[i] + bjets[1];
            if (mass_study(t3.M(),t4.M(),t5.M(),t6.M(),false,id)) {
                ++saccess;
                t2 = working[i];
                best = i;
            }
        }
        if (saccess == 0)
            continue;
        if (saccess > 1)
            cerr << "    HOX" << endl;
            
        ++success;
        Wlepton->Fill(t1.M());
        Wjet->Fill(t2.M());
        double rbqval = (bjets[0].Pt() + bjets[1].Pt())/t2.Pt();
        Wjetc->Fill(rbqval*t2.M());
        if (id == 0) {
            tlepton->Fill( (t1+bjets[0]).M() );
            tjet->Fill( (t2+bjets[1]).M() );
        } else {
            tlepton->Fill( (t1+bjets[1]).M() );
            tjet->Fill( (t2+bjets[0]).M() );
        }
        rbq->Fill( rbqval );
        bboth->Fill( bjets[0].M() );
        bboth->Fill( bjets[1].M() );
        
//         cout << "Lepton W:" << t1.M() << " Jet W:" << t2.M() << endl;
//         cout << "Light flavours: " << flavours[working_idx[best].first] << " " << flavours[working_idx[best].second] << endl << endl;
    }
    TCanvas *c1 = new TCanvas("c1");
    Wlepton->SetYTitle("events");
    Wlepton->SetXTitle("m (GeV)");
    Wlepton->Draw();
    TCanvas *c2 = new TCanvas("c2");
    Wjet->SetYTitle("events");
    Wjet->SetXTitle("m (GeV)");
    Wjet->Draw();
    TCanvas *c3 = new TCanvas("c3");
    tlepton->SetYTitle("events");
    tlepton->SetXTitle("m (GeV)");
    tlepton->Draw();
    TCanvas *c4 = new TCanvas("c4");
    tjet->SetYTitle("events");
    tjet->SetXTitle("m (GeV)");
    tjet->Draw();
    TCanvas *c5 = new TCanvas("c5");
    bboth->SetYTitle("events");
    bboth->SetXTitle("m (GeV)");
    bboth->Draw();
    TCanvas *c6 = new TCanvas("c6");
    rbq->SetYTitle("events");
    rbq->SetXTitle("rbq");
    rbq->Draw();
    TCanvas *c7 = new TCanvas("c7");
    Wjetc->SetYTitle("events");
    Wjetc->SetXTitle("m (GeV)");
    Wjetc->Draw();
    cout << success << endl;
    cout << nonb1 << " " << nonb0 << " " << fluu << endl;
}
