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

#include "tdrstyle_mod14.C"

using std::cout;
using std::endl;
using std::string;
using std::vector;

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

inline bool mass_study(double m1, double m2, double n1, double n2, bool noisy) {
    double sum_1 = m1+n2, sum_2 = m2+n1;
    double diff_1 = fabs(m1-n2), diff_2 = fabs(m2-n1);
    
    unsigned success_count = 0;
    if (compatibility(sum_1,diff_1)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m1 << "Jet t " << n2 << endl;
    }
    if (compatibility(sum_2,diff_2)) {
        ++success_count;
        if (noisy)
            cout << "Lepton t " << m2 << " Jet t " << n1 << endl;
    }
    
    if (success_count > 0)
        return true;
    return false;
}


void Fracs(string file) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    
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
    int success = 0;
    for(size_t x=0; x != mN; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>mJets);

        vector<unsigned> flavours;
        vector<TLorentzVector> bjets, ljets;
        TLorentzVector MET, lepton;
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
                    ljets.push_back(tmpVec);
                    flavour.push_back(mFlav[i]);
                }
            }
        }

        if (bjets.size()!=2)
            continue;

        TLorentzVector t1, t2, t3, t4, t5, t6;
        t1 = MET + lepton;
        if (t1.M() < 60 || t1.M() > 110) {
            return false;
        }

        vector<TLorentzVector> working;
        for (auto i = 0u; i < ljets.size()-1; ++i) {
            if (flavours[i]==0) continue;
            for (auto j = i+1; j < ljets.size(); ++j) {
                if (flavours[j]==0) continue;
                t2 = ljets[i] + ljets[j];
                if (t2.M() > 60 && t2.M() < 110)
                    working.push_back(t2);
            }
        }
        if (working.size() == 0)
            return false;
        if (working.size() > 1)
            cout << "     HOX!!!!" << endl;
        t3 = t1 + bjets[0];
        t4 = t1 + bjets[1];
        t5 = t2 + bjets[0];
        t6 = t2 + bjets[1];
        if (!mass_study(t3.M(),t4.M(),t5.M(),t6.M(),true))
            continue
        ++success;
        
        cout << "Lepton W:" << t1.m() << " Jet W:" << t2.m() << endl;
        for (auto i = 0u; i < ljets.size(); ++i) {
            cout << "Light flavour: " << flavours[i] << endl;
        }
    }
}
