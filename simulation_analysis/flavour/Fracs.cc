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
const double dr_limit = 0.5;

void Fracs(string file, string writefile) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");

    vector<TProfile> hists;
    hists.push_back( TProfile("","",ptBins,ptRange) );
    hists.push_back( TProfile("","",ptBins,ptRange) );
    hists.push_back( TProfile("","",ptBins,ptRange) );
    hists.push_back( TProfile("","",ptBins,ptRange) );
    hists.push_back( TProfile("","",ptBins,ptRange) );

    static const Int_t kMaxfJets = 100;

    Int_t           mJets;
    Double_t        mX[kMaxfJets];   //[mJets]
    Double_t        mY[kMaxfJets];   //[mJets]
    Double_t        mZ[kMaxfJets];   //[mJets]
    Double_t        mT[kMaxfJets];   //[mJets]
    Double_t        mAlpha[kMaxfJets]; //[mJets]
    Double_t        mDPhi[kMaxfJets]; //[mJets]
    Double_t        mDR[kMaxfJets]; //[mJets]

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
    jetTree->SetBranchAddress("fJets.fAlpha",mAlpha);
    jetTree->SetBranchAddress("fJets.fDPhi",mDPhi);
    jetTree->SetBranchAddress("fJets.fDR",mDR);

    /* event loop */
    std::size_t mCount = 0;
    std::size_t mN = jetTree->GetEntries();
    for(size_t x=0; x != mN; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>mJets);

        double xDiff = 0;
        double yDiff = 0;

        for (int i = 0; i < mJets; ++i) {
            TLorentzVector tmpVec(mX[i]+xDiff,mY[i]+yDiff,mZ[i],mT[i]);
            // Upper limit safety cut, eta
            //if (fabs(tmpVec.Eta())>1.3) continue;
            //if (mAlpha[i] > 0.2) continue;
            //if (mDPhi[i] < 2.8) continue;

            int fl = abs(mFlav[i]);
            //if (mDR[i] > dr_limit) fl = -1;
            if (fl == 6 || fl == 7) fl -= 2;
            if (fl > 5 && fl != 21) continue;

            hists[0].Fill(tmpVec.Pt(), (fl == 1 || fl == 2 || fl == 3)? 1:0, mWeight);
            hists[1].Fill(tmpVec.Pt(), fl == 21 ? 1:0, mWeight);
            hists[2].Fill(tmpVec.Pt(), fl == 4 ? 1:0, mWeight);
            hists[3].Fill(tmpVec.Pt(), fl == 5 ? 1:0, mWeight);
            hists[4].Fill(tmpVec.Pt(), fl <= 0 ? 1:0, mWeight);
        }
    }

    TFile* writerHandle = new TFile(writefile.c_str(),"RECREATE");
    vector<TH1D*> dummies;
    dummies.push_back( hists[0].ProjectionX("Quark","") );
    dummies.push_back( hists[1].ProjectionX("Gluon","") );
    dummies.push_back( hists[2].ProjectionX("Charm","") );
    dummies.push_back( hists[3].ProjectionX("Bottom","") );
    dummies.push_back( hists[4].ProjectionX("Unmatched","") );
    dummies[0]->SetLineColor(kRed);
    dummies[1]->SetLineColor(kBlack);
    dummies[2]->SetLineColor(kBlue);
    dummies[3]->SetLineColor(kGreen);
    dummies[4]->SetLineColor(kGray);
    for (unsigned i = 0; i < 5; ++i)
        dummies[i]->Write();

    writerHandle->Close();
}
