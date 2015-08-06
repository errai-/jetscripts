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

/* Switches for turning of plots */
const bool ptPlots = true;
const bool uptPlots = false;
const bool etaPlots = true;
const bool debugPlots = false;
const bool indicatorPlots = true;
const bool separatorPlots = false;

const int sampleType=1;
const bool hard=false; // Pythia6 logistics

const int ptBins = 48.;
const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

// Backup cuts for eta and pt
const double ptUpper = 2000;
const double ptLower = 30;

const int etaBins = 40;
const double etaMax = 2.5;//1.3;

const double drRange[]=
    {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.3,2.6,2.9,3.2,3.6,4,4.5,5};
const int drBins = 100;
const double drLower = 0.0;
const double drUpper = 1.0;

const int alphaBins = 100;
const double alphaLower = 0.0;
const double alphaUpper = 0.31;

const int dphiBins = 100;
const double dphiLower = 2.79;
const double dphiUpper = TMath::Pi();

struct HistStore
{
    vector<TProfile> ptFracs;
    vector<TProfile> etaFracs;
    vector<TProfile> uptFracs;
    
    vector<TH1D*> ptProjections;
    vector<TH1D*> etaProjections;
    vector<TH1D*> uptProjections;
    
    THStack* ptStack;
    THStack* etaStack;
    THStack* uptStack;
    
    TH1D* dR;
    TH1D* alpha;
    TH1D* dPhi;
    
    vector<TH1D*> mult;
    vector<TH1D*> ptd;
    vector<TH1D*> s2;
};

bool ReadTree(string file, HistStore& hists, TH1D* pt, TProfile* weighting)
{
    static const Int_t kMaxfJets = 100;

    Int_t           mJets;
    Double_t        mX[kMaxfJets];   //[mJets]
    Double_t        mY[kMaxfJets];   //[mJets]
    Double_t        mZ[kMaxfJets];   //[mJets]
    Double_t        mT[kMaxfJets];   //[mJets]

    Double_t        mWeight;
    Int_t           mFlav[kMaxfJets];   //[mJets]

    Int_t           mConstituents[kMaxfJets];   //[mJets]
    Double_t        mPTD[kMaxfJets];   //[mJets]
    Double_t        mSigma2[kMaxfJets];   //[mJets]

    Double_t        mDR[kMaxfJets];   //[mJets]
    Double_t        mAlpha[kMaxfJets];   //[mJets]
    Double_t        mDPhi[kMaxfJets];   //[mJets]

    /* Tree setup */
    TTree* jetTree;

    TChain* jetChain = new TChain("JetTree","");
    jetChain->Add(file.c_str()); jetTree = jetChain;

    if (!jetTree) return false;
    jetTree->SetMakeClass(1);

    jetTree->SetBranchAddress("fJets", &mJets);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fX", mX);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fY", mY);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", mZ);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fT", mT);
    jetTree->SetBranchAddress("fWeight", &mWeight);
    jetTree->SetBranchAddress("fJets.fFlav", mFlav);
    jetTree->SetBranchAddress("fJets.fConstituents", mConstituents);
    jetTree->SetBranchAddress("fJets.fPTD", mPTD);
    jetTree->SetBranchAddress("fJets.fSigma2", mSigma2);
    jetTree->SetBranchAddress("fJets.fDR",mDR);
    jetTree->SetBranchAddress("fJets.fAlpha",mAlpha);
    jetTree->SetBranchAddress("fJets.fDPhi",mDPhi);

    /* event loop */
    std::size_t mCount = 0;
    std::size_t mN = jetTree->GetEntries();
    for(size_t x=0; x != mN; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>mJets);

        double xDiff = 0;
        double yDiff = 0;
        if (hard) {
            xDiff = -( mX[0]+mX[1] )/2;
            yDiff = -( mY[0]+mY[1] )/2;
        }

        for (int i = 0; i < mJets; ++i) {
            TLorentzVector tmpVec(mX[i]+xDiff,mY[i]+yDiff,mZ[i],mT[i]);
            // Upper limit safety cut, eta
            if (fabs(tmpVec.Eta())>etaMax) continue;

            //mWeight=1;
            hists.ptFracs[0].Fill(tmpVec.Pt(), (mFlav[i] == 21)? 1:0, mWeight);
            hists.ptFracs[1].Fill(tmpVec.Pt(), (mFlav[i] == 1 || mFlav[i] == 2)? 1:0, mWeight);
            hists.ptFracs[2].Fill(tmpVec.Pt(), (mFlav[i] == 3)? 1:0, mWeight);
            hists.ptFracs[3].Fill(tmpVec.Pt(), (mFlav[i] == 4)? 1:0, mWeight);
            hists.ptFracs[4].Fill(tmpVec.Pt(), (mFlav[i] == 5)? 1:0, mWeight);
            hists.ptFracs[5].Fill(tmpVec.Pt(), (mFlav[i] < 0)? 1:0, mWeight);
            
            if (mFlav[i]<0) {
                hists.uptFracs[0].Fill(tmpVec.Pt(), (mFlav[i] == -21)? 1:0, mWeight);
                hists.uptFracs[1].Fill(tmpVec.Pt(), (mFlav[i] == -1 || mFlav[i] == -2)? 1:0, mWeight);
                hists.uptFracs[2].Fill(tmpVec.Pt(), (mFlav[i] == -3)? 1:0, mWeight);
                hists.uptFracs[3].Fill(tmpVec.Pt(), (mFlav[i] == -4)? 1:0, mWeight);
                hists.uptFracs[4].Fill(tmpVec.Pt(), (mFlav[i] == -5)? 1:0, mWeight);
            }
            
            hists.dR->Fill(mDR[i],mWeight);
            hists.alpha->Fill(mAlpha[i],mWeight);
            hists.dPhi->Fill(mDPhi[i],mWeight);

            pt->Fill(tmpVec.Pt(),mWeight);
            weighting->Fill(tmpVec.Pt(),mWeight);

            // Upper and lower limit safety cut, pt
            if(tmpVec.Pt()>ptUpper || tmpVec.Pt()<ptLower) continue;

            ++mCount;
            hists.etaFracs[0].Fill(tmpVec.Eta(), (mFlav[i] == 21)? 1:0, mWeight);
            hists.etaFracs[1].Fill(tmpVec.Eta(), (mFlav[i] == 1 || mFlav[i] == 2)? 1:0, mWeight);
            hists.etaFracs[2].Fill(tmpVec.Eta(), (mFlav[i] == 3)? 1:0, mWeight);
            hists.etaFracs[3].Fill(tmpVec.Eta(), (mFlav[i] == 4)? 1:0, mWeight);
            hists.etaFracs[4].Fill(tmpVec.Eta(), (mFlav[i] == 5)? 1:0, mWeight);
            hists.etaFracs[5].Fill(tmpVec.Eta(), (mFlav[i] < 0)? 1:0, mWeight);

            if(mFlav[i] == 21) {
                hists.mult[0]->Fill(mConstituents[i],mWeight);
                hists.ptd [0]->Fill(mPTD[i],mWeight);
                hists.s2  [0]->Fill(mSigma2[i],mWeight);
            } else if(mFlav[i] == 0) {
                hists.mult[2]->Fill(mConstituents[i],mWeight);
                hists.ptd [2]->Fill(mPTD[i],mWeight);
                hists.s2  [2]->Fill(mSigma2[i],mWeight);
            } else {
                hists.mult[1]->Fill(mConstituents[i],mWeight);
                hists.ptd [1]->Fill(mPTD[i],mWeight);
                hists.s2  [1]->Fill(mSigma2[i],mWeight);
            }
        }
    }
    cout << mCount << " jets" << endl;
    return true;
}

void stacker(HistStore& hists, int p6p8hpp) {
    
    for (std::size_t i = 0; i < 6; ++i) {
        hists.ptProjections.push_back(hists.ptFracs[i].ProjectionX("",""));
        hists.etaProjections.push_back(hists.etaFracs[i].ProjectionX("",""));
        hists.uptProjections.push_back(hists.uptFracs[i].ProjectionX("",""));
    }
    
    hists.ptStack = new THStack("","");
    hists.etaStack = new THStack("","");
    hists.uptStack = new THStack("","");

    if (p6p8hpp==0) {
        tdrDraw(hists.ptProjections[5],"",kOpenStar,kWhite,kSolid,-1,1001,kGray);
        tdrDraw(hists.ptProjections[0],"",kOpenCross,kCyan,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.ptProjections[1],"",kOpenDiamond,kOrange+7,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.ptProjections[2],"",kOpenSquare,kMagenta,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.ptProjections[3],"",kOpenTriangleUp,kSpring,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.ptProjections[4],"",kOpenTriangleDown,kPink,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.etaProjections[5],"",kOpenStar,kWhite,kSolid,-1,1001,kGray);
        tdrDraw(hists.etaProjections[0],"",kOpenCross,kCyan,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.etaProjections[1],"",kOpenDiamond,kOrange+7,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.etaProjections[2],"",kOpenSquare,kMagenta,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.etaProjections[3],"",kOpenTriangleUp,kSpring,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.etaProjections[4],"",kOpenTriangleDown,kPink,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.uptProjections[5],"",kOpenStar,kWhite,kSolid,-1,1001,kGray);
        tdrDraw(hists.uptProjections[0],"",kOpenCross,kCyan,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.uptProjections[1],"",kOpenDiamond,kOrange+7,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.uptProjections[2],"",kOpenSquare,kMagenta,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.uptProjections[3],"",kOpenTriangleUp,kSpring,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.uptProjections[4],"",kOpenTriangleDown,kPink,kSolid,-1,1001,kRed-9);
    } else if (p6p8hpp==1) {
        tdrDraw(hists.ptProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.ptProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.ptProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.ptProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.ptProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.ptProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.etaProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.etaProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.etaProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.etaProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.etaProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.etaProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.uptProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.uptProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.uptProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.uptProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.uptProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.uptProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
    } else if (p6p8hpp==2) {
        tdrDraw(hists.ptProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.ptProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.ptProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.ptProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.ptProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.ptProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.etaProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.etaProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.etaProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.etaProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.etaProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.etaProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
        tdrDraw(hists.uptProjections[5],"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
        tdrDraw(hists.uptProjections[0],"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
        tdrDraw(hists.uptProjections[1],"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
        tdrDraw(hists.uptProjections[2],"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
        tdrDraw(hists.uptProjections[3],"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
        tdrDraw(hists.uptProjections[4],"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
    }
    
    //light_quarks->Add(strange);
    hists.ptStack->Add(hists.ptProjections[4]);
    hists.ptStack->Add(hists.ptProjections[3]);
    hists.ptStack->Add(hists.ptProjections[2]);
    hists.ptStack->Add(hists.ptProjections[1]);
    hists.ptStack->Add(hists.ptProjections[0]);
    hists.ptStack->Add(hists.ptProjections[5]);
    
    //light_quarks->Add(strange);
    hists.etaStack->Add(hists.etaProjections[4]);
    hists.etaStack->Add(hists.etaProjections[3]);
    hists.etaStack->Add(hists.etaProjections[2]);
    hists.etaStack->Add(hists.etaProjections[1]);
    hists.etaStack->Add(hists.etaProjections[0]);
    hists.etaStack->Add(hists.etaProjections[5]);
    
    //light_quarks->Add(strange);
    hists.uptStack->Add(hists.uptProjections[4]);
    hists.uptStack->Add(hists.uptProjections[3]);
    hists.uptStack->Add(hists.uptProjections[2]);
    hists.uptStack->Add(hists.uptProjections[1]);
    hists.uptStack->Add(hists.uptProjections[0]);
    hists.uptStack->Add(hists.uptProjections[5]);
}

void HppP6P8(string herwFile, string p6File, string p8File) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    cout << "Order: herwig++, pythia6, pythia8" << endl;

    /* Herwig */
    HistStore hppHists;
    for (std::size_t i = 0; i < 6; ++i) { 
        hppHists.ptFracs.push_back( TProfile("","",ptBins,ptRange) );
        hppHists.etaFracs.push_back( TProfile("","",etaBins,-etaMax,etaMax) );
        hppHists.uptFracs.push_back( TProfile("","",ptBins,ptRange) );
    }
    
    for (std::size_t i = 0; i < 3; ++i) { 
        hppHists.mult.push_back(new TH1D("","",60,0,60));
        hppHists.ptd.push_back(new TH1D("","",100,0,1));
        hppHists.s2.push_back(new TH1D("","",100,0,0.2));
    }
    
    for (TH1D* mult : hppHists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : hppHists.ptd) { ptd->Sumw2(); }
    for (TH1D* s2 : hppHists.s2) { s2->Sumw2(); }
    
    hppHists.mult[0]->SetLineColor(kRed); hppHists.mult[2]->SetLineColor(kGreen);
    hppHists.ptd[0]->SetLineColor(kRed); hppHists.ptd[2]->SetLineColor(kGreen);
    hppHists.s2[0]->SetLineColor(kRed); hppHists.s2[2]->SetLineColor(kGreen);

    hppHists.dR = new TH1D("",";dR;Events",drBins,drLower,drUpper); 
    hppHists.alpha = new TH1D("",";alpha;Events",alphaBins,alphaLower,alphaUpper); 
    hppHists.dPhi = new TH1D("",";dPhi;Events",dphiBins,dphiLower,dphiUpper);
    
    TH1D* hppPt = new TH1D("hpp Pt",";pT;Events",ptBins,ptRange);
    TProfile* hppWeighting = new TProfile("hppw","",ptBins,ptRange);
    
    if (!ReadTree(herwFile, hppHists, hppPt, hppWeighting)) {
        cout << "Could not load Herwig data" << endl;
        return;
    }

    /* Pythia6 */
    HistStore p6Hists;
    for (std::size_t i = 0; i < 6; ++i) { 
        p6Hists.ptFracs.push_back( TProfile("","",ptBins,ptRange) );
        p6Hists.etaFracs.push_back( TProfile("","",etaBins,-etaMax,etaMax) );
        p6Hists.uptFracs.push_back( TProfile("","",ptBins,ptRange) );
    }
    
    for (std::size_t i = 0; i < 3; ++i) { 
        p6Hists.mult.push_back(new TH1D("","",60,0,60));
        p6Hists.ptd.push_back(new TH1D("","",100,0,1));
        p6Hists.s2.push_back(new TH1D("","",100,0,0.2));
    }
    
    for (TH1D* mult : p6Hists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : p6Hists.ptd) { ptd->Sumw2(); }
    for (TH1D* sigma2 : p6Hists.s2) { sigma2->Sumw2(); }

    p6Hists.mult[0]->SetLineColor(kRed); p6Hists.mult[2]->SetLineColor(kGreen);
    p6Hists.ptd[0]->SetLineColor(kRed); p6Hists.ptd[2]->SetLineColor(kGreen);
    p6Hists.s2[0]->SetLineColor(kRed); p6Hists.s2[2]->SetLineColor(kGreen);

    p6Hists.dR = new TH1D("",";dR;Events",drBins,drLower,drUpper); 
    p6Hists.alpha = new TH1D("",";alpha;Events",alphaBins,alphaLower,alphaUpper); 
    p6Hists.dPhi = new TH1D("",";dPhi;Events",dphiBins,dphiLower,dphiUpper);
    
    TH1D* p6Pt = new TH1D("p6 Pt",";pT;Events",ptBins,ptRange);
    TProfile* p6Weighting = new TProfile("p6w","",ptBins,ptRange);
    
    if (!ReadTree(p6File, p6Hists, p6Pt, p6Weighting)) {
        cout << "Could not load Pythia6 data" << endl;
        return;
    }
    
    /* Pythia8 plot initialization */
    HistStore p8Hists;
    for (std::size_t i = 0; i < 6; ++i) { 
        p8Hists.ptFracs.push_back( TProfile("","",ptBins,ptRange) );
        p8Hists.etaFracs.push_back( TProfile("","",etaBins,-etaMax,etaMax) );
        p8Hists.uptFracs.push_back( TProfile("","",ptBins,ptRange) );
    }
    
    for (std::size_t i = 0; i < 3; ++i) { 
        p8Hists.mult.push_back(new TH1D("","",60,0,60));
        p8Hists.ptd.push_back(new TH1D("","",100,0,1));
        p8Hists.s2.push_back(new TH1D("","",100,0,0.2));
    }
    
    for (TH1D* mult : p8Hists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : p8Hists.ptd) { ptd->Sumw2(); }
    for (TH1D* sigma2 : p8Hists.s2) { sigma2->Sumw2(); }
    
    p8Hists.mult[0]->SetLineColor(kRed); p8Hists.mult[2]->SetLineColor(kGreen);
    p8Hists.ptd[0]->SetLineColor(kRed); p8Hists.ptd[2]->SetLineColor(kGreen);
    p8Hists.s2[0]->SetLineColor(kRed); p8Hists.s2[2]->SetLineColor(kGreen);

    p8Hists.dR = new TH1D("",";dR;Events",drBins,drLower,drUpper); 
    p8Hists.alpha = new TH1D("",";alpha;Events",alphaBins,alphaLower,alphaUpper); 
    p8Hists.dPhi = new TH1D("",";dPhi;Events",dphiBins,dphiLower,dphiUpper);
    
    TH1D* p8Pt = new TH1D("p8 Pt",";pT;Events",ptBins,ptRange);
    TProfile* p8Weighting = new TProfile("p8w","",ptBins,ptRange);
    if (!ReadTree(p8File, p8Hists, p8Pt, p8Weighting)) {
        cout << "Could not load Pythia6 data" << endl;
        return;
    }

    /* Transform data to stacks */
    stacker( p6Hists, 0 );
    stacker( p8Hists, 1 );
    stacker( hppHists, 2 );

    if (ptPlots) {
        TH1D *h1 = new TH1D("h",";p_{T} (GeV);Fraction",ptBins,ptRange);
        setTDRStyle();
        gROOT->ForceStyle();
        gStyle->SetOptStat(kFALSE); //removes old legend
        gStyle->SetAxisColor(1, "XYZ");
        gStyle->SetStripDecimals(kTRUE);
        gStyle->SetTickLength(0.03, "XYZ");
        gStyle->SetNdivisions(510, "XYZ");
        gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
        gStyle->SetPadTickY(1);

        TCanvas *canv = tdrCanvas("c1",h1,12,0,1);
        h1->GetYaxis()->SetTitleOffset(1.25);
        h1->GetXaxis()->SetTitleOffset(1.0);
        h1->GetXaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetLabelSize(0.045);
        h1->GetXaxis()->SetTitleSize(0.045);
        h1->GetYaxis()->SetTitleSize(0.045);
        gPad->SetLogx();
        
        p6Hists.ptStack->SetHistogram(h1);
        p8Hists.ptStack->SetHistogram(h1);
        hppHists.ptStack->SetHistogram(h1);
        gStyle->SetOptStat(kFALSE); //removes old legend
        
        /* Generic stuff: */
        p6Hists.ptStack->SetMaximum(0.95);
        p8Hists.ptStack->SetMaximum(0.95);
        hppHists.ptStack->SetMaximum(0.95);

        /* Fraction plotting */
        h1->GetXaxis()->SetRange(5,48);
        //h1->GetYaxis()->SetRangeUser(-0.001,1.001);
        p8Hists.ptStack->Draw("");
        hppHists.ptStack->Draw("sameP");
        p6Hists.ptStack->Draw("sameP");

        /* Fraction legends */

        //TLegend *leg = tdrLeg(0.2,0.41,0.5,0.91);
        TLegend *leg = tdrLeg(0.3,0.4,0.75,0.85);
        TLegend *heading = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
        TLegend *sample = tdrLeg(0.65,0.50+0.05,0.75,0.505+0.05);
        //TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);
        TLegend *etacut = tdrLeg(0.65,0.50,0.75,0.505);
    
        if (sampleType==1) {
            sample->SetHeader("dijet sample");
        } else if (sampleType==2) {
            sample->SetHeader("#gamma+jet sample");
        } else if (sampleType==3) {
            sample->SetHeader("Z#mu#mu+jet sample");
        }
        
        //alphacut->SetHeader("#alpha<0.3");
        etacut->SetHeader("#left|#eta#right|< 2.5");

        hppHists.ptStack->GetXaxis()->SetNoExponent();
        hppHists.ptStack->GetXaxis()->SetMoreLogLabels(kTRUE);
        hppHists.ptStack->GetXaxis()->SetTitle("p_{T} (GeV)");
        hppHists.ptStack->GetYaxis()->SetTitle("Flavor fraction");
        
        leg->SetNColumns(3);
        //leg->SetHeader("CTEQ6L1/MRST LO**");
        leg->SetHeader("P8 P6 HW");

        leg->AddEntry(p8Hists.ptProjections[5]," ","f"); leg->AddEntry(p6Hists.ptProjections[5]," ","p"); leg->AddEntry(hppHists.ptProjections[5],"None ","p");
        leg->AddEntry(p8Hists.ptProjections[0]," ","f"); leg->AddEntry(p6Hists.ptProjections[0]," ","p"); leg->AddEntry(hppHists.ptProjections[0],"Gluon","p");
        leg->AddEntry(p8Hists.ptProjections[1]," ","f"); leg->AddEntry(p6Hists.ptProjections[1]," ","p"); leg->AddEntry(hppHists.ptProjections[1],"Light","p");
        leg->AddEntry(p8Hists.ptProjections[2]," ","f"); leg->AddEntry(p6Hists.ptProjections[2]," ","p"); leg->AddEntry(hppHists.ptProjections[2],"Strange","p");
        leg->AddEntry(p8Hists.ptProjections[3]," ","f"); leg->AddEntry(p6Hists.ptProjections[3]," ","p"); leg->AddEntry(hppHists.ptProjections[3],"Charm","p");
        leg->AddEntry(p8Hists.ptProjections[4]," ","f"); leg->AddEntry(p6Hists.ptProjections[4]," ","p"); leg->AddEntry(hppHists.ptProjections[4],"Bottom","p");
    }

    if (uptPlots) {
        TH1D *h1u = new TH1D("hu",";p_{T} (GeV);Fraction",ptBins,ptRange);
        setTDRStyle();
        gROOT->ForceStyle();
        gStyle->SetOptStat(kFALSE); //removes old legend
        gStyle->SetAxisColor(1, "XYZ");
        gStyle->SetStripDecimals(kTRUE);
        gStyle->SetTickLength(0.03, "XYZ");
        gStyle->SetNdivisions(510, "XYZ");
        gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
        gStyle->SetPadTickY(1);

        TCanvas *canv = tdrCanvas("c1u",h1u,12,0,1);
        h1u->GetYaxis()->SetTitleOffset(1.25);
        h1u->GetXaxis()->SetTitleOffset(1.0);
        h1u->GetXaxis()->SetLabelSize(0.045);
        h1u->GetYaxis()->SetLabelSize(0.045);
        h1u->GetXaxis()->SetTitleSize(0.045);
        h1u->GetYaxis()->SetTitleSize(0.045);
        gPad->SetLogx();
        
        p6Hists.uptStack->SetHistogram(h1u);
        p8Hists.uptStack->SetHistogram(h1u);
        hppHists.uptStack->SetHistogram(h1u);
        gStyle->SetOptStat(kFALSE); //removes old legend
        
        /* Generic stuff: */
        p6Hists.uptStack->SetMaximum(0.95);
        p8Hists.uptStack->SetMaximum(0.95);
        hppHists.uptStack->SetMaximum(0.95);

        /* Fraction plotting */
        h1u->GetXaxis()->SetRange(9,47);
        //h1->GetYaxis()->SetRangeUser(-0.001,1.001);
        p8Hists.uptStack->Draw("");
        hppHists.uptStack->Draw("sameP");
        p6Hists.uptStack->Draw("sameP");

        /* Fraction legends */

        //TLegend *leg = tdrLeg(0.2,0.41,0.5,0.91);
        TLegend *uleg = tdrLeg(0.3,0.4,0.75,0.85);
        TLegend *uheading = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
        TLegend *usample = tdrLeg(0.65,0.50+0.05,0.75,0.505+0.05);
        //TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);
        TLegend *uetacut = tdrLeg(0.65,0.50,0.75,0.505);
    
        if (sampleType==1) {
            usample->SetHeader("dijet sample");
        } else if (sampleType==2) {
            usample->SetHeader("#gamma+jet sample");
        } else if (sampleType==3) {
            usample->SetHeader("Z#mu#mu+jet sample");
        }
        
        //alphacut->SetHeader("#alpha<0.3");
        uetacut->SetHeader("#left|#eta#right|< 2.5");

        hppHists.uptStack->GetXaxis()->SetNoExponent();
        hppHists.uptStack->GetXaxis()->SetMoreLogLabels(kTRUE);
        hppHists.uptStack->GetXaxis()->SetTitle("p_{T} (GeV)");
        hppHists.uptStack->GetYaxis()->SetTitle("Flavor fraction");
        
        uleg->SetNColumns(3);
        //leg->SetHeader("CTEQ6L1/MRST LO**");
        uleg->SetHeader("P8 P6 HW");

        uleg->AddEntry(p8Hists.ptProjections[5]," ","f"); uleg->AddEntry(p6Hists.ptProjections[5]," ","p"); uleg->AddEntry(hppHists.ptProjections[5],"None ","p");
        uleg->AddEntry(p8Hists.ptProjections[0]," ","f"); uleg->AddEntry(p6Hists.ptProjections[0]," ","p"); uleg->AddEntry(hppHists.ptProjections[0],"Gluon","p");
        uleg->AddEntry(p8Hists.ptProjections[1]," ","f"); uleg->AddEntry(p6Hists.ptProjections[1]," ","p"); uleg->AddEntry(hppHists.ptProjections[1],"Light","p");
        uleg->AddEntry(p8Hists.ptProjections[2]," ","f"); uleg->AddEntry(p6Hists.ptProjections[2]," ","p"); uleg->AddEntry(hppHists.ptProjections[2],"Strange","p");
        uleg->AddEntry(p8Hists.ptProjections[3]," ","f"); uleg->AddEntry(p6Hists.ptProjections[3]," ","p"); uleg->AddEntry(hppHists.ptProjections[3],"Charm","p");
        uleg->AddEntry(p8Hists.ptProjections[4]," ","f"); uleg->AddEntry(p6Hists.ptProjections[4]," ","p"); uleg->AddEntry(hppHists.ptProjections[4],"Bottom","p");
    }

    if (etaPlots) {
        TH1D *h1eta = new TH1D("h eta",";p_{T} (GeV);Fraction",40,-etaMax,etaMax);
        setTDRStyle();
        gROOT->ForceStyle();
        gStyle->SetOptStat(kFALSE); //removes old legend
        gStyle->SetAxisColor(1, "XYZ");
        gStyle->SetStripDecimals(kTRUE);
        gStyle->SetTickLength(0.03, "XYZ");
        gStyle->SetNdivisions(510, "XYZ");
        gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
        gStyle->SetPadTickY(1);
        gStyle->SetOptLogx(0);

        TCanvas *canvEta = tdrCanvas("c1 eta",h1eta,12,0,1);
        h1eta->GetYaxis()->SetTitleOffset(1.25);
        h1eta->GetXaxis()->SetTitleOffset(1.0);
        h1eta->GetXaxis()->SetLabelSize(0.045);
        h1eta->GetYaxis()->SetLabelSize(0.045);
        h1eta->GetXaxis()->SetTitleSize(0.045);
        h1eta->GetYaxis()->SetTitleSize(0.045);

        hppHists.etaStack->SetHistogram(h1eta);
        p6Hists.etaStack->SetHistogram(h1eta);
        p8Hists.etaStack->SetHistogram(h1eta);

        gStyle->SetOptStat(kFALSE); //removes old legend

        p6Hists.etaStack->SetMaximum(0.95);
        p8Hists.etaStack->SetMaximum(0.95);
        hppHists.etaStack->SetMaximum(0.95);

        /* Fraction plotting */
        //h1eta->GetXaxis()->SetRange(9,47);
        //h1->GetYaxis()->SetRangeUser(-0.001,1.001);
        p8Hists.etaStack->Draw("");
        hppHists.etaStack->Draw("sameP");
        p6Hists.etaStack->Draw("sameP");

        /* Fraction legends */
        //heading->AddEntry()
        //TLegend *leg = tdrLeg(0.2,0.41,0.5,0.91);
        TLegend *legEta = tdrLeg(0.3,0.4,0.75,0.85);
        TLegend *headingEta = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
        TLegend *sampleEta = tdrLeg(0.65,0.50+0.05,0.75,0.505+0.05);
        //TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);

        if (sampleType==1) {
            sampleEta->SetHeader("dijet sample");
        } else if (sampleType==2) {
            sampleEta->SetHeader("#gamma+jet sample");
        } else if (sampleType==3) {
            sampleEta->SetHeader("Z#mu#mu+jet sample");
        }
        //heading->SetHeader("Pythia8 Simulation (4C Tune)");
        //alphacut->SetHeader("#alpha<0.3");

        hppHists.etaStack->GetXaxis()->SetNoExponent();
        hppHists.etaStack->GetXaxis()->SetTitle("#eta");
        hppHists.etaStack->GetYaxis()->SetTitle("Flavor fraction");

        legEta->SetNColumns(3);
        //legEta->SetHeader("CTEQ6L1/MRST LO**");
        legEta->SetHeader("P8 P6 Hw");

        legEta->AddEntry(p8Hists.etaProjections[5]," ","f"); legEta->AddEntry(p6Hists.etaProjections[5]," ","p"); legEta->AddEntry(hppHists.etaProjections[5],"None ","p");
        legEta->AddEntry(p8Hists.etaProjections[0]," ","f"); legEta->AddEntry(p6Hists.etaProjections[0]," ","p"); legEta->AddEntry(hppHists.etaProjections[0],"Gluon","p");
        legEta->AddEntry(p8Hists.etaProjections[1]," ","f"); legEta->AddEntry(p6Hists.etaProjections[1]," ","p"); legEta->AddEntry(hppHists.etaProjections[1],"Light","p");
        legEta->AddEntry(p8Hists.etaProjections[2]," ","f"); legEta->AddEntry(p6Hists.etaProjections[2]," ","p"); legEta->AddEntry(hppHists.etaProjections[2],"Strange","p");
        legEta->AddEntry(p8Hists.etaProjections[3]," ","f"); legEta->AddEntry(p6Hists.etaProjections[3]," ","p"); legEta->AddEntry(hppHists.etaProjections[3],"Charm","p");
        legEta->AddEntry(p8Hists.etaProjections[4]," ","f"); legEta->AddEntry(p6Hists.etaProjections[4]," ","p"); legEta->AddEntry(hppHists.etaProjections[4],"Bottom","p");
    }
    
    if (indicatorPlots) {
        /* Scaling */
        for (std::size_t i = 0; i < 3; ++i) {
            p6Hists.mult[i]->Scale(1/p6Hists.mult[i]->Integral());
            p6Hists.ptd[i]->Scale(1/p6Hists.ptd[i]->Integral());
            p6Hists.s2[i]->Scale(1/p6Hists.s2[i]->Integral());
            p8Hists.mult[i]->Scale(1/p8Hists.mult[i]->Integral());
            p8Hists.ptd[i]->Scale(1/p8Hists.ptd[i]->Integral());
            p8Hists.s2[i]->Scale(1/p8Hists.s2[i]->Integral());
            hppHists.mult[i]->Scale(1/hppHists.mult[i]->Integral());
            hppHists.ptd[i]->Scale(1/hppHists.ptd[i]->Integral());
            hppHists.s2[i]->Scale(1/hppHists.s2[i]->Integral());
        }
        
        /* Multiplicity */
        TH1D *h2 = new TH1D("h2",";Number of constituents;Events",60,0,60);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c2 = tdrCanvas("c2",h2,0,33,1);

        h2->SetMinimum(0);
        h2->SetMaximum(0.115);
        h2->GetYaxis()->SetNoExponent();
        h2->GetXaxis()->SetNoExponent();
        h2->GetXaxis()->SetRangeUser(0,60);
            
        tdrDraw(p8Hists.mult[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p6Hists.mult[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(hppHists.mult[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(p8Hists.mult[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.mult[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.mult[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        TLegend *multLeg = tdrLeg(0.5,0.6,0.8,0.9);
        multLeg->SetNColumns(3);

        multLeg->AddEntry(p8Hists.mult[0],"    ","l");
        multLeg->AddEntry(p6Hists.mult[0],"    ","f");
        multLeg->AddEntry(hppHists.mult[0],"    Gluon","p");
        multLeg->AddEntry(p8Hists.mult[1],"    ","l");
        multLeg->AddEntry(p6Hists.mult[1],"    ","f");
        multLeg->AddEntry(hppHists.mult[1],"    Quark","p");
        multLeg->SetHeader("P8 P6 HW");

        /* PTD */
        TH1D *h3 = new TH1D("h3",";p_{T}D;Events",100,0,1);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c3 = tdrCanvas("c3",h3,0,33,1);

        h3->SetMaximum(0.075);
        h3->GetYaxis()->SetNoExponent();
        h3->GetXaxis()->SetNoExponent();
        h3->GetXaxis()->SetRangeUser(0,30);

        tdrDraw(p8Hists.ptd[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p6Hists.ptd[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(hppHists.ptd[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(p8Hists.ptd[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.ptd[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.ptd[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        TLegend *ptdLeg = tdrLeg(0.5,0.6,0.8,0.9);
        ptdLeg->SetNColumns(3);

        ptdLeg->AddEntry(p8Hists.ptd[0],"    ","l");
        ptdLeg->AddEntry(p6Hists.ptd[0],"    ","f");
        ptdLeg->AddEntry(hppHists.ptd[0],"    Gluon","p");
        ptdLeg->AddEntry(p8Hists.ptd[1],"    ","l");
        ptdLeg->AddEntry(p6Hists.ptd[1],"    ","f");
        ptdLeg->AddEntry(hppHists.ptd[1],"    Quark","p");
        ptdLeg->SetHeader("P8 P6 HW");

        /* Sigma2 */
        TH1D *h4 = new TH1D("h4",";#sigma_{2};Events",100,0,0.2);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c4 = tdrCanvas("c4",h4,0,33,1);
        h4->SetMinimum(0);
        h4->SetMaximum(0.085);
        h4->GetYaxis()->SetNoExponent();
        h4->GetXaxis()->SetNoExponent();
        h4->GetXaxis()->SetRangeUser(0,0.2);

        tdrDraw(p8Hists.s2[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p6Hists.s2[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(hppHists.s2[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(p8Hists.s2[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.s2[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.s2[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);

        TLegend *s2Leg = tdrLeg(0.5,0.6,0.8,0.9);
        s2Leg->SetNColumns(3);

        s2Leg->AddEntry(p8Hists.s2[0],"    ","l");
        s2Leg->AddEntry(p6Hists.s2[0],"    ","f");
        s2Leg->AddEntry(hppHists.s2[0],"    Gluon","p");
        s2Leg->AddEntry(p8Hists.s2[1],"    ","l");
        s2Leg->AddEntry(p6Hists.s2[1],"    ","f");
        s2Leg->AddEntry(hppHists.s2[1],"    Quark","p");
        s2Leg->SetHeader("P8 P6 HW");
    }

    if (debugPlots) {
        /* Pt */
        TH1D *h5 = new TH1D("h5","pT;pT;Events",ptBins,ptRange);
        gStyle->SetOptLogy(1);
        TCanvas *c5 = tdrCanvas("c5",h5,0,33,1);
        h5->SetMinimum(0.0000000000001);
        h5->SetMaximum(1000);
        h5->GetYaxis()->SetNoExponent();
        h5->GetXaxis()->SetNoExponent();
        //     h5->GetYaxis()->SetMoreLogLabels();
        h5->GetXaxis()->SetMoreLogLabels();

        p8Pt->Scale(1/p8Pt->Integral());
        p6Pt->Scale(1/p6Pt->Integral());
        hppPt->Scale(1/hppPt->Integral());
        tdrDraw(p8Pt,"HIST",kFullTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        tdrDraw(p6Pt,"HIST P",kOpenTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        tdrDraw(hppPt,"HIST P",kFullTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        TLegend *ptLeg = tdrLeg(0.5,0.6,0.8,0.9);

        ptLeg->AddEntry(p8Pt,"Pythia8","f");
        ptLeg->AddEntry(p6Pt,"Pythia6","p");
        ptLeg->AddEntry(hppPt,"Herwig","p");

        /* Weights */
        TH1D *h6 = new TH1D("h6","weights;pT;weight",ptBins,ptRange);
        TCanvas *c6 = tdrCanvas("c6",h6,0,33,1);
        h6->SetMaximum(0.04);
        h6->SetMinimum(0.0000000001);

        tdrDraw(hppWeighting->ProjectionX(),"E",kFullTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        tdrDraw(p6Weighting->ProjectionX(),"E",kOpenTriangleUp,kRed,kSolid,-1,3003,kRed);
        tdrDraw(p8Weighting->ProjectionX(),"E",kFullTriangleUp,kRed,kSolid,-1,3003,kRed);
        gPad->SetLogx();
        gPad->SetLogy();
    }

    if (separatorPlots) {
        /* DR */
        TH1D *h7 = new TH1D("h7","dR;dR;Events",drBins,drLower,drUpper);
        gStyle->SetOptLogy(1);
        gStyle->SetOptLogx(0);
        TCanvas *c7 = tdrCanvas("c7",h7,0,33,1);
        h7->SetMinimum(0.0000001);
        h7->SetMaximum(1);

        p8Hists.dR->Scale(1/p8Hists.dR->Integral());
        p6Hists.dR->Scale(1/p6Hists.dR->Integral());
        hppHists.dR->Scale(1/hppHists.dR->Integral());
        tdrDraw(p8Hists.dR,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.dR,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.dR,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        TLegend *drLeg = tdrLeg(0.5,0.6,0.8,0.9);
        
        drLeg->AddEntry(p8Hists.dR,"Pythia 8","l");
        drLeg->AddEntry(p6Hists.dR,"Pythia 6","l");
        drLeg->AddEntry(hppHists.dR,"Herwig++","l");
        
        cout << "p8 over 0.3: " << p8Hists.dR->Integral(30,500)/p8Hists.dR->Integral() << endl;
        cout << "p6 over 0.3: " << p6Hists.dR->Integral(30,500)/p6Hists.dR->Integral() << endl;
        cout << "hwpp over 0.3: " << hppHists.dR->Integral(30,500)/hppHists.dR->Integral() << endl;
        
        cout << "p8 over 0.5: " << p8Hists.dR->Integral(50,500)/p8Hists.dR->Integral() << endl;
        cout << "p6 over 0.5: " << p6Hists.dR->Integral(50,500)/p6Hists.dR->Integral() << endl;
        cout << "hwpp over 0.5: " << hppHists.dR->Integral(50,500)/hppHists.dR->Integral() << endl;
        
        /* Alpha */
        TH1D *h8 = new TH1D("h8",";alpha;Events",alphaBins,alphaLower,alphaUpper);
        gStyle->SetOptLogy(0);
        gStyle->SetOptLogx(0);
        TCanvas *c8 = tdrCanvas("c8",h8,0,33,1);
        h8->SetMinimum(0.0);
        h8->SetMaximum(0.012);

        p8Hists.alpha->Scale(1/p8Hists.alpha->Integral());
        p6Hists.alpha->Scale(1/p6Hists.alpha->Integral());
        hppHists.alpha->Scale(1/hppHists.alpha->Integral());
        tdrDraw(p8Hists.alpha,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.alpha,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.alpha,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        TLegend *alphaLeg = tdrLeg(0.5,0.6,0.8,0.9);
        
        alphaLeg->AddEntry(p8Hists.alpha,"Pythia 8","l");
        alphaLeg->AddEntry(p6Hists.alpha,"Pythia 6","l");
        alphaLeg->AddEntry(hppHists.alpha,"Herwig++","l");
        
        /* dPhi */
        TH1D *h9 = new TH1D("h9",";dphi;Events",dphiBins,dphiLower,dphiUpper);
        gStyle->SetOptLogy(0);
        gStyle->SetOptLogx(0);
        TCanvas *c9 = tdrCanvas("c9",h9,0,33,1);
        h9->SetMinimum(0.0);
        h9->SetMaximum(0.03);

        p8Hists.dPhi->Scale(1/p8Hists.dPhi->Integral());
        p6Hists.dPhi->Scale(1/p6Hists.dPhi->Integral());
        hppHists.dPhi->Scale(1/hppHists.dPhi->Integral());
        tdrDraw(p8Hists.dPhi,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(p6Hists.dPhi,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.dPhi,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        TLegend *dphiLeg = tdrLeg(0.5,0.6,0.8,0.9);
        
        dphiLeg->AddEntry(p8Hists.dPhi,"Pythia 8","l");
        dphiLeg->AddEntry(p6Hists.dPhi,"Pythia 6","l");
        dphiLeg->AddEntry(hppHists.dPhi,"Herwig++","l");
    }
}
