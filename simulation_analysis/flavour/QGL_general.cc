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
const bool indicatorPlots = true;
const bool separatorPlots = true;

const int sampleType=1;
const bool hard=false; // Pythia6 logistics

const int ptBins = 48.;
const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

// Backup cuts for eta and pt
const double ptUpper = 100;
const double ptLower = 80;

const double etaMax = 1.3;

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
    TH1D* dR;
    TH1D* alpha;
    TH1D* dPhi;

    vector<TH1D*> mult;
    vector<TH1D*> ptd;
    vector<TH1D*> s2;
};

bool ReadTree(string file, HistStore& hists)
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

            hists.dR->Fill(mDR[i],mWeight);
            hists.alpha->Fill(mAlpha[i],mWeight);
            hists.dPhi->Fill(mDPhi[i],mWeight);

            // Upper and lower limit safety cut, pt
            if(tmpVec.Pt()>ptUpper || tmpVec.Pt()<ptLower) continue;

            ++mCount;

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

void QGL_general() {
    float x_min[] = {0, 0.15, 0}; 
    float x_max[] = {0.16, 1, 60};
    float y_min[] = {0, 0, 0}; 
    float y_max[] = {0.06, 0.06, 0.07};
    int bins[] = {100, 100, 60};

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    cout << "Order: herwig++, pythia6, pythia8" << endl;

    vector<string> HWFiles;
    HWFiles.push_back("cteq6l1/herwig_dijet_physics_1000000.root");
    HWFiles.push_back("cteq6l1/herwig_gammajet_physics_1000000.root");
    HWFiles.push_back("cteq6l1/herwig_Zjet_physics_1000000.root");
    vector<string> P6Files;
    P6Files.push_back("cteq6l1/pythia6_dijet_physics_1000000.root");
    P6Files.push_back("cteq6l1/pythia6_gammajet_physics_1000000.root");
    P6Files.push_back("cteq6l1/pythia6_Zjet_physics_1000000.root");
    vector<string> P8Files;
    P8Files.push_back("cteq6l1/pythia8_dijet_physics_1000000.root");
    P8Files.push_back("cteq6l1/pythia8_gammajet_physics_1000000.root");
    P8Files.push_back("cteq6l1/pythia8_Zjet_physics_1000000.root");

    /* Herwig */
    HistStore hppHists;
    
    for (std::size_t i = 0; i < 3; ++i) {
        hppHists.mult.push_back(new TH1D("","",bins[2],x_min[2],x_max[2]));
        hppHists.ptd.push_back(new TH1D("","",bins[1],x_min[1],x_max[1]));
        hppHists.s2.push_back(new TH1D("","",bins[0],x_min[0],x_max[0]));
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
    
    for ( auto file : HWFiles ) {
        if (!ReadTree(file, hppHists)) {
            cout << "Could not load Herwig data" << endl;
            return;
        }
    }

    /* Pythia6 */
    HistStore p6Hists;
    
    for (std::size_t i = 0; i < 3; ++i) { 
        p6Hists.mult.push_back(new TH1D("","",bins[2],x_min[2],x_max[2]));
        p6Hists.ptd.push_back(new TH1D("","",bins[1],x_min[1],x_max[1]));
        p6Hists.s2.push_back(new TH1D("","",bins[0],x_min[0],x_max[0]));
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
    
    for ( auto file : P6Files ) {
        if (!ReadTree(file, p6Hists)) {
            cout << "Could not load Pythia6 data" << endl;
            return;
        }
    }
    
    /* Pythia8 plot initialization */
    HistStore p8Hists;

    for (std::size_t i = 0; i < 3; ++i) { 
        p8Hists.mult.push_back(new TH1D("","",bins[2],x_min[2],x_max[2]));
        p8Hists.ptd.push_back(new TH1D("","",bins[1],x_min[1],x_max[1]));
        p8Hists.s2.push_back(new TH1D("","",bins[0],x_min[0],x_max[0]));
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
    
    for ( auto file : P8Files ) {
        if (!ReadTree(file, p8Hists)) {
            cout << "Could not load Pythia6 data" << endl;
            return;
        }
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
        tdrDraw(p6Hists.dR,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.dR,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        tdrDraw(p8Hists.dR,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
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
        tdrDraw(p6Hists.alpha,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.alpha,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        tdrDraw(p8Hists.alpha,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
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
        tdrDraw(p6Hists.dPhi,"HIST",kOpenTriangleDown,kRed,kSolid,-1,300,kRed);
        tdrDraw(hppHists.dPhi,"HIST",kFullTriangleDown,kGreen+3,kSolid,-1,3003,kGreen+3);
        tdrDraw(p8Hists.dPhi,"HIST",kFullTriangleDown,kBlue,kSolid,-1,300,kBlue);
        TLegend *dphiLeg = tdrLeg(0.5,0.6,0.8,0.9);
        
        dphiLeg->AddEntry(p8Hists.dPhi,"Pythia 8","l");
        dphiLeg->AddEntry(p6Hists.dPhi,"Pythia 6","l");
        dphiLeg->AddEntry(hppHists.dPhi,"Herwig++","l");
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
        TH1D *h2 = new TH1D("h2",";Number of constituents;Events",bins[2],x_min[2],x_max[2]);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c2 = tdrCanvas("c2",h2,0,33,1);

        h2->SetMinimum(y_min[2]);
        h2->SetMaximum(y_max[2]);
        h2->GetYaxis()->SetNoExponent();
        h2->GetXaxis()->SetNoExponent();

        tdrDraw(p6Hists.mult[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(p6Hists.mult[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.mult[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(hppHists.mult[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        tdrDraw(p8Hists.mult[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p8Hists.mult[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        TLegend *multLeg = tdrLeg(0.65,0.7,0.85,0.9);
        multLeg->SetNColumns(3);

        multLeg->AddEntry(p8Hists.mult[0],"    ","l");
        multLeg->AddEntry(p6Hists.mult[0],"    ","f");
        multLeg->AddEntry(hppHists.mult[0],"    Gluon","p");
        multLeg->AddEntry(p8Hists.mult[1],"    ","l");
        multLeg->AddEntry(p6Hists.mult[1],"    ","f");
        multLeg->AddEntry(hppHists.mult[1],"    Quark","p");
        multLeg->SetHeader("P8 P6 HW");

        /* PTD */
        TH1D *h3 = new TH1D("h3",";p_{T}D;Events",bins[1],x_min[1],x_max[1]);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c3 = tdrCanvas("c3",h3,0,33,1);

        h3->SetMinimum(y_min[1]);
        h3->SetMaximum(y_max[1]);
        h3->GetYaxis()->SetNoExponent();
        h3->GetXaxis()->SetNoExponent();

        tdrDraw(p6Hists.ptd[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(p6Hists.ptd[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.ptd[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(hppHists.ptd[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        tdrDraw(p8Hists.ptd[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p8Hists.ptd[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
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
        TH1D *h4 = new TH1D("h4",";#sigma_{2};Events",bins[0],x_min[0],x_max[0]);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c4 = tdrCanvas("c4",h4,0,33,1);
        h4->SetMinimum(y_min[0]);
        h4->SetMaximum(y_max[0]);
        h4->GetYaxis()->SetNoExponent();
        h4->GetXaxis()->SetNoExponent();

        tdrDraw(p6Hists.s2[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(p6Hists.s2[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(hppHists.s2[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(hppHists.s2[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        tdrDraw(p8Hists.s2[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(p8Hists.s2[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
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
}
