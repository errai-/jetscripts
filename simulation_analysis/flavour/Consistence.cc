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
const bool uptPlots = true;
const bool etaPlots = false;
const bool debugPlots = false;
const bool indicatorPlots = false;
const bool separatorPlots = false;

const int sampleType=3;
const bool hard=false; // Pythia6 logistics

const int ptBins = 48.;
const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

// Backup cuts for eta and pt
const double ptUpper = 2000;
const double ptLower = 0;

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


void Consistence(string diFile, string gFile, string ZFile) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    cout << "Order: herwig++, pythia6, pythia8" << endl;

    /* Dijet data */
    HistStore diHists;

    for (std::size_t i = 0; i < 3; ++i) { 
        diHists.mult.push_back(new TH1D("","",60,0,60));
        diHists.ptd.push_back(new TH1D("","",100,0,1));
        diHists.s2.push_back(new TH1D("","",100,0,0.2));
    }

    for (TH1D* mult : hppHists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : hppHists.ptd) { ptd->Sumw2(); }
    for (TH1D* s2 : hppHists.s2) { s2->Sumw2(); }

    diHists.mult[0]->SetLineColor(kRed); diHists.mult[2]->SetLineColor(kGreen);
    diHists.ptd[0]->SetLineColor(kRed); diHists.ptd[2]->SetLineColor(kGreen);
    diHists.s2[0]->SetLineColor(kRed); diHists.s2[2]->SetLineColor(kGreen);

    if (!ReadTree(diFile, diHists)) {
        cout << "Could not load Herwig data" << endl;
        return;
    }

    /* Gammajet data */
    HistStore gHists;

    for (std::size_t i = 0; i < 3; ++i) { 
        gHists.mult.push_back(new TH1D("","",60,0,60));
        gHists.ptd.push_back(new TH1D("","",100,0,1));
        gHists.s2.push_back(new TH1D("","",100,0,0.2));
    }

    for (TH1D* mult : gHists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : gHists.ptd) { ptd->Sumw2(); }
    for (TH1D* sigma2 : gHists.s2) { sigma2->Sumw2(); }

    gHists.mult[0]->SetLineColor(kRed); gHists.mult[2]->SetLineColor(kGreen);
    gHists.ptd[0]->SetLineColor(kRed); gHists.ptd[2]->SetLineColor(kGreen);
    gHists.s2[0]->SetLineColor(kRed); gHists.s2[2]->SetLineColor(kGreen);

    if (!ReadTree(gFile, gHists)) {
        cout << "Could not load Pythia6 data" << endl;
        return;
    }

    /* Zjet data */
    HistStore zHists;

    for (std::size_t i = 0; i < 3; ++i) { 
        zHists.mult.push_back(new TH1D("","",60,0,60));
        zHists.ptd.push_back(new TH1D("","",100,0,1));
        zHists.s2.push_back(new TH1D("","",100,0,0.2));
    }

    for (TH1D* mult : zHists.mult) { mult->Sumw2(); }
    for (TH1D* ptd : zHists.ptd) { ptd->Sumw2(); }
    for (TH1D* sigma2 : zHists.s2) { sigma2->Sumw2(); }

    zHists.mult[0]->SetLineColor(kRed); zHists.mult[2]->SetLineColor(kGreen);
    zHists.ptd[0]->SetLineColor(kRed); zHists.ptd[2]->SetLineColor(kGreen);
    zHists.s2[0]->SetLineColor(kRed); zHists.s2[2]->SetLineColor(kGreen);

    if (!ReadTree(zFile, zHists)) {
        cout << "Could not load Pythia6 data" << endl;
        return;
    }

    if (indicatorPlots) {
        /* Scaling */
        for (std::size_t i = 0; i < 3; ++i) {
            diHists.mult[i]->Scale(1/diHists.mult[i]->Integral());
            diHists.ptd[i]->Scale(1/diHists.ptd[i]->Integral());
            diHists.s2[i]->Scale(1/diHists.s2[i]->Integral());
            gHists.mult[i]->Scale(1/gHists.mult[i]->Integral());
            gHists.ptd[i]->Scale(1/gHists.ptd[i]->Integral());
            gHists.s2[i]->Scale(1/gHists.s2[i]->Integral());
            zHists.mult[i]->Scale(1/zHists.mult[i]->Integral());
            zHists.ptd[i]->Scale(1/zHists.ptd[i]->Integral());
            zHists.s2[i]->Scale(1/zHists.s2[i]->Integral());
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

        tdrDraw(diHists.mult[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(gHists.mult[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(zHists.mult[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(diHists.mult[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(gHists.mult[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(zHists.mult[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        TLegend *multLeg = tdrLeg(0.5,0.6,0.8,0.9);
        multLeg->SetNColumns(3);

        multLeg->AddEntry(diHists.mult[0],"    ","l");
        multLeg->AddEntry(gHists.mult[0],"    ","f");
        multLeg->AddEntry(zHists.mult[0],"    Gluon","p");
        multLeg->AddEntry(diHists.mult[1],"    ","l");
        multLeg->AddEntry(gHists.mult[1],"    ","f");
        multLeg->AddEntry(zHists.mult[1],"    Quark","p");
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

        tdrDraw(diHists.ptd[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(gHists.ptd[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(zHists.ptd[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(diHists.ptd[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(gHists.ptd[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(zHists.ptd[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        TLegend *ptdLeg = tdrLeg(0.5,0.6,0.8,0.9);
        ptdLeg->SetNColumns(3);

        ptdLeg->AddEntry(diHists.ptd[0],"    ","l");
        ptdLeg->AddEntry(gHists.ptd[0],"    ","f");
        ptdLeg->AddEntry(zHists.ptd[0],"    Gluon","p");
        ptdLeg->AddEntry(diHists.ptd[1],"    ","l");
        ptdLeg->AddEntry(gHists.ptd[1],"    ","f");
        ptdLeg->AddEntry(zHists.ptd[1],"    Quark","p");
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

        tdrDraw(diHists.s2[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,300,kRed-3);
        tdrDraw(gHists.s2[0],"HIST",kOpenTriangleDown,kCyan,kSolid,-1,3003,kCyan);
        tdrDraw(zHists.s2[0],"HIST P",kFullTriangleDown,kGreen+3,kSolid,-1,300,kGreen+3);
        tdrDraw(diHists.s2[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,300,kBlue);
        tdrDraw(gHists.s2[1],"HIST",kOpenTriangleUp,kOrange+7,kSolid,-1,3005,kOrange+7);
        tdrDraw(zHists.s2[1],"HIST P",kFullTriangleUp,kMagenta,kSolid,-1,300,kMagenta);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);

        TLegend *s2Leg = tdrLeg(0.5,0.6,0.8,0.9);
        s2Leg->SetNColumns(3);

        s2Leg->AddEntry(diHists.s2[0],"    ","l");
        s2Leg->AddEntry(gHists.s2[0],"    ","f");
        s2Leg->AddEntry(zHists.s2[0],"    Gluon","p");
        s2Leg->AddEntry(diHists.s2[1],"    ","l");
        s2Leg->AddEntry(gHists.s2[1],"    ","f");
        s2Leg->AddEntry(zHists.s2[1],"    Quark","p");
        s2Leg->SetHeader("P8 P6 HW");
    }

}
