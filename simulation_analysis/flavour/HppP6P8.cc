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

const int ptBins = 48.;
const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};
const double ptUpper = 2000;
const double ptLower = 0;
const double etaMax = 2.5;//1.3;
const int sampleType=1;
const int isHerwig=0;

void HppP6P8(string herwFile, string p6File, string p8File, bool debugPlots = false, bool indicatorPlots = false) {
/* Init branches */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    
    static const Int_t kMaxfJets = 100;
    
    Int_t           hppJets_;
    Double_t        hppX[kMaxfJets];   //[hppJets_]
    Double_t        hppY[kMaxfJets];   //[hppJets_]
    Double_t        hppZ[kMaxfJets];   //[hppJets_]
    Double_t        hppT[kMaxfJets];   //[hppJets_]
    Int_t           p6Jets_;
    Double_t        p6X[kMaxfJets];   //[p6Jets_]
    Double_t        p6Y[kMaxfJets];   //[p6Jets_]
    Double_t        p6Z[kMaxfJets];   //[p6Jets_]
    Double_t        p6T[kMaxfJets];   //[p6Jets_]
    Int_t           p8Jets_;
    Double_t        p8X[kMaxfJets];   //[p8Jets_]
    Double_t        p8Y[kMaxfJets];   //[p8Jets_]
    Double_t        p8Z[kMaxfJets];   //[p8Jets_]
    Double_t        p8T[kMaxfJets];   //[p8Jets_]
    
    Double_t        hppWeight;
    Double_t        p6Weight;
    Double_t        p8Weight;
    Int_t           hppFlav[kMaxfJets];   //[hppJets_]
    Int_t           p6Flav[kMaxfJets];   //[p6Jets_]
    Int_t           p8Flav[kMaxfJets];   //[p8Jets_]
    
    Int_t           hppConstituents[kMaxfJets];   //[hppJets_]
    Double_t        hppPTD[kMaxfJets];   //[hppJets_]
    Double_t        hppSigma2[kMaxfJets];   //[hppJets_]
    Int_t           p6Constituents[kMaxfJets];   //[p6Jets_]
    Double_t        p6PTD[kMaxfJets];   //[p6Jets_]
    Double_t        p6Sigma2[kMaxfJets];   //[p6Jets_]
    Int_t           p8Constituents[kMaxfJets];   //[p8Jets_]
    Double_t        p8PTD[kMaxfJets];   //[p8Jets_]
    Double_t        p8Sigma2[kMaxfJets];   //[p8Jets_]
    
/* Tree setup */
    TTree* hppTree; TTree* p6Tree; TTree* p8Tree;
    
    TChain * hppChain = new TChain("JetTree","");
    hppChain->Add(herwFile.c_str()); hppTree = hppChain;
    TChain * p6Chain = new TChain("JetTree","");
    p6Chain->Add(p6File.c_str()); p6Tree = p6Chain;
    TChain * p8Chain = new TChain("JetTree","");
    p8Chain->Add(p8File.c_str()); p8Tree = p8Chain;
    
    if (!p8Tree || !p6Tree || !hppTree) return;
    hppTree->SetMakeClass(1); p6Tree->SetMakeClass(1); p8Tree->SetMakeClass(1);
    

    hppTree->SetBranchAddress("fJets", &hppJets_);
    hppTree->SetBranchAddress("fJets.fP4.fCoordinates.fX", hppX);
    hppTree->SetBranchAddress("fJets.fP4.fCoordinates.fY", hppY);
    hppTree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", hppZ);
    hppTree->SetBranchAddress("fJets.fP4.fCoordinates.fT", hppT);
    hppTree->SetBranchAddress("fWeight", &hppWeight);
    hppTree->SetBranchAddress("fJets.fFlav", hppFlav);
    hppTree->SetBranchAddress("fJets.fConstituents", hppConstituents);
    hppTree->SetBranchAddress("fJets.fPTD", hppPTD);
    hppTree->SetBranchAddress("fJets.fSigma2", hppSigma2);

    p6Tree->SetBranchAddress("fJets", &p6Jets_);
    p6Tree->SetBranchAddress("fJets.fP4.fCoordinates.fX", p6X);
    p6Tree->SetBranchAddress("fJets.fP4.fCoordinates.fY", p6Y);
    p6Tree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", p6Z);
    p6Tree->SetBranchAddress("fJets.fP4.fCoordinates.fT", p6T);
    p6Tree->SetBranchAddress("fWeight", &p6Weight);
    p6Tree->SetBranchAddress("fJets.fFlav", p6Flav);
    p6Tree->SetBranchAddress("fJets.fConstituents", p6Constituents);
    p6Tree->SetBranchAddress("fJets.fPTD", p6PTD);
    p6Tree->SetBranchAddress("fJets.fSigma2", p6Sigma2);

    p8Tree->SetBranchAddress("fJets", &p8Jets_);
    p8Tree->SetBranchAddress("fJets.fP4.fCoordinates.fX", p8X);
    p8Tree->SetBranchAddress("fJets.fP4.fCoordinates.fY", p8Y);
    p8Tree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", p8Z);
    p8Tree->SetBranchAddress("fJets.fP4.fCoordinates.fT", p8T);
    p8Tree->SetBranchAddress("fWeight", &p8Weight);
    p8Tree->SetBranchAddress("fJets.fFlav", p8Flav);
    p8Tree->SetBranchAddress("fJets.fConstituents", p8Constituents);
    p8Tree->SetBranchAddress("fJets.fPTD", p8PTD);
    p8Tree->SetBranchAddress("fJets.fSigma2", p8Sigma2);
    
/* Herwig plot initialization */
    TProfile hppGluonFrac("hpp g","hpp g",ptBins,ptRange);
    TProfile hppLightquarkFrac("hpp lq","hpp lq",ptBins,ptRange);
    TProfile hppStrangeFrac("hpp s","hpp s",ptBins,ptRange);
    TProfile hppCharmFrac("hpp c","hpp c",ptBins,ptRange);
    TProfile hppBottomFrac("hpp b","hpp b",ptBins,ptRange);
    TProfile hppUnmatchedFrac("hpp unmatched","hpp unmatched",ptBins,ptRange);
    
    TProfile hppEtaGluonFrac("hpp eta g","hpp eta g",40,-etaMax,etaMax);
    TProfile hppEtaLightquarkFrac("hpp eta lq","hpp eta lq",40,-etaMax,etaMax);
    TProfile hppEtaStrangeFrac("hpp eta s","hpp eta s",40,-etaMax,etaMax);
    TProfile hppEtaCharmFrac("hpp eta c","hpp eta c",40,-etaMax,etaMax);
    TProfile hppEtaBottomFrac("hpp eta b","hpp eta b",40,-etaMax,etaMax);
    TProfile hppEtaUnmatchedFrac("hpp eta unmatched","hpp eta unmatched",40,-etaMax,etaMax);
    
    vector<TH1D*> hppMultHist;
    hppMultHist.push_back(new TH1D("hpp multiplicity g","hpp multiplicity g",60,0,60));
    hppMultHist.push_back(new TH1D("hpp multiplicity q","hpp multiplicity q",60,0,60)); 
    hppMultHist.push_back(new TH1D("hpp multiplicity u","hpp multiplicity u",60,0,60)); 
    for (TH1D* mult : hppMultHist) { mult->Sumw2(); }
    hppMultHist[0]->SetLineColor(kRed); hppMultHist[2]->SetLineColor(kGreen);

    vector<TH1D*> hppPTDHist;
    hppPTDHist.push_back(new TH1D("hpp pTD g","hpp pTD g",100,0,1));
    hppPTDHist.push_back(new TH1D("hpp pTD q","hpp pTD q",100,0,1));
    hppPTDHist.push_back(new TH1D("hpp pTD u","hpp pTD u",100,0,1));
    for (TH1D* ptd : hppPTDHist) { ptd->Sumw2(); }
    hppPTDHist[0]->SetLineColor(kRed); hppPTDHist[2]->SetLineColor(kGreen);

    vector<TH1D*> hppS2Hist;
    hppS2Hist.push_back(new TH1D("hpp sigma2 g","hpp sigma2 g",100,0,0.2));
    hppS2Hist.push_back(new TH1D("hpp sigma2 q","hpp sigma2 q",100,0,0.2));
    hppS2Hist.push_back(new TH1D("hpp sigma2 u","hpp sigma2 u",100,0,0.2));
    for (TH1D* sigma2 : hppS2Hist) { sigma2->Sumw2(); }
    hppS2Hist[0]->SetLineColor(kRed); hppS2Hist[2]->SetLineColor(kGreen);
    
    TH1D* hppPt = new TH1D("hpp Pt",";pT;Events",ptBins,ptRange);
    
/* Pythia6 plot initialization */
    TProfile p6GluonFrac("p6 g","",ptBins,ptRange);
    TProfile p6LightquarkFrac("p6 lq","",ptBins,ptRange);
    TProfile p6StrangeFrac("p6 s","",ptBins,ptRange);
    TProfile p6CharmFrac("p6 c","",ptBins,ptRange);
    TProfile p6BottomFrac("p6 b","",ptBins,ptRange);
    TProfile p6UnmatchedFrac("p6 unmatched","",ptBins,ptRange);
    
    TProfile p6EtaGluonFrac("p6 eta g","",40,-etaMax,etaMax);
    TProfile p6EtaLightquarkFrac("p6 eta lq","",40,-etaMax,etaMax);
    TProfile p6EtaStrangeFrac("p6 eta s","",40,-etaMax,etaMax);
    TProfile p6EtaCharmFrac("p6 eta c","",40,-etaMax,etaMax);
    TProfile p6EtaBottomFrac("p6 eta b","",40,-etaMax,etaMax);
    TProfile p6EtaUnmatchedFrac("p6 eta unmatched","",40,-etaMax,etaMax);
    
    vector<TH1D*> p6MultHist;
    p6MultHist.push_back(new TH1D("p6 multiplicity g","",60,0,60));
    p6MultHist.push_back(new TH1D("p6 multiplicity q","",60,0,60));
    p6MultHist.push_back(new TH1D("p6 multiplicity_u","",60,0,60));
    for (TH1D* mult : p6MultHist) { mult->Sumw2(); }
    p6MultHist[0]->SetLineColor(kRed); p6MultHist[2]->SetLineColor(kGreen);

    vector<TH1D*> p6PTDHist;
    p6PTDHist.push_back(new TH1D("p6 PTD g","p6 PTD g",100,0,1));
    p6PTDHist.push_back(new TH1D("p6 pTD q","p6 pTD q",100,0,1));
    p6PTDHist.push_back(new TH1D("p6 pTD u","p6 pTD u",100,0,1));
    for (TH1D* ptd : p6PTDHist) { ptd->Sumw2(); }
    p6PTDHist[0]->SetLineColor(kRed); p6PTDHist[2]->SetLineColor(kGreen);

    vector<TH1D*> p6S2Hist;
    p6S2Hist.push_back(new TH1D("p6 sigma2 g","p6 sigma2 g",100,0,0.2));
    p6S2Hist.push_back(new TH1D("p6 sigma2 q","p6 sigma2 q",100,0,0.2));
    p6S2Hist.push_back(new TH1D("p6 sigma2 u","p6 sigma2 u",100,0,0.2));
    for (TH1D* sigma2 : p6S2Hist) { sigma2->Sumw2(); }
    p6S2Hist[0]->SetLineColor(kRed); p6S2Hist[2]->SetLineColor(kGreen);
    
    TH1D* p6Pt = new TH1D("p6 Pt",";pT;Events",ptBins,ptRange);
    
/* Pythia8 plot initialization */
    TProfile p8GluonFrac("p8 g","p8 g",ptBins,ptRange);
    TProfile p8LightquarkFrac("p8 lq","p8 lq",ptBins,ptRange);
    TProfile p8StrangeFrac("p8 s","p8 s",ptBins,ptRange);
    TProfile p8CharmFrac("p8 c","p8 c",ptBins,ptRange);
    TProfile p8BottomFrac("p8 b","p8 b",ptBins,ptRange);
    TProfile p8UnmatchedFrac("p8 unmatched","p8 unmatched",ptBins,ptRange);
    
    TProfile p8EtaGluonFrac("p8 eta g","p8 eta g",40,-etaMax,etaMax);
    TProfile p8EtaLightquarkFrac("p8 eta lq","p8 eta lq",40,-etaMax,etaMax);
    TProfile p8EtaStrangeFrac("p8 eta s","p8 eta s",40,-etaMax,etaMax);
    TProfile p8EtaCharmFrac("p8 eta c","p8 eta c",40,-etaMax,etaMax);
    TProfile p8EtaBottomFrac("p8 eta b","p8 eta b",40,-etaMax,etaMax);
    TProfile p8EtaUnmatchedFrac("p8 eta unmatched","p8 eta unmatched",40,-etaMax,etaMax);
    
    vector<TH1D*> p8MultHist;
    p8MultHist.push_back(new TH1D("p8 multiplicity g","p8 multiplicity g",60,0,60));
    p8MultHist.push_back(new TH1D("p8 multiplicity q","p8 multiplicity q",60,0,60));
    p8MultHist.push_back(new TH1D("multiplicity_u","multiplicity_u",60,0,60));
    for (TH1D* mult : p8MultHist) { mult->Sumw2(); }
    p8MultHist[0]->SetLineColor(kRed); p8MultHist[2]->SetLineColor(kGreen);

    vector<TH1D*> p8PTDHist;
    p8PTDHist.push_back(new TH1D("p8 PTD g","p8 PTD g",100,0,1));
    p8PTDHist.push_back(new TH1D("p8 pTD q","p8 pTD q",100,0,1));
    p8PTDHist.push_back(new TH1D("p8 pTD u","p8 pTD u",100,0,1));
    for (TH1D* ptd : p8PTDHist) { ptd->Sumw2(); }
    p8PTDHist[0]->SetLineColor(kRed); p8PTDHist[2]->SetLineColor(kGreen);

    vector<TH1D*> p8S2Hist;
    p8S2Hist.push_back(new TH1D("p8 sigma2 g","p8 sigma2 g",100,0,0.2));
    p8S2Hist.push_back(new TH1D("p8 sigma2 q","p8 sigma2 q",100,0,0.2));
    p8S2Hist.push_back(new TH1D("p8 sigma2 u","p8 sigma2 u",100,0,0.2));
    for (TH1D* sigma2 : p8S2Hist) { sigma2->Sumw2(); }
    p8S2Hist[0]->SetLineColor(kRed); p8S2Hist[2]->SetLineColor(kGreen);
    
    TH1D* p8Pt = new TH1D("p8 Pt",";pT;Events",ptBins,ptRange);

/* Herwig++ event loop */
    std::size_t hppCount = 0;
    std::size_t hppN = hppTree->GetEntries();
    TProfile* hppWeighting = new TProfile("hppw","",ptBins,ptRange);
    for(size_t x=0; x != hppN; ++x) {
        hppTree->GetEntry(x);
        assert(kMaxfJets>hppJets_);

        for (int i = 0; i < hppJets_; ++i) {
            TLorentzVector tmpVec(hppX[i],hppY[i],hppZ[i],hppT[i]);
            if (fabs(tmpVec.Eta())>etaMax) continue;
            ++hppCount;

            //hppWeight=1;
            hppGluonFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 21)? 1:0, hppWeight);
            hppLightquarkFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 1 || hppFlav[i] == 2)? 1:0, hppWeight);
            hppStrangeFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 3)? 1:0, hppWeight);
            hppCharmFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 4)? 1:0, hppWeight);
            hppBottomFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 5)? 1:0, hppWeight);
            hppUnmatchedFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 0)? 1:0, hppWeight);

            hppPt->Fill(tmpVec.Pt(),hppWeight);
            hppWeighting->Fill(tmpVec.Pt(),hppWeight);

            if(tmpVec.Pt()>ptUpper || tmpVec.Pt()<ptLower) continue;
            hppEtaGluonFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 21)? 1:0, hppWeight);
            hppEtaLightquarkFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 1 || hppFlav[i] == 2)? 1:0, hppWeight);
            hppEtaStrangeFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 3)? 1:0, hppWeight);
            hppEtaCharmFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 4)? 1:0, hppWeight);
            hppEtaBottomFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 5)? 1:0, hppWeight);
            hppEtaUnmatchedFrac.Fill(tmpVec.Eta(), (hppFlav[i] == 0)? 1:0, hppWeight);

            if(hppFlav[i] == 21) {
                hppMultHist[0]->Fill(hppConstituents[i],hppWeight);
                hppPTDHist [0]->Fill(hppPTD[i],hppWeight);
                hppS2Hist  [0]->Fill(hppSigma2[i],hppWeight);
            } else if(hppFlav[i] == 0) {
                hppMultHist[2]->Fill(hppConstituents[i],hppWeight);
                hppPTDHist [2]->Fill(hppPTD[i],hppWeight);
                hppS2Hist  [2]->Fill(hppSigma2[i],hppWeight);
            } else {
                hppMultHist[1]->Fill(hppConstituents[i],hppWeight);
                hppPTDHist [1]->Fill(hppPTD[i],hppWeight);
                hppS2Hist  [1]->Fill(hppSigma2[i],hppWeight);
            }
        }
    }
    cout << hppCount << " herwig jets" << endl;

/* Pythia6 event loop */
    std::size_t p6Count = 0;
    std::size_t p6N  = p6Tree->GetEntries();
    TProfile* p6Weighting = new TProfile("p6w","",ptBins,ptRange);
    for(size_t x=0; x != p6N; ++x) {
        p6Tree->GetEntry(x);
        assert(kMaxfJets>p6Jets_);
        double xDiff = 0;// -( p6X[0]+p6X[1] )/2;
        double yDiff = 0;// -( p6Y[0]+p6Y[1] )/2;

        for (int i = 0; i < p6Jets_; ++i) {
            TLorentzVector tmpVec(p6X[i]+xDiff,p6Y[i]+yDiff,p6Z[i],p6T[i]);
            if (fabs(tmpVec.Eta())>etaMax) continue;
            ++p6Count;
            //p6Weight = 1;

            p6GluonFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 21)? 1:0, p6Weight);
            p6LightquarkFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 1 || p6Flav[i] == 2)? 1:0, p6Weight);
            p6StrangeFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 3)? 1:0, p6Weight);
            p6CharmFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 4)? 1:0, p6Weight);
            p6BottomFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 5)? 1:0, p6Weight);
            p6UnmatchedFrac.Fill(tmpVec.Pt(), (p6Flav[i] == 0)? 1:0, p6Weight);

            p6Pt->Fill(tmpVec.Pt(),p6Weight);
            p6Weighting->Fill(tmpVec.Pt(),p6Weight);

            if(tmpVec.Pt()>ptUpper || tmpVec.Pt()<ptLower) continue;
            p6EtaGluonFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 21)? 1:0, p6Weight);
            p6EtaLightquarkFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 1 || p6Flav[i] == 2)? 1:0, p6Weight);
            p6EtaStrangeFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 3)? 1:0, p6Weight);
            p6EtaCharmFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 4)? 1:0, p6Weight);
            p6EtaBottomFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 5)? 1:0, p6Weight);
            p6EtaUnmatchedFrac.Fill(tmpVec.Eta(), (p6Flav[i] == 0)? 1:0, p6Weight);

            if(p6Flav[i] == 21) {
                p6MultHist[0]->Fill(p6Constituents[i],p6Weight);
                p6PTDHist [0]->Fill(p6PTD[i],p6Weight);
                p6S2Hist  [0]->Fill(p6Sigma2[i],p6Weight);
            } else if(p6Flav[i] == 0) {
                p6MultHist[2]->Fill(p6Constituents[i],p6Weight);
                p6PTDHist [2]->Fill(p6PTD[i],p6Weight);
                p6S2Hist  [2]->Fill(p6Sigma2[i],p6Weight);
            } else {
                p6MultHist[1]->Fill(p6Constituents[i],p6Weight);
                p6PTDHist [1]->Fill(p6PTD[i],p6Weight);
                p6S2Hist  [1]->Fill(p6Sigma2[i],p6Weight);
            }
        }
    }
    cout << p6Count << " pythia6 jets" << endl;

/* Pythia8 event loop */
    std::size_t p8Count = 0;
    std::size_t p8N  = p8Tree->GetEntries();
    TProfile* p8Weighting = new TProfile("p8w","",ptBins,ptRange);
    for(size_t x=0; x != p8N; ++x) {
        p8Tree->GetEntry(x);
        assert(kMaxfJets>p8Jets_);

        for (int i = 0; i < p8Jets_; ++i) {
            TLorentzVector tmpVec(p8X[i],p8Y[i],p8Z[i],p8T[i]);
            if (fabs(tmpVec.Eta())>etaMax) continue;
            ++p8Count;
            //p8Weight = 1;

            p8GluonFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 21)? 1:0, p8Weight);
            p8LightquarkFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 1 || p8Flav[i] == 2)? 1:0, p8Weight);
            p8StrangeFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 3)? 1:0, p8Weight);
            p8CharmFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 4)? 1:0, p8Weight);
            p8BottomFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 5)? 1:0, p8Weight);
            p8UnmatchedFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 0)? 1:0, p8Weight);

            p8Pt->Fill(tmpVec.Pt(),p8Weight);
            p8Weighting->Fill(tmpVec.Pt(),p8Weight);

            if(tmpVec.Pt()>ptUpper || tmpVec.Pt()<ptLower) continue;
            p8EtaGluonFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 21)? 1:0, p8Weight);
            p8EtaLightquarkFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 1 || p8Flav[i] == 2)? 1:0, p8Weight);
            p8EtaStrangeFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 3)? 1:0, p8Weight);
            p8EtaCharmFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 4)? 1:0, p8Weight);
            p8EtaBottomFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 5)? 1:0, p8Weight);
            p8EtaUnmatchedFrac.Fill(tmpVec.Eta(), (p8Flav[i] == 0)? 1:0, p8Weight);

            if(p8Flav[i] == 21) {
                p8MultHist[0]->Fill(p8Constituents[i],p8Weight);
                p8PTDHist [0]->Fill(p8PTD[i],p8Weight);
                p8S2Hist  [0]->Fill(p8Sigma2[i],p8Weight);
            } else if(p8Flav[i] == 0) {
                p8MultHist[2]->Fill(p8Constituents[i],p8Weight);
                p8PTDHist [2]->Fill(p8PTD[i],p8Weight);
                p8S2Hist  [2]->Fill(p8Sigma2[i],p8Weight);
            } else {
                p8MultHist[1]->Fill(p8Constituents[i],p8Weight);
                p8PTDHist [1]->Fill(p8PTD[i],p8Weight);
                p8S2Hist  [1]->Fill(p8Sigma2[i],p8Weight);         
            }
        }
    }
    cout << p8Count << " pythia8 jets" << endl;

/* Herwig++ fraction histograms */
    TH1D *hppLightquarks = hppLightquarkFrac.ProjectionX("hpp light quarks","");
    TH1D *hppGluons = hppGluonFrac.ProjectionX("hpp gluons","");
    TH1D *hppStrange = hppStrangeFrac.ProjectionX("hpp strange","");
    TH1D *hppCharm = hppCharmFrac.ProjectionX("hpp charm","");
    TH1D *hppBottom = hppBottomFrac.ProjectionX("hpp bottom","");
    TH1D *hppUnmatched = hppUnmatchedFrac.ProjectionX("hpp unmatch","");

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

    THStack *hppHs  = new THStack("hpp hs","");
    TCanvas *canv = tdrCanvas("c1",h1,12,0,1);
    hppHs->SetHistogram(h1);
    h1->GetYaxis()->SetTitleOffset(1.25);
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleSize(0.045);
    gPad->SetLogx();

    tdrDraw(hppUnmatched,"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
    tdrDraw(hppGluons,"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(hppLightquarks,"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(hppStrange,"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(hppCharm,"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(hppBottom,"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
    
    //light_quarks->Add(strange);
    hppHs->Add(hppBottom);
    hppHs->Add(hppCharm);
    hppHs->Add(hppStrange);
    hppHs->Add(hppLightquarks);
    hppHs->Add(hppGluons);
    hppHs->Add(hppUnmatched);

/* Pythia6 fraction histograms */    
    TH1D *p6Lightquarks = p6LightquarkFrac.ProjectionX("p6 light quarks","");
    TH1D *p6Gluons = p6GluonFrac.ProjectionX("p6 gluons","");
    TH1D *p6Strange = p6StrangeFrac.ProjectionX("p6 strange","");
    TH1D *p6Charm = p6CharmFrac.ProjectionX("p6 charm","");
    TH1D *p6Bottom = p6BottomFrac.ProjectionX("p6 bottom","");
    TH1D *p6Unmatched = p6UnmatchedFrac.ProjectionX("p6 unmatch","");
    
    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(p6Unmatched,"",kOpenStar,kWhite,kSolid,-1,1001,kGray);
    tdrDraw(p6Gluons,"",kOpenCross,kCyan,kSolid,-1,1001,kBlue-9);
    tdrDraw(p6Lightquarks,"",kOpenDiamond,kOrange+7,kSolid,-1,1001,kYellow-9);
    tdrDraw(p6Strange,"",kOpenSquare,kMagenta,kSolid,-1,1001,kAzure-8);
    tdrDraw(p6Charm,"",kOpenTriangleUp,kSpring,kSolid,-1,1001,kGreen-9);
    tdrDraw(p6Bottom,"",kOpenTriangleDown,kPink,kSolid,-1,1001,kRed-9);

    THStack *p6Hs  = new THStack("p6 hs","");
    p6Hs->SetHistogram(h1);
    
    //light_quarks->Add(strange);
    p6Hs->Add(p6Bottom);
    p6Hs->Add(p6Charm);
    p6Hs->Add(p6Strange);
    p6Hs->Add(p6Lightquarks);
    p6Hs->Add(p6Gluons);
    p6Hs->Add(p6Unmatched);

/* Pythia8 fraction histograms */    
    TH1D *p8Lightquarks = p8LightquarkFrac.ProjectionX("p8 light quarks","");
    TH1D *p8Gluons = p8GluonFrac.ProjectionX("p8 gluons","");
    TH1D *p8Strange = p8StrangeFrac.ProjectionX("p8 strange","");
    TH1D *p8Charm = p8CharmFrac.ProjectionX("p8 charm","");
    TH1D *p8Bottom = p8BottomFrac.ProjectionX("p8 bottom","");
    TH1D *p8Unmatched = p8UnmatchedFrac.ProjectionX("p8 unmatch","");
    
    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(p8Unmatched,"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
    tdrDraw(p8Gluons,"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(p8Lightquarks,"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(p8Strange,"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(p8Charm,"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(p8Bottom,"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);

    THStack *p8Hs  = new THStack("p8 hs","");
    p8Hs->SetHistogram(h1);
    
    //light_quarks->Add(strange);
    p8Hs->Add(p8Bottom);
    p8Hs->Add(p8Charm);
    p8Hs->Add(p8Strange);
    p8Hs->Add(p8Lightquarks);
    p8Hs->Add(p8Gluons);
    p8Hs->Add(p8Unmatched);
    
/* Generic stuff: */
    p6Hs->SetMaximum(0.95);
    p8Hs->SetMaximum(0.95);
    hppHs->SetMaximum(0.95);

/* Fraction plotting */
    h1->GetXaxis()->SetRange(9,47);
    //h1->GetYaxis()->SetRangeUser(-0.001,1.001);
    p8Hs->Draw("");
    hppHs->Draw("sameP");
    p6Hs->Draw("sameP");

/* Fraction legends */
    //heading->AddEntry()
    //TLegend *leg = tdrLeg(0.2,0.41,0.5,0.91);
    TLegend *leg = tdrLeg(0.3,0.4,0.75,0.85);
    TLegend *heading = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
    TLegend *sample = tdrLeg(0.65,0.50+0.05,0.75,0.505+0.05);
    //TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);
    TLegend *etacut = tdrLeg(0.65,0.50,0.75,0.505);

//     gPad->Update();
//     gPad->GetCanvas()->Modified();
  
    if (sampleType==1) {
        sample->SetHeader("dijet sample");
    } else if (sampleType==2) {
        sample->SetHeader("#gamma+jet sample");
    } else if (sampleType==3) {
        sample->SetHeader("Z#mu#mu+jet sample");
    }
    //heading->SetHeader("Pythia8 Simulation (4C Tune)");
    //alphacut->SetHeader("#alpha<0.3");
    etacut->SetHeader("#left|#eta#right|< 2.5");

    hppHs->GetXaxis()->SetNoExponent();
    hppHs->GetXaxis()->SetMoreLogLabels(kTRUE);
    hppHs->GetXaxis()->SetTitle("p_{T} (GeV)");
    hppHs->GetYaxis()->SetTitle("Flavor fraction");
    
    leg->SetNColumns(3);
    //leg->SetHeader("CTEQ6L1/MRST LO**");
    leg->SetHeader("P8 P6 HW");

    leg->AddEntry(p8Unmatched," ","f");
    leg->AddEntry(p6Unmatched," ","p");
    leg->AddEntry(hppUnmatched,"None ","p");
    leg->AddEntry(p8Gluons," ","f");
    leg->AddEntry(p6Gluons," ","p");
    leg->AddEntry(hppGluons,"Gluon","p");
    leg->AddEntry(p8Lightquarks," ","f");
    leg->AddEntry(p6Lightquarks," ","p");
    leg->AddEntry(hppLightquarks,"Light","p");
    leg->AddEntry(p8Strange," ","f");
    leg->AddEntry(p6Strange," ","p");
    leg->AddEntry(hppStrange,"Strange","p");
    leg->AddEntry(p8Charm," ","f");
    leg->AddEntry(p6Charm," ","p");
    leg->AddEntry(hppCharm,"Charm","p");
    leg->AddEntry(p8Bottom," ","f");
    leg->AddEntry(p6Bottom," ","p");
    leg->AddEntry(hppBottom,"Bottom","p");
    //gPad->RedrawAxis();
    //canv->Print("fracs.pdf");

/* Eta plotting */

    // Herwigpp fractions
    TH1D *hppEtaLightquarks = hppEtaLightquarkFrac.ProjectionX("hpp eta light quarks","");
    TH1D *hppEtaGluons = hppEtaGluonFrac.ProjectionX("hpp eta gluons","");
    TH1D *hppEtaStrange = hppEtaStrangeFrac.ProjectionX("hpp eta strange","");
    TH1D *hppEtaCharm = hppEtaCharmFrac.ProjectionX("hpp eta charm","");
    TH1D *hppEtaBottom = hppEtaBottomFrac.ProjectionX("hpp eta bottom","");
    TH1D *hppEtaUnmatched = hppEtaUnmatchedFrac.ProjectionX("hpp eta unmatch","");
    
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

    THStack *hppEtaHs  = new THStack("hpp hs","");
    TCanvas *canvEta = tdrCanvas("c1 eta",h1eta,12,0,1);
    hppEtaHs->SetHistogram(h1eta);
    h1eta->GetYaxis()->SetTitleOffset(1.25);
    h1eta->GetXaxis()->SetTitleOffset(1.0);
    h1eta->GetXaxis()->SetLabelSize(0.045);
    h1eta->GetYaxis()->SetLabelSize(0.045);
    h1eta->GetXaxis()->SetTitleSize(0.045);
    h1eta->GetYaxis()->SetTitleSize(0.045);

    tdrDraw(hppEtaUnmatched,"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
    tdrDraw(hppEtaGluons,"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(hppEtaLightquarks,"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(hppEtaStrange,"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(hppEtaCharm,"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(hppEtaBottom,"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);
    
    //light_quarks->Add(strange);
    hppEtaHs->Add(hppEtaBottom);
    hppEtaHs->Add(hppEtaCharm);
    hppEtaHs->Add(hppEtaStrange);
    hppEtaHs->Add(hppEtaLightquarks);
    hppEtaHs->Add(hppEtaGluons);
    hppEtaHs->Add(hppEtaUnmatched);

/* Pythia6 fraction histograms */
    TH1D *p6EtaLightquarks = p6EtaLightquarkFrac.ProjectionX("p6 eta light quarks","");
    TH1D *p6EtaGluons = p6EtaGluonFrac.ProjectionX("p6 eta gluons","");
    TH1D *p6EtaStrange = p6EtaStrangeFrac.ProjectionX("p6 eta strange","");
    TH1D *p6EtaCharm = p6EtaCharmFrac.ProjectionX("p6 eta charm","");
    TH1D *p6EtaBottom = p6EtaBottomFrac.ProjectionX("p6 eta bottom","");
    TH1D *p6EtaUnmatched = p6EtaUnmatchedFrac.ProjectionX("p6 eta unmatch","");

    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(p6EtaUnmatched,"",kOpenStar,kWhite,kSolid,-1,1001,kGray);
    tdrDraw(p6EtaGluons,"",kOpenCross,kCyan,kSolid,-1,1001,kBlue-9);
    tdrDraw(p6EtaLightquarks,"",kOpenDiamond,kOrange+7,kSolid,-1,1001,kYellow-9);
    tdrDraw(p6EtaStrange,"",kOpenSquare,kMagenta,kSolid,-1,1001,kAzure-8);
    tdrDraw(p6EtaCharm,"",kOpenTriangleUp,kSpring,kSolid,-1,1001,kGreen-9);
    tdrDraw(p6EtaBottom,"",kOpenTriangleDown,kPink,kSolid,-1,1001,kRed-9);

    THStack *p6EtaHs  = new THStack("p6 hs","");
    p6EtaHs->SetHistogram(h1eta);

    //light_quarks->Add(strange);
    p6EtaHs->Add(p6EtaBottom);
    p6EtaHs->Add(p6EtaCharm);
    p6EtaHs->Add(p6EtaStrange);
    p6EtaHs->Add(p6EtaLightquarks);
    p6EtaHs->Add(p6EtaGluons);
    p6EtaHs->Add(p6EtaUnmatched);

/* Pythia8 fraction histograms */
    TH1D *p8EtaLightquarks = p8EtaLightquarkFrac.ProjectionX("p8 eta light quarks","");
    TH1D *p8EtaGluons = p8EtaGluonFrac.ProjectionX("p8 eta gluons","");
    TH1D *p8EtaStrange = p8EtaStrangeFrac.ProjectionX("p8 eta strange","");
    TH1D *p8EtaCharm = p8EtaCharmFrac.ProjectionX("p8 eta charm","");
    TH1D *p8EtaBottom = p8EtaBottomFrac.ProjectionX("p8 eta bottom","");
    TH1D *p8EtaUnmatched = p8EtaUnmatchedFrac.ProjectionX("p8 eta unmatch","");

    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(p8EtaUnmatched,"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
    tdrDraw(p8EtaGluons,"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(p8EtaLightquarks,"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(p8EtaStrange,"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(p8EtaCharm,"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(p8EtaBottom,"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);

    THStack *p8EtaHs  = new THStack("p8 hs","");
    p8EtaHs->SetHistogram(h1eta);

    //light_quarks->Add(strange);
    p8EtaHs->Add(p8EtaBottom);
    p8EtaHs->Add(p8EtaCharm);
    p8EtaHs->Add(p8EtaStrange);
    p8EtaHs->Add(p8EtaLightquarks);
    p8EtaHs->Add(p8EtaGluons);
    p8EtaHs->Add(p8EtaUnmatched);

    p6EtaHs->SetMaximum(0.95);
    p8EtaHs->SetMaximum(0.95);
    hppEtaHs->SetMaximum(0.95);

/* Fraction plotting */
    //h1eta->GetXaxis()->SetRange(9,47);
    //h1->GetYaxis()->SetRangeUser(-0.001,1.001);
    p8EtaHs->Draw("");
    hppEtaHs->Draw("sameP");
    p6EtaHs->Draw("sameP");

/* Fraction legends */
    //heading->AddEntry()
    //TLegend *leg = tdrLeg(0.2,0.41,0.5,0.91);
    TLegend *legEta = tdrLeg(0.3,0.4,0.75,0.85);
    TLegend *headingEta = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
    TLegend *sampleEta = tdrLeg(0.65,0.50+0.05,0.75,0.505+0.05);
    //TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);

//     gPad->Update();
//     gPad->GetCanvas()->Modified();
  
    if (sampleType==1) {
        sampleEta->SetHeader("dijet sample");
    } else if (sampleType==2) {
        sampleEta->SetHeader("#gamma+jet sample");
    } else if (sampleType==3) {
        sampleEta->SetHeader("Z#mu#mu+jet sample");
    }
    //heading->SetHeader("Pythia8 Simulation (4C Tune)");
    //alphacut->SetHeader("#alpha<0.3");

    hppEtaHs->GetXaxis()->SetNoExponent();
    hppEtaHs->GetXaxis()->SetTitle("#eta");
    hppEtaHs->GetYaxis()->SetTitle("Flavor fraction");
    
    legEta->SetNColumns(3);
    //legEta->SetHeader("CTEQ6L1/MRST LO**");
    legEta->SetHeader("P8 P6 Hw");
    
    legEta->AddEntry(p8EtaUnmatched," ","f");
    legEta->AddEntry(p6EtaUnmatched," ","p");
    legEta->AddEntry(hppEtaUnmatched,"None ","p");
    legEta->AddEntry(p8EtaGluons," ","f");
    legEta->AddEntry(p6EtaGluons," ","p");
    legEta->AddEntry(hppEtaGluons,"Gluon","p");
    legEta->AddEntry(p8EtaLightquarks," ","f");
    legEta->AddEntry(p6EtaLightquarks," ","p");
    legEta->AddEntry(hppEtaLightquarks,"Light","p");
    legEta->AddEntry(p8EtaStrange," ","f");
    legEta->AddEntry(p6EtaStrange," ","p");
    legEta->AddEntry(hppEtaStrange,"Strange","p");
    legEta->AddEntry(p8EtaCharm," ","f");
    legEta->AddEntry(p6EtaCharm," ","p");
    legEta->AddEntry(hppEtaCharm,"Charm","p");
    legEta->AddEntry(p8EtaBottom," ","f");
    legEta->AddEntry(p6EtaBottom," ","p");
    legEta->AddEntry(hppEtaBottom,"Bottom","p");
    
    if (indicatorPlots) {
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
        tdrDraw(p8MultHist[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,3003,kRed-3);
        tdrDraw(hppMultHist[0],"HIST P",kFullTriangleDown,kRed-3,kSolid,-1,3003,kRed-3);
        tdrDraw(p8MultHist[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,3005,kBlue);
        tdrDraw(hppMultHist[1],"HIST P",kFullTriangleUp,kBlue,kSolid,-1,3005,kBlue);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        for (TH1D* mult : p8MultHist) { mult->Scale(1/mult->Integral()); }
        for (TH1D* mult : hppMultHist) { mult->Scale(1/mult->Integral()); }

        TLegend *multLeg = tdrLeg(0.5,0.6,0.8,0.9);
        multLeg->SetNColumns(2);

        multLeg->AddEntry(p8MultHist[0],"         ","f");
        multLeg->AddEntry(hppMultHist[0],"       Gluon","p");
        multLeg->AddEntry(p8MultHist[1],"         ","f");
        multLeg->AddEntry(hppMultHist[1],"       Quark","p");
        multLeg->SetHeader("Pythia8/Herwig++");

        /* PTD */
        TH1D *h3 = new TH1D("h3",";p_{T}D;Events",100,0,1);
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        TCanvas *c3 = tdrCanvas("c3",h3,0,33,1);

        h3->SetMaximum(0.075);
        h3->GetYaxis()->SetNoExponent();
        h3->GetXaxis()->SetNoExponent();
        h3->GetXaxis()->SetRangeUser(0,30);

        tdrDraw(p8PTDHist[0],"HIST",kFullTriangleDown,kRed-3,kSolid,-1,3003,kRed-3);
        tdrDraw(hppPTDHist[0],"HIST P",kFullTriangleDown,kRed-3,kSolid,-1,3003,kRed-3);
        tdrDraw(p8PTDHist[1],"HIST",kFullTriangleUp,kBlue,kSolid,-1,3005,kBlue);
        tdrDraw(hppPTDHist[1],"HIST P",kFullTriangleUp,kBlue,kSolid,-1,3005,kBlue);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        //     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);

        for (TH1D* ptd : p8PTDHist) { ptd->Scale(1/ptd->Integral()); }
        for (TH1D* ptd : hppPTDHist) { ptd->Scale(1/ptd->Integral()); }

        TLegend *ptdLeg = tdrLeg(0.5,0.6,0.8,0.9);
        ptdLeg->SetNColumns(2);

        ptdLeg->AddEntry(p8PTDHist[0],"         ","f");
        ptdLeg->AddEntry(hppPTDHist[0],"       Gluon","p");
        ptdLeg->AddEntry(p8PTDHist[1],"         ","f");
        ptdLeg->AddEntry(hppPTDHist[1],"       Quark","p");
        ptdLeg->SetHeader("Pythia8/Herwig++");

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

        tdrDraw(p8S2Hist[0],"HIST",kFullTriangleUp,kRed-3,kSolid,-1,3004,kRed-3);
        tdrDraw(hppS2Hist[0],"HIST P",kFullTriangleUp,kRed-3,kSolid,-1,3004,kRed-3);
        tdrDraw(p8S2Hist[1],"HIST",kFullTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        tdrDraw(hppS2Hist[1],"HIST P",kFullTriangleDown,kBlue,kSolid,-1,3003,kBlue);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);
        //     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);

        for (TH1D* sigma2 : p8S2Hist) { sigma2->Scale(1/sigma2->Integral()); }
        for (TH1D* sigma2 : hppS2Hist) { sigma2->Scale(1/sigma2->Integral()); }

        TLegend *s2Leg = tdrLeg(0.5,0.6,0.8,0.9);
        s2Leg->SetNColumns(2);

        s2Leg->AddEntry(p8S2Hist[0],"         ","f");
        s2Leg->AddEntry(hppS2Hist[0],"       Gluon","p");
        s2Leg->AddEntry(p8S2Hist[1],"         ","f");
        s2Leg->AddEntry(hppS2Hist[1],"       Quark","p");
        s2Leg->SetHeader("Pythia8/Herwig++");
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
}
