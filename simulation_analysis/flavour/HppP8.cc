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

#include "tdrstyle_mod14.C"

using std::cout;
using std::endl;
using std::string;

const int ptBins = 48.;
const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

void HppP8(string herwFile, string pythFile) {
    /* Init branches */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    
    static const Int_t kMaxfJets = 100;
    
    Int_t           hppJets_;
    Double_t        hppX[kMaxfJets];   //[hppJets_]
    Double_t        hppY[kMaxfJets];   //[hppJets_]
    Double_t        hppZ[kMaxfJets];   //[hppJets_]
    Double_t        hppT[kMaxfJets];   //[hppJets_]
    Int_t           p8Jets_;
    Double_t        p8X[kMaxfJets];   //[p8Jets_]
    Double_t        p8Y[kMaxfJets];   //[p8Jets_]
    Double_t        p8Z[kMaxfJets];   //[p8Jets_]
    Double_t        p8T[kMaxfJets];   //[p8Jets_]
    
    Double_t        hppWeight;
    Double_t        p8Weight;
    Int_t           hppFlav[kMaxfJets];   //[hppJets_]
    Int_t           p8Flav[kMaxfJets];   //[p8Jets_]
    
    Int_t           hppConstituents[kMaxfJets];   //[hppJets_]
    Double_t        hppPTD[kMaxfJets];   //[hppJets_]
    Double_t        hppSigma2[kMaxfJets];   //[hppJets_]
    Int_t           p8Constituents[kMaxfJets];   //[p8Jets_]
    Double_t        p8PTD[kMaxfJets];   //[p8Jets_]
    Double_t        p8Sigma2[kMaxfJets];   //[p8Jets_]
    
    /* Tree setup */
    TTree* hppTree; TTree* p8Tree;
    
    TChain * hppChain = new TChain("JetTree","");
    hppChain->Add(herwFile.c_str()); hppTree = hppChain;
    TChain * p8Chain = new TChain("JetTree","");
    p8Chain->Add(pythFile.c_str()); p8Tree = p8Chain;
    
    if (!p8Tree || !hppTree) return;
    hppTree->SetMakeClass(1); p8Tree->SetMakeClass(1);
    

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
    TProfile hppGluonFrac("g","g",ptBins,ptRange);
    TProfile hppLightquarkFrac("lq","lq",ptBins,ptRange);
    TProfile hppStrangeFrac("s","s",ptBins,ptRange);
    TProfile hppCharmFrac("c","c",ptBins,ptRange);
    TProfile hppBottomFrac("b","b",ptBins,ptRange);
    TProfile hppUnmatchedFrac("unmatched","unmatched",ptBins,ptRange);
    
    TH1D* multiplicity_g = new TH1D("multiplicity_g","multiplicity_g",60,0,60);
    TH1D* multiplicity_q = new TH1D("multiplicity_q","multiplicity_q",60,0,60);
    TH1D* multiplicity_u = new TH1D("multiplicity_u","multiplicity_u",60,0,60);
    multiplicity_q->Sumw2();
    multiplicity_g->Sumw2();
    multiplicity_u->Sumw2();
    multiplicity_g->SetLineColor(kRed);
    multiplicity_u->SetLineColor(kGreen);

    TH1D* pTD_g = new TH1D("pTD_g","pTD_g",100,0,1);
    TH1D* pTD_q = new TH1D("pTD_q","pTD_q",100,0,1);
    TH1D* pTD_u = new TH1D("pTD_u","pTD_u",100,0,1);
    pTD_u->Sumw2();
    pTD_g->Sumw2();
    pTD_q->Sumw2();
    pTD_g->SetLineColor(kRed);
    pTD_u->SetLineColor(kGreen);

    TH1D* sigma2_g = new TH1D("sigma2_g","sigma2_g",100,0,0.2);
    TH1D* sigma2_q = new TH1D("sigma2_q","sigma2_q",100,0,0.2);
    TH1D* sigma2_u = new TH1D("sigma2_u","sigma2_u",100,0,0.2);
    sigma2_u->Sumw2();
    sigma2_g->Sumw2();
    sigma2_q->Sumw2();
    sigma2_g->SetLineColor(kRed);
    sigma2_u->SetLineColor(kGreen);
     
    
/* Pythia8 plot initialization */
    TProfile p8GluonFrac("g","g",ptBins,ptRange);
    TProfile p8LightquarkFrac("lq","lq",ptBins,ptRange);
    TProfile p8StrangeFrac("s","s",ptBins,ptRange);
    TProfile p8CharmFrac("c","c",ptBins,ptRange);
    TProfile p8BottomFrac("b","b",ptBins,ptRange);
    TProfile p8UnmatchedFrac("unmatched","unmatched",ptBins,ptRange);
    
    TH1D* p8multiplicity_g = new TH1D("multiplicity_g","multiplicity_g",60,0,60);
    TH1D* p8multiplicity_q = new TH1D("multiplicity_q","multiplicity_q",60,0,60);
    TH1D* p8multiplicity_u = new TH1D("multiplicity_u","multiplicity_u",60,0,60);
    p8multiplicity_q->Sumw2();
    p8multiplicity_g->Sumw2();
    p8multiplicity_u->Sumw2();
    p8multiplicity_g->SetLineColor(kRed);
    p8multiplicity_u->SetLineColor(kGreen);

    TH1D* p8pTD_g = new TH1D("pTD_g","pTD_g",100,0,1);
    TH1D* p8pTD_q = new TH1D("pTD_q","pTD_q",100,0,1);
    TH1D* p8pTD_u = new TH1D("pTD_u","pTD_u",100,0,1);
    p8pTD_u->Sumw2();
    p8pTD_g->Sumw2();
    p8pTD_q->Sumw2();
    p8pTD_g->SetLineColor(kRed);
    p8pTD_u->SetLineColor(kGreen);

    TH1D* p8sigma2_g = new TH1D("sigma2_g","sigma2_g",100,0,0.2);
    TH1D* p8sigma2_q = new TH1D("sigma2_q","sigma2_q",100,0,0.2);
    TH1D* p8sigma2_u = new TH1D("sigma2_u","sigma2_u",100,0,0.2);
    p8sigma2_u->Sumw2();
    p8sigma2_g->Sumw2();
    p8sigma2_q->Sumw2();
    p8sigma2_g->SetLineColor(kRed);
    p8sigma2_u->SetLineColor(kGreen);
     

    std::size_t hppN = hppTree->GetEntries();
    std::size_t p8N  = p8Tree->GetEntries();
    
    for(size_t x=0; x != hppN; ++x) {
        hppTree->GetEntry(x);
        assert(kMaxfJets>hppJets_);
        
        for (int i = 0; i < hppJets_; ++i) {
            TLorentzVector tmpVec(hppX[i],hppY[i],hppZ[i],hppT[i]);
            
            hppGluonFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 21)? 1:0);
            hppLightquarkFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 1 || hppFlav[i] == 2)? 1:0);
            hppStrangeFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 3)? 1:0);
            hppCharmFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 4)? 1:0);
            hppBottomFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 5)? 1:0);
            hppUnmatchedFrac.Fill(tmpVec.Pt(), (hppFlav[i] == 0)? 1:0);
                                
            if(tmpVec.Pt()>100 || tmpVec.Pt()<80) continue;
            if(hppFlav[i] == 21) {
                multiplicity_g->Fill(hppConstituents[i],hppWeight);
                pTD_g->Fill(hppPTD[i],hppWeight);
                sigma2_g->Fill(hppSigma2[i],hppWeight);
            } else if(hppFlav[i] == 0) {
                multiplicity_u->Fill(hppConstituents[i],hppWeight);
                pTD_u->Fill(hppPTD[i],hppWeight);
                sigma2_u->Fill(hppSigma2[i],hppWeight);
            } else {
                multiplicity_q->Fill(hppConstituents[i],hppWeight);
                pTD_q->Fill(hppPTD[i],hppWeight);
                sigma2_q->Fill(hppSigma2[i],hppWeight);         
            }
        }
    }
    
    for(size_t x=0; x != p8N; ++x) {
        p8Tree->GetEntry(x);
        assert(kMaxfJets>p8Jets_);
        
        for (int i = 0; i < p8Jets_; ++i) {
            TLorentzVector tmpVec(p8X[i],p8Y[i],p8Z[i],p8T[i]);
            
            p8GluonFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 21)? 1:0);
            p8LightquarkFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 1 || p8Flav[i] == 2)? 1:0);
            p8StrangeFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 3)? 1:0);
            p8CharmFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 4)? 1:0);
            p8BottomFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 5)? 1:0);
            p8UnmatchedFrac.Fill(tmpVec.Pt(), (p8Flav[i] == 0)? 1:0);
                                
            if(tmpVec.Pt()>100 || tmpVec.Pt()<80) continue;
            if(p8Flav[i] == 21) {
                p8multiplicity_g->Fill(p8Constituents[i],p8Weight);
                p8pTD_g->Fill(p8PTD[i],p8Weight);
                p8sigma2_g->Fill(p8Sigma2[i],p8Weight);
            } else if(p8Flav[i] == 0) {
                p8multiplicity_u->Fill(p8Constituents[i],p8Weight);
                p8pTD_u->Fill(p8PTD[i],p8Weight);
                p8sigma2_u->Fill(p8Sigma2[i],p8Weight);
            } else {
                p8multiplicity_q->Fill(p8Constituents[i],p8Weight);
                p8pTD_q->Fill(p8PTD[i],p8Weight);
                p8sigma2_q->Fill(p8Sigma2[i],p8Weight);         
            }
        }
    }
    
    TH1D *light_quarks = hppLightquarkFrac.ProjectionX("light quarks","");
    TH1D *gluons = hppGluonFrac.ProjectionX("gluons","");
    TH1D *strange = hppStrangeFrac.ProjectionX("strange","");
    TH1D *charm = hppCharmFrac.ProjectionX("charm","");
    TH1D *bottom = hppBottomFrac.ProjectionX("bottom","");
    TH1D *unmatched = hppUnmatchedFrac.ProjectionX("unmatched","");
    
    TH1D *h1 = new TH1D("h",";p_{T} (GeV);Fraction",ptBins,ptRange);
    gROOT->ForceStyle();
    gStyle->SetOptStat(kFALSE); //removes old legend

    THStack *hs  = new THStack("hs","test stacked histograms");
    TCanvas *canv = tdrCanvas("c1",h1,12,0,1);
    hs->SetHistogram(h1);
    setTDRStyle();
    gStyle->SetOptLogx(1);
    h1->GetYaxis()->SetTitleOffset(1.25);
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleSize(0.045);
    gPad->SetLogx();

    tdrDraw(unmatched,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
    tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(strange,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);
    
    //light_quarks->Add(strange);
    hs->Add(bottom);
    hs->Add(charm);
    hs->Add(strange);
    hs->Add(light_quarks);
    hs->Add(gluons);
    hs->Add(unmatched);
    
    TH1D *p8light_quarks = p8LightquarkFrac.ProjectionX("p8 light quarks","");
    TH1D *p8gluons = p8GluonFrac.ProjectionX("p8 gluons","");
    TH1D *p8strange = p8StrangeFrac.ProjectionX("p8 strange","");
    TH1D *p8charm = p8CharmFrac.ProjectionX("p8 charm","");
    TH1D *p8bottom = p8BottomFrac.ProjectionX("p8 bottom","");
    TH1D *p8unmatched = p8UnmatchedFrac.ProjectionX("p8 unmatched","");
    
    TH1D *p8h = new TH1D("p8h",";p_{T} (GeV);Fraction",1000,20,2000);

    tdrDraw(p8unmatched,"",kFullStar,kGray+2,kSolid,-1,1001,kGray);
    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(p8gluons,"",kFullDotLarge,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(p8light_quarks,"",kFullDiamond,kYellow-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(p8strange,"",kFullSquare,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(p8charm,"",kFullTriangleUp,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(p8bottom,"",kFullTriangleDown,kRed-2,kSolid,-1,1001,kRed-9);

    THStack *p8hs  = new THStack("hs","test stacked histograms");
    p8hs->SetHistogram(h1);
    
    //light_quarks->Add(strange);
    p8hs->Add(p8bottom);
    p8hs->Add(p8charm);
    p8hs->Add(p8strange);
    p8hs->Add(p8light_quarks);
    p8hs->Add(p8gluons);
    p8hs->Add(p8unmatched);

    //heading->AddEntry()
    TLegend *leg = tdrLeg(0.5,0.82,0.175,0.50);
    TLegend *heading = tdrLeg(0.675-0.3,0.50+0.44,0.775-0.3,0.505+0.44);
    TLegend *sample = tdrLeg(0.675-0.05,0.50+0.05,0.775-0.05,0.505+0.05);
    TLegend *alphacut = tdrLeg(0.77,0.50,0.87,0.505);
    TLegend *etacut = tdrLeg(0.61,0.50,0.71,0.505);

//     gPad->Update();
//     gPad->GetCanvas()->Modified();
  
    sample->SetHeader("#gamma+jet sample");
    heading->SetHeader("Pythia8 Simulation (4C Tune)");
    alphacut->SetHeader("#alpha<0.3");
    etacut->SetHeader("#left|#eta#right|< 1.3,");

    hs->GetXaxis()->SetNoExponent();
    hs->GetXaxis()->SetMoreLogLabels(kTRUE);
    hs->GetXaxis()->SetTitle("p_{T} (GeV)");
    hs->GetYaxis()->SetTitle("Flavor fraction");
    
    leg->AddEntry(unmatched,"None","f");
    leg->AddEntry(gluons,"Gluon","f");
    leg->AddEntry(light_quarks,"Light","f");
    leg->AddEntry(strange,"Strange","f");
    leg->AddEntry(charm,"Charm","f");
    leg->AddEntry(bottom,"Bottom","f");
    
    hs->Draw("same");
    p8hs->Draw("sameP");
    leg->Draw("same");
    fixOverlay();
    hs->GetXaxis()->SetRange(9,47);
    hs->GetYaxis()->SetRangeUser(-0.01,1.01);

    //STACKED
    TH1D *h2 = new TH1D("h2",";p_{T}D;Events",100,0,1);
    gStyle->SetOptLogx(0);
    TCanvas *c3 = tdrCanvas("c3",h2,0,33);

    TH1F *uD = (TH1F*)pTD_u->Clone("uD");
    TH1F *gD = (TH1F*)pTD_g->Clone("gD");
    TH1F *qD = (TH1F*)pTD_q->Clone("qD");

    h2->SetMaximum(0.07);
    h2->GetYaxis()->SetNoExponent();
    h2->GetXaxis()->SetNoExponent();
    h2->GetXaxis()->SetRangeUser(0,30);
//     pTD_g->Add(pTD_u);
//     pTD_q->Add(pTD_g);
//     p8pTD_g->Add(p8pTD_u);
//     p8pTD_q->Add(p8pTD_g);

    tdrDraw(pTD_g,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
//     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    tdrDraw(pTD_q,"HIST",kDiamond,kBlue,kSolid,-1,3005,kBlue);
    tdrDraw(p8pTD_g,"P",kDiamond,kRed-3,kSolid,-1,3003,kRed-3);
//     tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    tdrDraw(p8pTD_q,"P",kCircle,kBlue,kSolid,-1,3005,kBlue);
    pTD_g->Scale(1/pTD_g->Integral());
    pTD_q->Scale(1/pTD_q->Integral());
    p8pTD_g->Scale(1/p8pTD_g->Integral());
    p8pTD_q->Scale(1/p8pTD_q->Integral());
//     h2->GetXaxis()->SetRangeUser(0.1,1);

    // Multiplicity
    TH1D *h4 = new TH1D("h4",";Number of constituents;Events",60,0,60);
    gStyle->SetOptLogx(0);
    TCanvas *c4 = tdrCanvas("c4",h4,0,33);

    h4->SetMinimum(0);
    h4->SetMaximum(0.07);
    h4->GetYaxis()->SetNoExponent();
    h4->GetXaxis()->SetNoExponent();
    h4->GetXaxis()->SetRangeUser(0,60);
//     multiplicity_g->Add(multiplicity_u);
//     multiplicity_q->Add(multiplicity_g);

    tdrDraw(multiplicity_g,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
//     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    tdrDraw(multiplicity_q,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);

    tdrDraw(p8multiplicity_g,"P",kDot,kRed-3,kSolid,-1,3003,kRed-3);
//     tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    tdrDraw(p8multiplicity_q,"P",kDot,kBlue,kSolid,-1,3005,kBlue);
    multiplicity_g->Scale(1/multiplicity_g->Integral());
    multiplicity_q->Scale(1/multiplicity_q->Integral());
    p8multiplicity_g->Scale(1/p8multiplicity_g->Integral());
    p8multiplicity_q->Scale(1/p8multiplicity_q->Integral());
   
    // Sigma2
    TH1D *h5 = new TH1D("h5",";#sigma_{2};Events",100,0,0.2);
    gStyle->SetOptLogx(0);
    TCanvas *c5 = tdrCanvas("c5",h5,0,33);
    h5->SetMinimum(0);
    h5->SetMaximum(0.06);
    h5->GetYaxis()->SetNoExponent();
    h5->GetXaxis()->SetNoExponent();
    h5->GetXaxis()->SetRangeUser(0,0.2);

    tdrDraw(sigma2_q,"HIST",kDot,kBlue,kSolid,-1,3003,kBlue);
    tdrDraw(sigma2_g,"HIST",kDot,kRed-3,kSolid,-1,3004,kRed-3);
//     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);

    tdrDraw(p8sigma2_q,"P",kDot,kBlue,kSolid,-1,3003,kBlue);
    tdrDraw(p8sigma2_g,"P",kDot,kRed-3,kSolid,-1,3004,kRed-3);
//     tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);
    sigma2_q->Scale(1/sigma2_q->Integral());
    sigma2_g->Scale(1/sigma2_g->Integral());
//                 gs->Scale(1/gs->Integral());
    p8sigma2_q->Scale(1/p8sigma2_q->Integral());
    p8sigma2_g->Scale(1/p8sigma2_g->Integral());
}