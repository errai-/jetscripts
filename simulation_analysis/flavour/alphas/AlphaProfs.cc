#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>

using namespace std;

void AlphaProfs(Int_t generator = 0)
{
	TDirectory *curdir = gDirectory;

	string files[3];

    files[0] = "../cteq6l1-cms/pythia8_dijet_physics_2000000.root";
    files[1] = "../cteq6l1-cms/pythia6_dijet_physics_1000000.root";
    files[2] = "../cteq6l1-cms/herwig_dijet_physics_2000000.root";

    Int_t ptBins = 40;
    Double_t ptRange[] = {40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840};

    vector<TProfile*> profs;
    vector<TH1D*> weights;
    profs.push_back(new TProfile("prof10","",ptBins,ptRange));
    weights.push_back(new TH1D("weights10","",ptBins,ptRange));
    profs.push_back(new TProfile("prof15","",ptBins,ptRange));
    weights.push_back(new TH1D("weights15","",ptBins,ptRange));
    profs.push_back(new TProfile("prof20","",ptBins,ptRange));
    weights.push_back(new TH1D("weights20","",ptBins,ptRange));
    profs.push_back(new TProfile("prof30","",ptBins,ptRange));
    weights.push_back(new TH1D("weights30","",ptBins,ptRange));

	stringstream tmpString("");
	tmpString << files[generator];
	TFile *f = new TFile(tmpString.str().c_str(),"READ");
	assert(f && !f->IsZombie());
	
	TTree *tree = (TTree*)f->Get("JetTree");
	assert(tree && !tree->IsZombie());

    static const Int_t kMaxfJets = 100;

    Int_t           fJets_;
    Double_t        fX[kMaxfJets];   //[fJets_]
    Double_t        fY[kMaxfJets];   //[fJets_]
    Double_t        fZ[kMaxfJets];   //[fJets_]
    Double_t        fT[kMaxfJets];   //[fJets_]
    Double_t        fWeight;
    Int_t           fFlav[kMaxfJets];   //[fJets_]
    Int_t           fConstituents[kMaxfJets];   //[fJets_]
    Double_t        fPTD[kMaxfJets];   //[fJets_]
    Double_t        fSigma2[kMaxfJets];   //[fJets_]
    Double_t        fAlpha[kMaxfJets];   //[fJets_]
    Double_t        fPT[kMaxfJets];   //[fJets_]

    tree->SetMakeClass(1);
	unsigned int N = (unsigned int)tree->GetEntries(); 
	
    tree->SetBranchAddress("fJets", &fJets_);
    tree->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
    tree->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
    tree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
    tree->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
    tree->SetBranchAddress("fWeight", &fWeight);
    tree->SetBranchAddress("fJets.fAlpha", fAlpha);
    tree->SetBranchAddress("fJets.fPartonPT", fPT);

	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);

        for (int i = 0; i < fJets_; ++i) {
            TLorentzVector p4(fX[i],fY[i],fZ[i],fT[i]);

		    if(fPT[i]>840 || fPT[i]<40) continue;
            if(abs(p4.Eta())>1.3) continue;

            if (fAlpha[i] < 0.1) {
                profs[0]->Fill(fPT[i],p4.Pt()/fPT[i],fWeight);
                weights[0]->Fill(fPT[i],fWeight);
            }
            if (fAlpha[i] < 0.15) {
                profs[1]->Fill(fPT[i],p4.Pt()/fPT[i],fWeight);
                weights[1]->Fill(fPT[i],fWeight);
            }
            if (fAlpha[i] < 0.2) {
                profs[2]->Fill(fPT[i],p4.Pt()/fPT[i],fWeight);
                weights[2]->Fill(fPT[i],fWeight);
            }
            if (fAlpha[i] < 0.3) {
                profs[3]->Fill(fPT[i],p4.Pt()/fPT[i],fWeight);
                weights[3]->Fill(fPT[i],fWeight);
            }
        }
	}
    f->Close();	
	
	stringstream tmpString2;
    tmpString2 << "alphafracs_";
    if (generator==0) {
        tmpString2 << "p8";
    } else if (generator==1) {
        tmpString2 << "p6";
    } else {
        tmpString2 << "hwpp";
    }
    tmpString2 << ".root";
    TFile *g = new TFile(tmpString2.str().c_str(),"RECREATE");
	assert(g && !g->IsZombie());

    for (std::size_t i = 0; i < 4; ++i) {
        profs[i]->Write();
    }

    for (std::size_t i = 0; i < 4; ++i) {
        weights[i]->Scale(1./weights[i]->Integral());
        weights[i]->Write();
    }
}

