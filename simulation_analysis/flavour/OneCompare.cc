#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>

using namespace std;

void OneCompare(Int_t generator = 0)
{
	TDirectory *curdir = gDirectory;

	string files[3];
    if (generator == 0) {
        files[0] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia8_Zjet_physics_1000000.root";
        files[1] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia8_gammajet_physics_1000000.root";
        files[2] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia8_dijet_physics_1000000.root";
    } else if (generator == 1) {
        files[0] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/herwig_Zjet_physics_1000000.root";
        files[1] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/herwig_gammajet_physics_1000000.root";
        files[2] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/herwig_dijet_physics_1000000.root";
    } else if (generator == 2) {
        files[0] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia6_Zjet_physics_1000000.root";
        files[1] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia6_gammajet_physics_1000000.root";
        files[2] = "/home/hannu/Cern/jet_analysis/simulation_analysis/flavour/sim_dir/cteq6l1/pythia6_dijet_physics_1000000.root";
    }
    int markers[] = {kFullSquare,kFullCircle,kFullTriangleUp,kOpenSquare,kOpenCircle,kOpenTriangleUp};
    int color[] = {kRed-3,kBlue,kGreen-6};
    string variable[] = {"sigma2","pTD","multiplicity"};
    string partons[] = {"gluon","quark","unmatched"};
    float x_min[] = {0, 0, 0};
    float x_max[] = {0.2, 1, 60};
    float y_min[] = {0, 0, 0};
    float y_max[] = {0.08, 0.08, 0.1};
    int bins[] = {100, 100, 60};

    vector<vector<vector<TH1D*> > > plots;
	for(unsigned int q=0; q != 3; ++q) { //loop over files
		std::stringstream tmpString("");
		tmpString << files[q];
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

        tree->SetMakeClass(1);
		unsigned int N = (unsigned int)tree->GetEntries(); 
		
        tree->SetBranchAddress("fJets", &fJets_);
        tree->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
        tree->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
        tree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
        tree->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
        tree->SetBranchAddress("fWeight", &fWeight);
        tree->SetBranchAddress("fJets.fFlav", fFlav);
        tree->SetBranchAddress("fJets.fConstituents", fConstituents);
        tree->SetBranchAddress("fJets.fPTD", fPTD);
        tree->SetBranchAddress("fJets.fSigma2", fSigma2);

		vector<vector<TH1D*> > one_sample;
		for(unsigned int p = 0; p < 3; ++p) { //book histograms for gluons, quarks, and unmatched
			vector<TH1D*> s;
			for(unsigned int k=0; k < 3; ++k) { //QGL variables
				std::stringstream temp("");
				temp << variable[k] << "_" << partons[p];
				s.push_back(new TH1D(temp.str().c_str(),temp.str().c_str(),bins[k],x_min[k],x_max[k]));
			}
			one_sample.push_back(s);
		}

		plots.push_back(one_sample);
		for(unsigned int x=0; x != N; ++x)
		{
			tree->GetEntry(x);

            for (int i = 0; i < fJets_; ++i) {
                TLorentzVector p4(fX[i],fY[i],fZ[i],fT[i]);

			    if(p4.Pt()>100 || p4.Pt()<80) continue;
			    if(fFlav[i] == 21) {
				    one_sample[0][0]->Fill(fSigma2[i],fWeight);
				    one_sample[0][1]->Fill(fPTD[i],fWeight);
				    one_sample[0][2]->Fill(fConstituents[i],fWeight);
			    } else if(fFlav[i] <= 0) {
				    one_sample[2][0]->Fill(fSigma2[i],fWeight);
				    one_sample[2][1]->Fill(fPTD[i],fWeight);
				    one_sample[2][2]->Fill(fConstituents[i],fWeight);
			    } else {
				    one_sample[1][0]->Fill(fSigma2[i],fWeight);
				    one_sample[1][1]->Fill(fPTD[i],fWeight);
				    one_sample[1][2]->Fill(fConstituents[i],fWeight);
			    }
            }
		}
	}
	
	setTDRStyle();
	
	vector<TH1D*> xkcd;
	TH1D *h1 = new TH1D("h1",";#sigma_{2};Events",bins[0],x_min[0],x_max[0]);
	TH1D *h2 = new TH1D("h2",";p_{T}D;Events",bins[1],x_min[1],x_max[1]);
	TH1D *h3 = new TH1D("h3",";multiplicity;Events",bins[2],x_min[2],x_max[2]);
	xkcd.push_back(h1);
	xkcd.push_back(h2);
	xkcd.push_back(h3);
		
	for(unsigned int i = 0; i < 3; ++i)
	{
		xkcd[i] ->SetMinimum(y_min[i]);
		xkcd[i] ->SetMaximum(y_max[i]);
		xkcd[i] ->GetYaxis()->SetNoExponent();
		xkcd[i] ->GetXaxis()->SetNoExponent();
		xkcd[i] ->GetXaxis()->SetRangeUser(x_min[i],x_max[i]);
	}

	vector<TCanvas*> c(3);

	for(int l = 0; l<3 ; ++l)
	{
		std::stringstream c_head("");

		for(int t = 0; t < 3; ++t)
		{
			if(t%3==0) 
			{
				c_head << variable[l];
				c[t] = tdrCanvas(c_head.str().c_str(),xkcd[l],0,33);
				c_head.str("");
			}
			tdrDraw(plots[t][0][l] ,"P",markers[t%3] ,color[t%3]);
			tdrDraw(plots[t][1][l],"P",markers[3+t%3],color[t%3]);
			plots[t][0][l]->Scale(1/plots[t][0][l]->Integral());
			plots[t][1][l]->Scale(1/plots[t][1][l]->Integral());
			
			if(t%3==2)
			{
				TLegend* leg = tdrLeg(0.88+0.05,0.72,0.675+0.05,0.30+0.045*4);
				leg->SetHeader("gluon");
				leg->AddEntry(plots[t-2][0][l],"Z+jet","P");
				leg->AddEntry(plots[t-1][0][l],"#gamma+jet","P");
				leg->AddEntry(plots[t][0][l],"dijet","P");
				leg->Draw();

				TLegend* gel = tdrLeg(0.88+0.05,0.48,0.675+0.05,0.24);
				gel->SetHeader("quark");
				gel->AddEntry(plots[t-2][1][l],"Z+jet","P");
				gel->AddEntry(plots[t-1][1][l],"#gamma+jet","P");
				gel->AddEntry(plots[t][1][l],"dijet","P");
				gel->Draw();				
			}
		}
	}
}

