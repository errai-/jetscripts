#include "common_ttbar.h"

void FetchData(string fileName, vector<TH1F*>& mult, vector<TH1F*>& ptd, vector<TH1F*>& s2, unsigned idx) {
    unsigned place = 3*idx;
	mult.push_back( new TH1F("","",60,0,60) );
	mult.push_back( new TH1F("","",60,0,60) );
	mult.push_back( new TH1F("","",60,0,60) );
	ptd.push_back( new TH1F("","",100,0,1) );
	ptd.push_back( new TH1F("","",100,0,1) );
	ptd.push_back( new TH1F("","",100,0,1) );
	s2.push_back( new TH1F("","",100,0,0.2) );
	s2.push_back( new TH1F("","",100,0,0.2) );
	s2.push_back( new TH1F("","",100,0,0.2) );

    static const Int_t kMaxfJets = 100;

    Int_t           fJets;
    Double_t        fX[kMaxfJets];   //[fJets]
    Double_t        fY[kMaxfJets];   //[fJets]
    Double_t        fZ[kMaxfJets];   //[fJets]
    Double_t        fT[kMaxfJets];   //[fJets]

    Double_t        fWeight;
    Int_t          fFlav[kMaxfJets];   //[fJets]
    Int_t          fConstituents[kMaxfJets];   //[fJets_]
    Float_t        fPTD[kMaxfJets];   //[fJets_]
    Float_t        fSigma2[kMaxfJets];   //[fJets_]

    /* Tree setup */
    TTree* jetTree;

    TChain* jetChain = new TChain("JetTree","");
    jetChain->Add(fileName.c_str()); jetTree = jetChain;

    jetTree->SetMakeClass(1);

    jetTree->SetBranchAddress("fJets", &fJets);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
    jetTree->SetBranchAddress("fWeight", &fWeight);
    jetTree->SetBranchAddress("fJets.fFlav", fFlav);
    jetTree->SetBranchAddress("fJets.fConstituents", fConstituents);
    jetTree->SetBranchAddress("fJets.fPTD", fPTD);
    jetTree->SetBranchAddress("fJets.fSigma2", fSigma2);

    unsigned N = jetTree->GetEntries();

    for (unsigned i = 0; i < 3; ++i) {
        mult[place+i]->Sumw2();
        ptd[place+i]->Sumw2();
        s2[place+i]->Sumw2();
    }
	mult[place]->SetLineColor(kRed);
	mult[place+2]->SetLineColor(kGreen);
	ptd[place]->SetLineColor(kRed);
	ptd[place+2]->SetLineColor(kGreen);
	s2[place]->SetLineColor(kRed);
	s2[place+2]->SetLineColor(kGreen);

    unsigned acceptance = 0;
    for(size_t x=0; x != N; ++x)
    {
        jetTree->GetEntry(x);
        assert(kMaxfJets>fJets);

        if (idx == 100) {
            vector<unsigned> flavours;
            vector<TLorentzVector> bjets, ljets;
            TLorentzVector MET, lepton, neutrino;
            unsigned flav_count = 0;
            for (int i = 0; i < fJets; ++i) {
                int fl = fFlav[i];
                TLorentzVector tmpVec(fX[i],fY[i],fZ[i],fT[i]);
                if (fl==10)
                    MET = tmpVec;
                else if (fl==11||fl==13||fl==15)
                    lepton = tmpVec;
                else if (fl==12||fl==14||fl==16)
                    neutrino = tmpVec;
                else {
                    if (tmpVec.Pt() < 30)
                        continue;

                    if (fl==5)
                        bjets.push_back(tmpVec);
                    else {
                        if (fl!=0)
                            ++flav_count;
                        ljets.push_back(tmpVec);
                        flavours.push_back(fl);
                    }
                }
            }
            if (bjets.size()!=2) {
                continue;
            }
            if (flav_count<2) {
                continue;
            }

            // Containers for combined lorentz vectors
            TLorentzVector t1, t1_alt, t2, t3, t4, t5, t6; 

            // Reconstruct W mass from MET and the lepton
            t1 = neutrino + lepton;
            if (t1.M() < 60 || t1.M() > 110) {
                continue;
            }

            // Find jet pairs that correspond to the W
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
            if (working.size() == 0) {
                continue;
            }
        }

        for (int i = 0; i < fJets; ++i) {
            int id = fFlav[i];
            if ( id>=10 && id<20 )
                continue;
            if ( id==0 )
                continue;
            TLorentzVector tmpVec(fX[i],fY[i],fZ[i],fT[i]);
		
            if (fabs(tmpVec.Eta())>1.3)
                continue;

            if(tmpVec.Pt()>100 || tmpVec.Pt()<80)
                continue;

            ++acceptance;
		    if(id == 21) {
			    mult[place]->Fill(fConstituents[i],fWeight);
			    ptd[place]->Fill(fPTD[i],fWeight);
			    s2[place]->Fill(fSigma2[i],fWeight);
		    } else if(id == 5) {
			    mult[place+2]->Fill(fConstituents[i],fWeight);
			    ptd[place+2]->Fill(fPTD[i],fWeight);
			    s2[place+2]->Fill(fSigma2[i],fWeight);
		    } else if (id > 0) {
			    mult[place+1]->Fill(fConstituents[i],fWeight);
			    ptd[place+1]->Fill(fPTD[i],fWeight);
			    s2[place+1]->Fill(fSigma2[i],fWeight);	
		    }
        }
    }
    cout << "Accepted: " << acceptance << endl;
}

void Compare(string ttbarFile, string diFile)
{
    /* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");

    vector<TH1F*> mult;
    vector<TH1F*> ptd;
    vector<TH1F*> s2;
    FetchData(ttbarFile, mult, ptd, s2, 0);
    FetchData(diFile, mult, ptd, s2, 1);

	/*****************Constituents*****************/

	setTDRStyle();
	TH1D *h1 = new TH1D("h1",";Number of constituents;Events",60,0,60);
	h1->SetMinimum(0);
	h1->SetMaximum(1.05);
	h1->GetYaxis()->SetNoExponent();
	h1->GetXaxis()->SetNoExponent();
	h1->GetXaxis()->SetRangeUser(0,40);

	TCanvas *c1 = tdrCanvas("c1",h1,0,33);
    if (mult[0]->Integral()>0.0f) {
        mult[0]->Smooth(3);
	    mult[0]->Scale(1/mult[0]->GetMaximum());
	    tdrDraw(mult[0],"HIST P",kFullSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (mult[2]->Integral()>0.0f) {
        mult[2]->Smooth(3);
	    mult[2]->Scale(1/mult[2]->GetMaximum());
	    tdrDraw(mult[2],"HIST P",kFullTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (mult[1]->Integral()>0.0f) {
        mult[1]->Smooth(3);
	    mult[1]->Scale(1/mult[1]->GetMaximum());
	    tdrDraw(mult[1],"HIST P",kFullCircle,kBlue,kSolid,-1,3005,kBlue);
    }

    if (mult[3]->Integral()>0.0f) {
        mult[3]->Smooth(3);
	    mult[3]->Scale(1/mult[3]->GetMaximum());
	    tdrDraw(mult[3],"HIST P",kOpenSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (mult[5]->Integral()>0.0f) {
        mult[5]->Smooth(3);
	    mult[5]->Scale(1/mult[5]->GetMaximum());
	    tdrDraw(mult[5],"HIST P",kOpenTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (mult[4]->Integral()>0.0f) {
        mult[4]->Smooth(3);
	    mult[4]->Scale(1/mult[4]->GetMaximum());
	    tdrDraw(mult[4],"HIST P",kOpenCircle,kBlue,kSolid,-1,3005,kBlue);
    }
    TLegend *multLeg = tdrLeg(0.15,0.6,0.35,0.9);
    multLeg->AddEntry(mult[0],"gluon ttbar","p");      
    multLeg->AddEntry(mult[3],"gluon dijet","p");      
    multLeg->AddEntry(mult[2],"b quark ttbar","p");    
    multLeg->AddEntry(mult[5],"b quark ttbar","p");    
    multLeg->AddEntry(mult[1],"non-b quark ttbar","p");
    multLeg->AddEntry(mult[4],"non-b quark dijet","p");
	
		
	/****************pTD****************/

	setTDRStyle();
	TH1D *h2 = new TH1D("h2",";p_{T}D;Events",100,0.15,0.85);
    h2->SetMinimum(0.0);
	h2->SetMaximum(1.05);
	h2->GetYaxis()->SetNoExponent();
	h2->GetXaxis()->SetNoExponent();

	TCanvas *c2 = tdrCanvas("c2",h2,0,33);
    if (ptd[0]->Integral()>0.0f) {
        ptd[0]->Smooth(3);
	    ptd[0]->Scale(1/ptd[0]->GetMaximum());
	    tdrDraw(ptd[0],"HIST P",kFullSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (ptd[2]->Integral()>0.0f) {
        ptd[2]->Smooth(3);
	    ptd[2]->Scale(1/ptd[2]->GetMaximum());
	    tdrDraw(ptd[2],"HIST P",kFullTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (ptd[1]->Integral()>0.0f) {
        ptd[1]->Smooth(3);
	    ptd[1]->Scale(1/ptd[1]->GetMaximum());
	    tdrDraw(ptd[1],"HIST P",kFullCircle,kBlue,kSolid,-1,3005,kBlue);
    }

    if (ptd[3]->Integral()>0.0f) {
        ptd[3]->Smooth(3);
	    ptd[3]->Scale(1/ptd[3]->GetMaximum());
	    tdrDraw(ptd[3],"HIST P",kOpenSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (ptd[5]->Integral()>0.0f) {
        ptd[5]->Smooth(3);
	    ptd[5]->Scale(1/ptd[5]->GetMaximum());
	    tdrDraw(ptd[5],"HIST P",kOpenTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (ptd[4]->Integral()>0.0f) {
        ptd[4]->Smooth(3);
	    ptd[4]->Scale(1/ptd[4]->GetMaximum());
	    tdrDraw(ptd[4],"HIST P",kOpenCircle,kBlue,kSolid,-1,3005,kBlue);
    }

    TLegend *ptdLeg = tdrLeg(0.15,0.6,0.35,0.9);
    ptdLeg->AddEntry(mult[0],"gluon ttbar","p");      
    ptdLeg->AddEntry(mult[3],"gluon dijet","p");      
    ptdLeg->AddEntry(mult[2],"b quark ttbar","p");    
    ptdLeg->AddEntry(mult[5],"b quark ttbar","p");    
    ptdLeg->AddEntry(mult[1],"non-b quark ttbar","p");
    ptdLeg->AddEntry(mult[4],"non-b quark dijet","p");

	/****************sigma2****************/
	
    setTDRStyle();
	TH1D *h3 = new TH1D("h3",";#sigma_{2};Events",100,0,0.14);
	h3->SetMinimum(0);
	h3->SetMaximum(1.05);
	h3->GetYaxis()->SetNoExponent();
	h3->GetXaxis()->SetNoExponent();
	h3->GetXaxis()->SetRangeUser(0,0.2);

	TCanvas *c3 = tdrCanvas("c3",h3,0,33);
    if (s2[0]->Integral()>0.0f) {
	    s2[0]->Smooth(3);
        s2[0]->Scale(1/s2[0]->GetMaximum());
	    tdrDraw(s2[0],"HIST P",kFullSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (s2[2]->Integral()>0.0f) {
        s2[2]->Smooth(3);
	    s2[2]->Scale(1/s2[2]->GetMaximum());
	    tdrDraw(s2[2],"HIST P",kFullTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (s2[1]->Integral()>0.0f) {
        s2[1]->Smooth(3);
	    s2[1]->Scale(1/s2[1]->GetMaximum());
	    tdrDraw(s2[1],"HIST P",kFullCircle,kBlue,kSolid,-1,3005,kBlue);
    }

    if (s2[3]->Integral()>0.0f) {
	    s2[3]->Smooth(3);
        s2[3]->Scale(1/s2[3]->GetMaximum());
	    tdrDraw(s2[3],"HIST P",kOpenSquare,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (s2[5]->Integral()>0.0f) {
        s2[5]->Smooth(3);
	    s2[5]->Scale(1/s2[5]->GetMaximum());
	    tdrDraw(s2[5],"HIST P",kOpenTriangleUp,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (s2[4]->Integral()>0.0f) {
        s2[4]->Smooth(3);
	    s2[4]->Scale(1/s2[4]->GetMaximum());
	    tdrDraw(s2[4],"HIST P",kOpenCircle,kBlue,kSolid,-1,3005,kBlue);
    }

    TLegend *s2Leg = tdrLeg(0.65,0.6,0.85,0.9);
    s2Leg->AddEntry(mult[0],"gluon ttbar","p");
    s2Leg->AddEntry(mult[3],"gluon dijet","p");
    s2Leg->AddEntry(mult[1],"b quark ttbar","p");
    s2Leg->AddEntry(mult[4],"b quark ttbar","p");
    s2Leg->AddEntry(mult[2],"non-b quark ttbar","p");
    s2Leg->AddEntry(mult[5],"non-b quark dijet","p");
}

