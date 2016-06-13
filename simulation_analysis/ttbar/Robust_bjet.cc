#include "common_ttbar.h"

bool unscaled = false;
double PTDMax = 0.05;
double S2Max = 0.07;
double MultMax = 0.08;

void Plot(string fileName)
{
/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");

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

    unsigned acceptance = 0;
    for(size_t x=0; x != N; ++x)
    {
        jetTree->GetEntry(x);
        assert(kMaxfJets>fJets);
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
			    multiplicity_g->Fill(fConstituents[i],fWeight);
			    pTD_g->Fill(fPTD[i],fWeight);
			    sigma2_g->Fill(fSigma2[i],fWeight);
		    } else if(id == 5) {
			    multiplicity_u->Fill(fConstituents[i],fWeight);
			    pTD_u->Fill(fPTD[i],fWeight);
			    sigma2_u->Fill(fSigma2[i],fWeight);
		    } else if (id > 0) {
			    multiplicity_q->Fill(fConstituents[i],fWeight);
			    pTD_q->Fill(fPTD[i],fWeight);
			    sigma2_q->Fill(fSigma2[i],fWeight);	
		    }
        }
    }
    cout << "Accepted: " << acceptance << endl;
	
	/*****************Constituents*****************/

	//STACKED
	setTDRStyle();
	TH1F *um = (TH1F*)multiplicity_u->Clone("um");
	TH1F *gm = (TH1F*)multiplicity_g->Clone("gm");
	TH1F *qm = (TH1F*)multiplicity_q->Clone("qm");

	TH1D *h4 = new TH1D("h4",";Number of constituents;Events",60,0,60);
	h4->SetMinimum(0);
	h4->SetMaximum(MultMax);
	h4->GetYaxis()->SetNoExponent();
	h4->GetXaxis()->SetNoExponent();
	h4->GetXaxis()->SetRangeUser(0,60);
	multiplicity_g->Add(multiplicity_u);
	multiplicity_q->Add(multiplicity_g);

    if (unscaled) {
	    TCanvas *c5 = tdrCanvas("c5",h4,0,33);
        if (multiplicity_g->Integral()>0.0f)
	        tdrDraw(multiplicity_g,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
        if (multiplicity_u->Integral()>0.0f)
	        tdrDraw(multiplicity_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        if (multiplicity_q->Integral()>0.0f)
	        tdrDraw(multiplicity_q,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
    }

	//SEPARATE
	setTDRStyle();
	TH1D *h1 = new TH1D("h1",";Number of constituents;Events",60,0,60);
	h1->SetMinimum(0);
	h1->SetMaximum(1.05);
	h1->GetYaxis()->SetNoExponent();
	h1->GetXaxis()->SetNoExponent();
	h1->GetXaxis()->SetRangeUser(0,60);

	TCanvas *c2 = tdrCanvas("c2",h1,0,33);
    if (gm->Integral()>0.0f) {
        gm->Smooth();
	    gm->Scale(1/gm->GetMaximum());
	    tdrDraw(gm,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (um->Integral()>0.0f) {
        um->Smooth();
	    um->Scale(1/um->GetMaximum());
	    tdrDraw(um,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (qm->Integral()>0.0f) {
        qm->Smooth();
	    qm->Scale(1/qm->GetMaximum());
	    tdrDraw(qm,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
    }
	
		
	/****************pTD****************/

	//STACKED
	setTDRStyle();
	TH1F *uD = (TH1F*)pTD_u->Clone("uD");
	TH1F *gD = (TH1F*)pTD_g->Clone("gD");
	TH1F *qD = (TH1F*)pTD_q->Clone("qD");

	TH1D *h2 = new TH1D("h2",";p_{T}D;Events",100,0,1);
	h2->SetMaximum(PTDMax);
	h2->GetYaxis()->SetNoExponent();
	h2->GetXaxis()->SetNoExponent();
	h2->GetXaxis()->SetRangeUser(0,30);
	pTD_g->Add(pTD_u);
	pTD_q->Add(pTD_g);

    if (unscaled) {
	    TCanvas *c3 = tdrCanvas("c3",h2,0,33);
        if (pTD_g->Integral()>0.0f)
	        tdrDraw(pTD_g,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
        if (pTD_u->Integral()>0.0f)
	        tdrDraw(pTD_u,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
        if (pTD_q->Integral()>0.0f)
	        tdrDraw(pTD_q,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
	}

	//SEPARATE
	setTDRStyle();
	TH1D *h5 = new TH1D("h5",";p_{T}D;Events",100,0,1);
	h5->SetMaximum(1.05);
	h5->GetYaxis()->SetNoExponent();
	h5->GetXaxis()->SetNoExponent();
	h5->GetXaxis()->SetRangeUser(0,30);

	TCanvas *c6 = tdrCanvas("c6",h5,0,33);
    if (gD->Integral()>0.0f) {
        gD->Smooth();
	    gD->Scale(1/gD->GetMaximum());
	    tdrDraw(gD,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (uD->Integral()>0.0f) {
        uD->Smooth();
	    uD->Scale(1/uD->GetMaximum());
	    tdrDraw(uD,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (qD->Integral()>0.0f) {
        qD->Smooth();
	    qD->Scale(1/qD->GetMaximum());
	    tdrDraw(qD,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
    }

	/****************sigma2****************/
	
	//STACKED
	setTDRStyle();
	TH1F *us = (TH1F*)sigma2_u->Clone("us");
	TH1F *gs = (TH1F*)sigma2_g->Clone("gs");
	TH1F *qs = (TH1F*)sigma2_q->Clone("qs");

	TH1D *h6 = new TH1D("h6",";#sigma_{2};Events",100,0,0.2);
	h6->SetMinimum(0);
	h6->SetMaximum(S2Max);
	h6->GetYaxis()->SetNoExponent();
	h6->GetXaxis()->SetNoExponent();
	h6->GetXaxis()->SetRangeUser(0,0.2);
	sigma2_g->Add(sigma2_u);
	sigma2_q->Add(sigma2_g);

    if (unscaled) {
	    TCanvas *c7 = tdrCanvas("c7",h6,0,33);
        if (sigma2_g->Integral()>0.0f)
	        tdrDraw(sigma2_g,"HIST",kDot,kRed-3,kSolid,-1,3004,kRed-3);
        if (sigma2_u->Integral()>0.0f)
	        tdrDraw(sigma2_u,"HIST",kDot,kGreen-1,kSolid,-1,3005,kGreen-1);
        if (sigma2_q->Integral()>0.0f)
	        tdrDraw(sigma2_q,"HIST",kDot,kBlue,kSolid,-1,3003,kBlue);
	}

	//SEPARATE
	setTDRStyle();
	TH1D *h3 = new TH1D("h3",";#sigma_{2};Events",100,0,0.2);
	h3->SetMinimum(0);
	h3->SetMaximum(1.05);
	h3->GetYaxis()->SetNoExponent();
	h3->GetXaxis()->SetNoExponent();
	h3->GetXaxis()->SetRangeUser(0,0.2);

	TCanvas *c4 = tdrCanvas("c4",h3,0,33);
    if (gs->Integral()>0.0f) {
	    gs->Smooth();
        gs->Scale(1/gs->GetMaximum());
	    tdrDraw(gs,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
    }
    if (us->Integral()>0.0f) {
        us->Smooth();
	    us->Scale(1/us->GetMaximum());
	    tdrDraw(us,"HIST",kDot,kGreen-1,kSolid,-1,3004,kGreen-1);
    }
    if (qs->Integral()>0.0f) {
        qs->Smooth();
	    qs->Scale(1/qs->GetMaximum());
	    tdrDraw(qs,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
    }
}

