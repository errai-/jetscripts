#include "common_ttbar.h"

void TTBar(string file) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    TH1D *Wlepton = new TH1D("","m_W lepton",noW,lowerW,upperW);
    TH1D *Wtlepton = new TH1D("","m_W alt lepton",noW,lowerW,upperW);
    TH1D *Wjet = new TH1D("","m_W jet",noW,lowerW,upperW);
    TH1D *Wjetc = new TH1D("","m_W jet",200,0.0,300.0);
    TH1D *tlepton = new TH1D("","m_t lepton",noT,lowerT,upperT);
    TH1D *ttlepton = new TH1D("","m_t alt lepton",noT,lowerT,upperT);
    TH1D *tjet = new TH1D("","m_t jet",noT,lowerT,upperT);
    TH1D *bboth = new TH1D("","m_b jet",100,3,23);
    TH1D *rbq = new TH1D("","pTb/pTW",100,0,3);
    TH1D *Wfrac = new TH1D("","frac",100,0.5,2);

    static const Int_t kMaxfJets = 100;

    Int_t           mJets;
    Double_t        mX[kMaxfJets];   //[mJets]
    Double_t        mY[kMaxfJets];   //[mJets]
    Double_t        mZ[kMaxfJets];   //[mJets]
    Double_t        mT[kMaxfJets];   //[mJets]

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

    /* event loop */
    std::size_t mCount = 0;
    std::size_t mN = jetTree->GetEntries();
    int success = 0, nonb1 = 0, nonb0 = 0, noflav = 0, nolw = 0, noqw = 0, nots = 0;
    for(size_t x=0; x != mN; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>mJets);

        vector<unsigned> flavours;
        vector<TLorentzVector> bjets, ljets;
        TLorentzVector MET, lepton, neutrino;
        unsigned flav_count = 0;
        for (int i = 0; i < mJets; ++i) {
            TLorentzVector tmpVec(mX[i],mY[i],mZ[i],mT[i]);
            if (mFlav[i]==10)
                MET = tmpVec;
            else if (mFlav[i]==11||mFlav[i]==13||mFlav[i]==15)
                lepton = tmpVec;
            else if (mFlav[i]==12||mFlav[i]==14||mFlav[i]==16)
                neutrino = tmpVec;
            else {
                if (tmpVec.Pt() < 30 || tmpVec.Eta() > 2.5)
                    continue;

                if (mFlav[i]==5)
                    bjets.push_back(tmpVec);
                else {
                    ljets.push_back(tmpVec);
                    flavours.push_back(mFlav[i]);
                }
            }
        }
        if (bjets.size()!=2) {
            if (bjets.size()==1)
                ++nonb1;
            else
                ++nonb0;
            continue;
        }
        // Interesting tests
        //MET.SetPz( neutrino.Pz() );
        //MET.SetPz( 0 );
        TLorentzVector MET_plain = MET;
        MET_plain.SetPz( 0 );
        MET_plain.SetE( MET_plain.P() );
        MET.SetPz( pz_calc(lepton,MET,neutrino) );
        MET.SetE( MET.P() );
        // Containers for combined lorentz vectors
        TLorentzVector t1, t1_alt, t2, t3, t4, t5, t6;

        // Reconstruct W mass from MET and the lepton
        t1_alt = neutrino + lepton;
        t1 = MET + lepton;
        if (t1.M() < 60 || t1.M() > 110) {
            ++nolw;
            continue;
        }

        // Find jet pairs that correspond to the W
        vector<TLorentzVector> working;
        vector< pair<unsigned,unsigned> > working_idx;
        for (auto i = 0u; i < ljets.size()-1; ++i) {
            for (auto j = i+1; j < ljets.size(); ++j) {
                t2 = ljets[i] + ljets[j];
                if (t2.M() > 60 && t2.M() < 110) {
                    working.push_back(t2);
                    working_idx.push_back( std::make_pair(i,j) );
                }
            }
        }
        if (working.size() == 0) {
            ++noqw;
            continue;
        }

        // Pairings of the lepton-W with bjets
        t3 = t1 + bjets[0];
        t4 = t1 + bjets[1];
        unsigned tmatch = 0;
        unsigned best;
        double best_diff = 10000;
        int id;

        // Find pairings of the quark-W with bjets
        for ( unsigned i = 0; i < working.size(); ++i ) {
            t5 = working[i] + bjets[0];
            t6 = working[i] + bjets[1];
            double loc_diff;
            unsigned loc_mult = mass_study(t3.M(),t4.M(),t5.M(),t6.M(),false,id,loc_diff);
            if (loc_mult>0) {
                tmatch += loc_mult;
                if ( loc_diff < best_diff ) {
                    t2 = working[i];
                    best = i;
                    best_diff = loc_diff;
                }
            }
        }
        // HOX: this is an important constrain in its philosophical nature
        // takes away one fifth of the samples
        if (tmatch != 1) {
            ++nots;
            continue;
        }

        ++success;
        Wtlepton->Fill(t1_alt.M());
        //Wlepton->Fill(t1.M());
        Wlepton->Fill( (MET_plain+lepton).M() );
        Wjet->Fill(t2.M());
        double rbqval = (bjets[0].Pt() + bjets[1].Pt())/t2.Pt();
        Wjetc->Fill(rbqval*t2.M());
        if (id == 0) {
            ttlepton->Fill( (t1_alt+bjets[0]).M() );
            tlepton->Fill( (t1+bjets[0]).M() );
            tjet->Fill( (t2+bjets[1]).M() );
        } else {
            ttlepton->Fill( (t1_alt+bjets[1]).M() );
            tlepton->Fill( (t1+bjets[1]).M() );
            tjet->Fill( (t2+bjets[0]).M() );
        }
        rbq->Fill( rbqval );
        Wfrac->Fill( MET.Pz()/neutrino.Pz() );
        bboth->Fill( bjets[0].M() );
        bboth->Fill( bjets[1].M() );

//         cout << "Lepton W:" << t1.M() << " Jet W:" << t2.M() << endl;
//         cout << "Light flavours: " << flavours[working_idx[best].first] << " " << flavours[working_idx[best].second] << endl << endl;
    }
    TCanvas *c1 = new TCanvas("c1");
    Wlepton->SetYTitle("events");
    Wlepton->SetXTitle("m (GeV)");
    Wlepton->Draw();
    TCanvas *c11 = new TCanvas("c11");
    Wtlepton->SetYTitle("events");
    Wtlepton->SetXTitle("m (GeV)");
    Wtlepton->Draw();
    TCanvas *c2 = new TCanvas("c2");
    Wjet->SetYTitle("events");
    Wjet->SetXTitle("m (GeV)");
    Wjet->Draw();
    TCanvas *c3 = new TCanvas("c3");
    tlepton->SetYTitle("events");
    tlepton->SetXTitle("m (GeV)");
    tlepton->Draw();
    TCanvas *c31 = new TCanvas("c31");
    ttlepton->SetYTitle("events");
    ttlepton->SetXTitle("m (GeV)");
    ttlepton->Draw();
    TCanvas *c4 = new TCanvas("c4");
    tjet->SetYTitle("events");
    tjet->SetXTitle("m (GeV)");
    tjet->Draw();
    TCanvas *c8 = new TCanvas("c8");
    Wfrac->SetYTitle("events");
    Wfrac->SetXTitle("fraction");
    Wfrac->Draw();
    TCanvas *c9 = new TCanvas("c9");
    tlepton->Scale(1/tlepton->GetMaximum());
    ttlepton->Scale(1/ttlepton->GetMaximum());
    tjet->Scale(1/tjet->GetMaximum());
    tdrDraw(ttlepton,"C",kFullCircle,kBlue,kSolid,-1,0,kBlue);
    tdrDraw(tlepton,"C",kFullCircle,kRed-2,kSolid,-1,0,kRed-9);
    tdrDraw(tjet,"C",kFullCircle,kGreen-1,kSolid,-1,0,kGreen-1);
    gStyle->SetOptStat(0);
    TLegend *tplots = tdrLeg(0.65,0.6,0.85,0.9);
    tplots->AddEntry(ttlepton,"lepton, ideal","p");
    tplots->AddEntry(tlepton,"lepton, reconstructed","p");
    tplots->AddEntry(tjet,"jet, reconstructed","p");
    
    TCanvas *c10 = new TCanvas("c10");
    Wlepton->Scale(1/Wlepton->GetMaximum());
    Wtlepton->Scale(1/Wtlepton->GetMaximum());
    Wjet->Scale(1/Wjet->GetMaximum());
    tdrDraw(Wtlepton,"C",kFullCircle,kBlue,kSolid,-1,0,kBlue);
    tdrDraw(Wlepton,"C",kFullCircle,kRed-2,kSolid,-1,0,kRed-9);
    tdrDraw(Wjet,"C",kFullCircle,kGreen-1,kSolid,-1,0,kGreen-1);
    gStyle->SetOptStat(0);
    TLegend *Wplots = tdrLeg(0.65,0.6,0.85,0.9);
    Wplots->AddEntry(Wtlepton,"lepton, ideal","p");
    Wplots->AddEntry(Wlepton,"lepton, MET+lepton","p");
    Wplots->AddEntry(Wjet,"jet, reconstructed","p");
    
    cout << success << " successful matches." << endl;
    cout << nonb1 << " missing second b." << endl;
    cout << nonb0 << " missing first b." << endl;
    cout << noflav << " missing light flav." << endl;
    cout << nolw << " no lepton w match." << endl;
    cout << noqw << " no quark w match." << noqw << endl;
    cout << nots << " TTbar matching failed." << endl;
    cout << "Good: " << goods << " Bad: " << bads << endl;

    // Rather uninteresting canvases:

    //TCanvas *c5 = new TCanvas("c5");
    //bboth->SetYTitle("events");
    //bboth->SetXTitle("m (GeV)");
    //bboth->Draw();
    //TCanvas *c6 = new TCanvas("c6");
    //rbq->SetYTitle("events");
    //rbq->SetXTitle("rbq");
    //rbq->Draw();
    //TCanvas *c7 = new TCanvas("c7");
    //Wjetc->SetYTitle("events");
    //Wjetc->SetXTitle("m (GeV)");
    //Wjetc->Draw();
}
