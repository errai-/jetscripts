#include "common_ttbar.h"

bool merger = false;
int falsetag = 0;
int etaBins = 40;
double etaMax = 2.5;

void Stack(string file) {

/* Initialization */
    gROOT->ProcessLine(".L sim_dir/lib/libJetEvent.so");
    TProfile gluonFrac("g","g",ptBins,ptRange);
    TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
    TProfile strangeFrac("s","s",ptBins,ptRange);
    TProfile charmFrac("c","c",ptBins,ptRange);
    TProfile bottomFrac("b","b",ptBins,ptRange);
    TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);
    TProfile gEtaFrac("ge","ge",etaBins,-etaMax,etaMax);
    TProfile lqEtaFrac("lqe","lqe",etaBins,-etaMax,etaMax);
    TProfile sEtaFrac("se","se",etaBins,-etaMax,etaMax);
    TProfile cEtaFrac("ce","ce",etaBins,-etaMax,etaMax);
    TProfile bEtaFrac("be","be",etaBins,-etaMax,etaMax);
    TProfile uEtaFrac("ue","ue",etaBins,-etaMax,etaMax);

    static const Int_t kMaxfJets = 100;

    Int_t           fJets;
    Double_t        fX[kMaxfJets];   //[mJets]
    Double_t        fY[kMaxfJets];   //[mJets]
    Double_t        fZ[kMaxfJets];   //[mJets]
    Double_t        fT[kMaxfJets];   //[mJets]

    Double_t        fWeight;
    Int_t           fFlav[kMaxfJets];   //[mJets]

    /* Tree setup */
    TTree* jetTree;

    TChain* jetChain = new TChain("JetTree","");
    jetChain->Add(file.c_str()); jetTree = jetChain;

    jetTree->SetMakeClass(1);

    jetTree->SetBranchAddress("fJets", &fJets);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fX", fX);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fY", fY);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fZ", fZ);
    jetTree->SetBranchAddress("fJets.fP4.fCoordinates.fT", fT);
    jetTree->SetBranchAddress("fWeight", &fWeight);
    jetTree->SetBranchAddress("fJets.fFlav", fFlav);

    /* event loop */
    unsigned counter = 0;
    unsigned N = jetTree->GetEntries();
    int success = 0, nonb1 = 0, nonb0 = 0, noflav = 0, nolw = 0, noqw = 0, nots = 0;
    for(unsigned x=0; x != N; ++x) {
        jetTree->GetEntry(x);

        assert(kMaxfJets>fJets);

        vector<unsigned> flavours;
        vector<TLorentzVector> bjets, ljets;
        TLorentzVector MET, lepton, neutrino;

        bool stop = false;
        int bcount = 0;
        int flavors = 0;
        int jets = 0;
        for (int i = 0; i < fJets; ++i) {
            TLorentzVector tmpVec(fX[i],fY[i],fZ[i],fT[i]);
            int fl = abs(fFlav[i]);
            if ( fl >= 10 && fl < 20 )
                continue;
            if ( i<4 && (fabs(tmpVec.Eta()) > 2.5 || tmpVec.Pt()<30) ) {
                stop = true;
                break;
            }
            ++jets;
            if (fl!=0)
                ++flavors;
            if (fl == 5)
                ++bcount;
        }

        if (stop)
            continue;
        //if (bcount!=2)
        //    continue;
        //if (flavors > 4)
        //    ++falsetag;
        //if (flavors!=4)
        //    continue;
        ++counter;

        vector<double> pt_store;
        for (int i = 0; i < fJets; ++i) {
            TLorentzVector tmpVec(fX[i],fY[i],fZ[i],fT[i]);
            int fl = abs(fFlav[i]);

            if (fl == 6 || fl == 7) fl -= 2;
            if (fl > 5 && fl != 21) continue;
            //if (fl == 0) {
            //    if (i < 4)
            //        pt_store.push_back( tmpVec.Pt() );
            //    continue;
            //}

            gluonFrac.Fill(tmpVec.Pt(), (fl == 21)? 1:0, fWeight);
            lightquarkFrac.Fill(tmpVec.Pt(), (fl == 1 || fl == 2)? 1:0, fWeight);
            strangeFrac.Fill(tmpVec.Pt(), (fl == 3)? 1:0, fWeight);
            charmFrac.Fill(tmpVec.Pt(), (fl == 4)? 1:0, fWeight);
            bottomFrac.Fill(tmpVec.Pt(), (fl == 5)? 1:0, fWeight);
            unmatchedFrac.Fill(tmpVec.Pt(), (fl <= 0)? 1:0, fWeight);

            if (tmpVec.Pt()<30)
                continue;

            gEtaFrac.Fill(tmpVec.Eta(), (fl == 21)? 1:0, fWeight);
            lqEtaFrac.Fill(tmpVec.Eta(), (fl == 1 || fl == 2)? 1:0, fWeight);
            sEtaFrac.Fill(tmpVec.Eta(), (fl == 3)? 1:0, fWeight);
            cEtaFrac.Fill(tmpVec.Eta(), (fl == 4)? 1:0, fWeight);
            bEtaFrac.Fill(tmpVec.Eta(), (fl == 5)? 1:0, fWeight);
            uEtaFrac.Fill(tmpVec.Eta(), (fl <= 0)? 1:0, fWeight);
        }
        //unsigned pos = 0;
        //for (unsigned i = 0; pos<pt_store.size() && i < 2-min(bcount,2); ++i, ++pos) {
        //    gluonFrac.Fill(pt_store[pos], 0, fWeight);
        //    lightquarkFrac.Fill(pt_store[pos], 0, fWeight);
        //    strangeFrac.Fill(pt_store[pos], 0, fWeight);
        //    charmFrac.Fill(pt_store[pos], 0, fWeight);
        //    bottomFrac.Fill(pt_store[pos], 0, fWeight);
        //    unmatchedFrac.Fill(pt_store[pos], 1, fWeight);
        //}
        //for (unsigned i = 0; pos<pt_store.size() && i < 2-(flavors-min(bcount,2)); ++i, ++pos) {
        //    gluonFrac.Fill(pt_store[pos], 1, fWeight);
        //    lightquarkFrac.Fill(pt_store[pos], 0, fWeight);
        //    strangeFrac.Fill(pt_store[pos], 0, fWeight);
        //    charmFrac.Fill(pt_store[pos], 0, fWeight);
        //    bottomFrac.Fill(pt_store[pos], 0, fWeight);
        //    unmatchedFrac.Fill(pt_store[pos], 0, fWeight);
        //}

        //for (int i = 0; i < mJets; ++i) {
        //    TLorentzVector tmpVec(mX[i],mY[i],mZ[i],mT[i]);
        //    if (mFlav[i]==10)
        //        MET = tmpVec;
        //    else if (mFlav[i]==11||mFlav[i]==13||mFlav[i]==15)
        //        lepton = tmpVec;
        //    else if (mFlav[i]==12||mFlav[i]==14||mFlav[i]==16)
        //        neutrino = tmpVec;
        //    else {
        //        if (tmpVec.Pt() < 30)
        //            continue;

        //        if (mFlav[i]==5)
        //            bjets.push_back(tmpVec);
        //        else {
        //            if (mFlav[i]!=0)
        //                ++flav_count;
        //            ljets.push_back(tmpVec);
        //            flavours.push_back(mFlav[i]);
        //        }
        //    }
        //}
        //if (bjets.size()!=2) {
        //    if (bjets.size()==1)
        //        ++nonb1;
        //    else
        //        ++nonb0;
        //    continue;
        //}
        //if (flav_count!=2) {
        //    ++noflav;
        //    continue;
        //}
        //// Ideal MET
        //MET.SetPz( neutrino.Pz() );
        //MET.SetE( MET.P() );

        //// Containers for combined lorentz vectors
        //TLorentzVector t1, t1_alt, t2, t3, t4, t5, t6;

        //// Reconstruct W mass from MET and the lepton
        //t1_alt = MET + lepton;
        //t1 = neutrino + lepton;
        //if (t1.M() < 60 || t1.M() > 110) {
        //    ++nolw;
        //    continue;
        //}

        //// Find jet pairs that correspond to the W
        //vector<TLorentzVector> working;
        //vector< pair<unsigned,unsigned> > working_idx;
        //for (auto i = 0u; i < ljets.size()-1; ++i) {
        //    if (flavours[i]==0) continue;
        //    for (auto j = i+1; j < ljets.size(); ++j) {
        //        if (flavours[j]==0) continue;
        //        t2 = ljets[i] + ljets[j];
        //        if (t2.M() > 60 && t2.M() < 110) {
        //            working.push_back(t2);
        //            working_idx.push_back( std::make_pair(i,j) );
        //        }
        //    }
        //}
        //if (working.size() == 0) {
        //    ++noqw;
        //    continue;
        //}

        //// Pairings of the lepton-W with bjets
        //t3 = t1 + bjets[0];
        //t4 = t1 + bjets[1];
        //unsigned tmatch = 0;
        //unsigned best;
        //int id;

        //// Find pairings of the quark-W with bjets
        //for ( unsigned i = 0; i < working.size(); ++i ) {
        //    t5 = working[i] + bjets[0];
        //    t6 = working[i] + bjets[1];
        //    if (mass_study(t3.M(),t4.M(),t5.M(),t6.M(),false,id)) {
        //        ++tmatch;
        //        t2 = working[i];
        //        best = i;
        //    }
        //}
        //if (tmatch == 0) {
        //    ++nots;
        //    continue;
        //}
    }

    cout << counter << " events analyzed" << endl;
    cout << falsetag << " false flavor tags" << endl;

    TH1D *light_quarks = lightquarkFrac.ProjectionX("light quarks","");
    TH1D *gluons = gluonFrac.ProjectionX("gluons","");
    TH1D *strange = strangeFrac.ProjectionX("strange","");
    if(merger) light_quarks->Add(strange);
    TH1D *charm = charmFrac.ProjectionX("charm","");
    TH1D *bottom = bottomFrac.ProjectionX("bottom","");
    TH1D *unmatched = unmatchedFrac.ProjectionX("unmatched","");
    
    TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);

    tdrDraw(unmatched,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(strange,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

    THStack *hs  = new THStack("hs","test stacked histograms");

    TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

    //light_quarks->Add(strange);
    hs->Add(bottom);
    hs->Add(charm);
    if(!merger)hs->Add(strange);
    hs->Add(light_quarks);
    hs->Add(gluons);
    hs->Add(unmatched);

    hs->Draw();

    hs->GetXaxis()->SetNoExponent();
    // hs->GetXaxis()->SetNoExponent(kTRUE);
    hs->GetXaxis()->SetRange(9,47);
    hs->GetXaxis()->SetMoreLogLabels(kTRUE);
    hs->GetXaxis()->SetTitle("p_{T} (GeV)");
    hs->GetYaxis()->SetTitle("Flavor fraction");
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetTitleOffset(1.2);
    hs->GetXaxis()->SetTitleSize(0.05);
    hs->GetXaxis()->SetTitleOffset(1);
    // hs->SetLogx();
    hs->SetMaximum(0.97);
    hs->GetXaxis()->UnZoom();

    double x0, y0;
    x0 = 0.4;
    y0 = 0.05;

    // //hadronic def
    // TLegend *leg = tdrLeg(0.642617,0.304878,0.968121,0.623693);              
    // TLegend *sample = tdrLeg(0.173658,0.3223-0.01,0.729027,0.54878-0.01);            //
    // TLegend *alphacut = tdrLeg(0.162752+0.16,0.336237,0.708054+0.16,0.409408);           //hadronic
    // TLegend *etacut = tdrLeg(0.166107,0.334495 ,0.709732,0.407666);          

    //QCDaware def
    //TLegend *leg = tdrLeg(0.5+0.5,0.82-0.2,0.175+0.5,0.50-0.2);           
    //TLegend *sample = tdrLeg(0.675-0.05-x0,0.50-y0,0.775-0.05-x0,0.505-y0);               //QCDaware
    //TLegend *alphacut = tdrLeg(0.77-x0,0.50-0.05-y0,0.87-x0,0.505-0.05-y0);               //goes
    //TLegend *etacut = tdrLeg(0.61-x0,0.50-0.05-y0,0.71-x0,0.505-0.05-y0);             //here

    ////physics def 
    TLegend *leg = tdrLeg(0.61,0.45,0.31,0.16);
    TLegend *sample = tdrLeg(0.175,0.86,0.275,0.865);         //
    TLegend *alphacut = tdrLeg(0.77,0.50-0.05,0.87,0.505-0.05);         //physics
    TLegend *etacut = tdrLeg(0.175,0.86-0.05,0.275,0.865-0.05);           //

    //sample->SetHeader("#gamma+jet sample");
    sample->SetHeader("ttbarlepton+jet sample");
    //TLegend *heading = tdrLeg(0.675-0.4,0.50+0.5,0.775-0.4,0.505+0.5);·   
    //heading->SetHeader("Hadronic Definition, #sqrt{s} = 8 TeV");
    alphacut->SetHeader("");
    etacut->SetHeader("#left|#eta#right|< 2.5");

    leg->AddEntry(unmatched,"None","f");    
    leg->AddEntry(gluons,"Gluon","f");
    leg->AddEntry(light_quarks,"Light","f");
    leg->AddEntry(strange,"Strange","f");
    leg->AddEntry(charm,"Charm","f");
    leg->AddEntry(bottom,"Bottom","f");

    gPad->SetLogx();

// Stop

    TH1D *lqe = lqEtaFrac.ProjectionX("light quarks e","");
    TH1D *ge = gEtaFrac.ProjectionX("gluons e","");
    TH1D *se = sEtaFrac.ProjectionX("strange e","");
    if(merger) lqe->Add(se);
    TH1D *ce = cEtaFrac.ProjectionX("charm e","");
    TH1D *be = bEtaFrac.ProjectionX("bottom e","");
    TH1D *ue = uEtaFrac.ProjectionX("unmatched e","");
    
    TH1D *he = new TH1D("he",";#eta;Fraction",etaBins,-etaMax,etaMax);

    tdrDraw(ue,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
    gStyle->SetOptStat(kFALSE); //removes old legend
    tdrDraw(ge,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
    tdrDraw(lqe,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
    tdrDraw(se,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
    tdrDraw(ce,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
    tdrDraw(be,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

    THStack *hse  = new THStack("hs","test stacked histograms");

    TCanvas *c2 = tdrCanvas("c2",he,2,0,kSquare);

    //light_quarks->Add(strange);
    hse->Add(be);
    hse->Add(ce);
    if(!merger)hse->Add(se);
    hse->Add(lqe);
    hse->Add(ge);
    hse->Add(ue);

    hse->Draw();

    hse->GetXaxis()->SetNoExponent();
    // hs->GetXaxis()->SetNoExponent(kTRUE);
    hse->GetXaxis()->SetRange(-etaMax,etaMax);
    hse->GetXaxis()->SetMoreLogLabels(kTRUE);
    hse->GetXaxis()->SetTitle("#eta");
    hse->GetYaxis()->SetTitle("Flavor fraction");
    hse->GetYaxis()->SetTitleSize(0.05);
    hse->GetYaxis()->SetTitleOffset(1.2);
    hse->GetXaxis()->SetTitleSize(0.05);
    hse->GetXaxis()->SetTitleOffset(1);
    // hs->SetLogx();
    hse->SetMaximum(0.97);
    hse->GetXaxis()->UnZoom();

    // //hadronic def
    // TLegend *leg = tdrLeg(0.642617,0.304878,0.968121,0.623693);              
    // TLegend *sample = tdrLeg(0.173658,0.3223-0.01,0.729027,0.54878-0.01);            //
    // TLegend *alphacut = tdrLeg(0.162752+0.16,0.336237,0.708054+0.16,0.409408);           //hadronic
    // TLegend *etacut = tdrLeg(0.166107,0.334495 ,0.709732,0.407666);          

    //QCDaware def
    //TLegend *leg = tdrLeg(0.5+0.5,0.82-0.2,0.175+0.5,0.50-0.2);           
    //TLegend *sample = tdrLeg(0.675-0.05-x0,0.50-y0,0.775-0.05-x0,0.505-y0);               //QCDaware
    //TLegend *alphacut = tdrLeg(0.77-x0,0.50-0.05-y0,0.87-x0,0.505-0.05-y0);               //goes
    //TLegend *etacut = tdrLeg(0.61-x0,0.50-0.05-y0,0.71-x0,0.505-0.05-y0);             //here

    ////physics def 
    TLegend *lege = tdrLeg(0.61,0.45,0.31,0.16);
    TLegend *samplee = tdrLeg(0.175,0.86,0.275,0.865);         //
    TLegend *alphacute = tdrLeg(0.77,0.50-0.05,0.87,0.505-0.05);         //physics
    TLegend *etacute = tdrLeg(0.175,0.86-0.05,0.275,0.865-0.05);           //

    //sample->SetHeader("#gamma+jet sample");
    samplee->SetHeader("ttbarlepton+jet sample");
    //TLegend *heading = tdrLeg(0.675-0.4,0.50+0.5,0.775-0.4,0.505+0.5);·   
    //heading->SetHeader("Hadronic Definition, #sqrt{s} = 8 TeV");
    alphacute->SetHeader("");
    etacute->SetHeader("#left|#eta#right|< 2.5");

    lege->AddEntry(unmatched,"None","f");    
    lege->AddEntry(gluons,"Gluon","f");
    lege->AddEntry(light_quarks,"Light","f");
    lege->AddEntry(strange,"Strange","f");
    lege->AddEntry(charm,"Charm","f");
    lege->AddEntry(bottom,"Bottom","f");
}
