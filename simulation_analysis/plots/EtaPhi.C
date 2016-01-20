#define EtaPhi_cxx
#include "etaphi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;

bool accept() {
    string ses;
    cin >> ses;
    if ( ses == "y" ) return true;
    else return false;
}

void EtaPhi::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);

//        if ( jentry == 0 ) {

        //if ( jentry == 53 || jentry == 71 || jentry == 163  jentry == 205 )
        //    continue;

            bool stahp = false;
            int count = 0;
            for ( int j = 0; j < fJets; ++j ) {
                int fl = fFlav[j];
                TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
                if (fl < 0 || fl == 10) {
                    continue;
                }

                //if ( fDR[j] > 0.7 ) {
                //    ++count;
                //}

                if ( fl != 0 )
                    ++count;

                if ( fabs(t.Eta()) > 3 ) {
                    stahp = true;
                    break;
                }
            }
            if ( stahp || count > 2 ) {
                continue;
            }
            cout << jentry << " " << count << endl;
            if ( !accept() )
                continue;

            TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
            c1->SetLogz(1);
            TH2D* grid1 = new TH2D("g1","g1",50,-5,5,50,-TMath::Pi(),TMath::Pi());
            grid1->SetMaximum(10000);
            grid1->SetMinimum(0);

            TCanvas* c2 = new TCanvas("c2", "c2", 1200, 1000);
            c2->SetLogz(1);
            TH2D* grid2 = new TH2D("g2","g2",50,-5,5,50,-TMath::Pi(),TMath::Pi());
            grid2->SetMaximum(10000);
            grid2->SetMinimum(0);

            TCanvas* c3 = new TCanvas("c3", "c3", 1200, 1000);
            c3->SetLogz(1);
            TH2D* grid3 = new TH2D("g3","g3",50,-5,5,50,-TMath::Pi(),TMath::Pi());
            grid3->SetMaximum(10000);
            grid3->SetMinimum(0);

            TCanvas* c4 = new TCanvas("c4", "c4", 1200, 1000);
            c4->SetLogz(1);
            TH2D* grid4 = new TH2D("g4","g4",50,-5,5,50,-TMath::Pi(),TMath::Pi());
            grid4->SetMaximum(10000);
            grid4->SetMinimum(0);

            TCanvas* c5 = new TCanvas("c5", "c5", 1200, 1000);
            c5->SetLogz(1);
            TH2D* grid5 = new TH2D("g5","g5",50,-5,5,50,-TMath::Pi(),TMath::Pi());
            grid5->SetMaximum(10000);
            grid5->SetMinimum(0);

            vector<int> counters(5,0);
            for ( int j = 0; j < fJets; ++j ) {
                int fl = fFlav[j];
                TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);

                if (fl == 10)
                    continue;

                fl *= -1;


                if ( fl == 7 || fl == 8 ) { // Final state
                    ++counters[1];
                    grid2->Fill(t.Eta(),t.Phi(),t.E());
                } else if ( fl == 9 ) { // Hard proc + correction
                    ++counters[3];
                    grid4->Fill(t.Eta(),t.Phi(),t.E());
                } else if ( fl == 10 || fl == 11 ) { //final parton state
                    ++counters[4];
                    grid5->Fill(t.Eta(),t.Phi(),t.E());
                } else if ( fl == 12 || fl == 13 ) { // hard proc again
                    ++counters[2];
                    grid3->Fill(t.Eta(),t.Phi(),t.E());
                }else if ( (fl > 0 && fl <= 6) || fl==21 ) { // Jet const.
                    ++counters[0];
                    grid1->Fill(t.Eta(),t.Phi(),t.E());
                } else if ( fl <= 0 ) { // jets
                    ++counters[0];
                    grid1->Fill(t.Eta(),t.Phi(),100*t.E());
                    cout << fl << " " << fDR[j] << " " << fAlpha[j] << " " << t.Pt() << endl;
                }
            }
            gStyle->SetOptStat(0);
            grid5->Draw("COLZ");
            c4->cd();
            gStyle->SetOptStat(0);
            grid4->Draw("COLZ");
            c3->cd();
            gStyle->SetOptStat(0);
            grid3->Draw("COLZ");
            c2->cd();
            gStyle->SetOptStat(0);
            grid2->Draw("COLZ");
            c1->cd();
            gStyle->SetOptStat(0);
            grid1->Draw("COLZ");
            for ( int i = 0; i < 5; ++i )
                cout << " " << counters[i] << endl;
            break;
//        }
   }

}
