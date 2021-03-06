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

double phif(double val, double center) {
    double corr = val - center;
    if ( corr < -TMath::Pi() )
        corr += 2*TMath::Pi();
    else if ( corr > TMath::Pi() )
        corr -= 2*TMath::Pi();
    return corr;
}

void EtaPhi::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    cout << "Choose event, -1 for no choice" << endl;
    int choice;
    cin >> choice;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);

        bool stahp = false;
        int count = 0;
        vector<int> jets;
        double phi_sum = 0;

        if ( choice >= 0 && choice != jentry )
            continue;

        cout << jentry << endl;
        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
            if (fl < 0 || fl == 10) {
                continue;
            }

            jets.push_back(j);

            cout << "  " << fl << endl;
            if ( fl != 0 )
                ++count;

            if ( fabs(t.Rapidity()) > 3 && jets.size() < 3 ) {
                stahp = true;
                break;
            }
            phi_sum += t.Phi();
        }

        if ( stahp || (count > 1 && choice == -1) ) {
            continue;
        }

        if ( !accept() )
            continue;
        TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
        TH1D* th = new TH1D("","",100,-5.0,5.0);
        gStyle->SetOptStat(0);
        th->SetMinimum(-TMath::Pi());
        th->SetMaximum(TMath::Pi());
        th->SetXTitle("y");
        th->SetYTitle("#phi");
        th->Draw("P");
        phi_sum /= jets.size();

        for ( unsigned j = 0; j < jets.size(); ++j ) {
            TLorentzVector t(fX[jets[j]],fY[jets[j]],fZ[jets[j]],fT[jets[j]]);
            TEllipse *el = new TEllipse(t.Rapidity(),phif(t.Phi(),phi_sum),0.5,0.5);
            Color_t ell_col = kRed+1;
            if (j == 1)
                ell_col = kRed;
            else if (j == 2)
                ell_col = kMagenta;
            else if (j == 3)
                ell_col = kCyan+2;
            else if (j > 3)
                ell_col = kBlue;
            el->SetLineColor(ell_col);
            el->SetFillStyle(0);
            el->Draw("same");
        }

        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);

            if (fl == 10)
                continue;

            fl *= -1;
            double saizu = TMath::Log10(t.Pt()/0.01);
            Color_t col;
            Style_t sty;

            if ( fl == 10 ) {
                col = kBlue;
                sty = kCircle;
            } else if ( fl == 11 ) {
                col = kAzure + 7;
                sty = kCircle;
            } else if ( fl == 12 ) {
                col = kRed;
                sty = kOpenSquare;
            } else if ( fl == 13 ) {
                col = kGreen-7;
                sty = kCircle;
            } else if ( fl == 14 ) {
                col = kYellow + 1;
                sty = kCircle;
            } else if ( fl == 15 ) {
                col = kGreen+2;
                sty = kOpenSquare;
            //} else if ( fl == 7 || fl == 8 ) {
            //    col = kBlack;
            //    sty = kOpenSquare;
            } else if ( fl == 9 ) {
                col = kMagenta;
                sty = kDiamond;
            } else {
                continue;
            }

            if (saizu < 0)
                continue;

            if ( fl == 13 || fl == 14 || fl == 15 || fl == 9 )
                saizu *= 2;

            TGraph *f = new TGraph(1);
            f->SetPoint(0,t.Rapidity(),phif(t.Phi(),phi_sum));
            f->SetMarkerStyle(sty);
            f->SetMarkerColor(col);
            f->SetMarkerSize(saizu);
            f->Draw("sameP");
        }
        break;
   }

}
