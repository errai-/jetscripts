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

const double radius2 = TMath::Power(0.5,2);
const bool lines = true;
const double pi2 = 2*TMath::Pi();

bool accept() {
    string ses;
    cin >> ses;
    if ( ses == "y" ) return true;
    else return false;
}

double phif(double val, double center) {
    double corr = val - center;
    while ( corr < -TMath::Pi() )
        corr += pi2;
    while ( corr > TMath::Pi() )
        corr -= pi2;
    return corr;
}

void Append(TH2D* th, double eta, double phi, double pt) {
    for (int phi_ind = 1, Np = th->GetNbinsY(); phi_ind <= Np; ++phi_ind) {
        double phi_pos = th->GetYaxis()->GetBinCenter(phi_ind);
        double phi_diff = fabs( phi - phi_pos );
        if ( phi_diff > TMath::Pi() )
            phi_diff = pi2 - phi_diff;
        for (int eta_ind = 1, Ne = th->GetNbinsX(); eta_ind <= Ne; ++eta_ind) {
            double eta_pos = th->GetXaxis()->GetBinCenter(eta_ind);
            double dist2 = TMath::Power(phi_diff, 2)+TMath::Power(eta - eta_pos, 2);

            th->Fill(eta_pos,phi_pos,pt*TMath::Exp(-dist2/radius2));
        }
    }
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
        double phi_sum = 0, phi_diff = 0, pt_sum = 0;
        vector<int> jets;

        if ( choice >= 0 && choice != jentry )
            continue;

        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
            if (fl < 0 || fl == 10)
                continue;
            jets.push_back(j);

            if ( fl != 0 )
                ++count;

            if ( fabs(t.Rapidity()) > 3 && jets.size() < 3 ) {
                stahp = true;
                break;
            }

            if ( jets.size() < 3 ) {
                phi_sum += t.Phi()*t.Pt();
                pt_sum += t.Pt();
                phi_diff += t.Phi()*TMath::Power(-1.0,int(jets.size()));
            }
        }
        if ( /*stahp ||*/ (count > 1 && choice == -1) ) {
            continue;
        }
        cout << jentry << " " << count << endl;
        if ( !accept() )
            continue;
        phi_sum /= pt_sum;
        if ( fabs(phi_diff) < TMath::Pi() )
            phi_sum += TMath::Pi();

        TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
        c1->SetLogz(1);
        TH2D* grid = new TH2D("","",400,-5,5,400,-TMath::Pi(),TMath::Pi());
        grid->SetMaximum(200);
        grid->SetMinimum(1);

        // Give all bins a nice color in log plot
        for (int eta_ind = 1, Ne = grid->GetNbinsX(); eta_ind <= Ne; ++eta_ind) {
            double eta_pos = grid->GetXaxis()->GetBinCenter(eta_ind);
            for (int phi_ind = 1, Np = grid->GetNbinsY(); phi_ind <= Np; ++phi_ind) {
                double phi_pos = grid->GetYaxis()->GetBinCenter(phi_ind);
                grid->Fill(eta_pos,phi_pos,1);
            }
        }

        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);

            if (fl == 10)
                continue;

            fl *= -1;

            if ( fl == 7 || fl == 8 ) { // Final state
                Append( grid, t.Rapidity(), phif(t.Phi(),phi_sum), t.Pt() );
            }
        }
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        grid->SetXTitle("y");
        grid->SetYTitle("#phi");
        grid->Draw("COLZ");

        for ( unsigned j = 0; j < jets.size(); ++j ) { 
            TLorentzVector t(fX[jets[j]],fY[jets[j]],fZ[jets[j]],fT[jets[j]]);
            TEllipse *el = new TEllipse(t.Rapidity(),phif(t.Phi(),phi_sum),0.5,0.5);
            Color_t ell_col = kRed+1;
            if (j == 1)
                ell_col = kRed;
            else if (j == 2)
                ell_col = kMagenta;
            else if (j == 3)
                ell_col = kCyan;
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
                col = kCyan;
                sty = kCircle;
            } else if ( fl == 12 ) {
                col = kRed;
                sty = kOpenSquare;
            } else if ( fl == 13 ) {
                col = kBlue+2;
                sty = kCircle;
            } else if ( fl == 14 ) {
                col = kCyan+2;
                sty = kCircle;
            } else if ( fl == 15 ) {
                col = kGreen+3;
                sty = kOpenSquare;
            //} else if ( fl == 7 || fl == 8 ) {
            //    col = kBlack;
            //    sty = kOpenSquare;
            } else if ( fl == 9 ) {
                col = kMagenta-3;
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

        if (lines) {
            TLine *line1 = new TLine(-5,0,5,0);
            line1->SetLineStyle(9);
            line1->Draw("same");

            TLine *line2 = new TLine(-5,-TMath::Pi()/2,5,-TMath::Pi()/2);
            line2->SetLineStyle(2);
            line2->SetLineColor(8);
            line2->Draw("same");

            TLine *line3 = new TLine(-5,TMath::Pi()/2,5,TMath::Pi()/2);
            line3->SetLineStyle(2);
            line3->SetLineColor(8);
            line3->Draw("same");
        }

        break;
   }

}
