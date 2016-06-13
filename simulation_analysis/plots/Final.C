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

void EtaPhi::QuarkPairs( vector<pair<int,int> >& smaller, vector<pair<int,int> >& larger ) {
    return;
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
        double pt_sum = 0;
        double phi_diff = 0;

        if ( choice >= 0 && choice != jentry )
            continue;

        cout << jentry << endl;
        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            if (fl < 0 || (fl >= 10 && fl < 20) )
                continue;
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
            jets.push_back(j);

            cout << "  " << fl << endl;
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

        if ( choice == -1 && ( stahp ) ) //(count > 1 && choice == -1) ) {
            continue;

        if ( !accept() ) {
            choice = -1;
            continue;
        }

        phi_sum /= pt_sum;
        if ( fabs(phi_diff) < TMath::Pi() )
            phi_sum += TMath::Pi();

        TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1200);
        TH1D* th = new TH1D("","",100,-5.0,5.0);
        gStyle->SetOptStat(0);
        th->SetMinimum(-TMath::Pi());
        th->SetMaximum((3/2.0)*TMath::Pi());
        th->SetXTitle("y");
        th->SetYTitle("#phi");
        th->Draw("P");

        for ( unsigned j = 0; j < jets.size(); ++j ) {
            TLorentzVector t(fX[jets[j]],fY[jets[j]],fZ[jets[j]],fT[jets[j]]);
            double loc_phi = phif(t.Phi(),phi_sum);
            TEllipse *el = new TEllipse(t.Rapidity(),loc_phi,0.5,0.5);
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

        vector<int> counters(5,0);
        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];

            if (fl >= 0)
                continue;

            fl *= -1;
            fl -= 1;

            int group = fConstituents[j];
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);

            double saizu = TMath::Log10(t.Pt()/0.01);
            Color_t col;
            Style_t sty;

            if ( group == 3 ) {
                if ( fl < 6 ) {
                    col = kBlue+2;
                    sty = kCircle;
                } else {
                    col = kGreen+3;
                    sty = kOpenSquare;
                }
            } else if ( group == 4 ) {
                col = kCyan;
                sty = kCircle;
            } else if ( group == 5 ) {
                col = kMagenta-3;
                sty = kDiamond;
            } else if ( group == 2 ) {
                if ( fl == 0 ) {
                    col = kRed;
                    sty = kOpenSquare;
                } else if ( fl == 1 ) {
                    col = kBlue;
                    sty = kOpenCircle;
                }
            } else if ( group == 9 ) {
                col = kRed;
                sty = kOpenSquare;
            } else {
                continue;
            }

            if ( group == 3 || group == 4 || group == 5 )
                saizu *= 2;

            if (saizu < 0)
                continue;

            if (t.Rapidity() > 4 && t.Rapidity() < 5 && t.Phi() < -1)
                cout << t.Rapidity() << " " << t.Phi() << " " << fl << " " << saizu << endl;

            TGraph *f = new TGraph(1);
            double loc_phi = phif(t.Phi(),phi_sum);
            f->SetPoint(0,t.Rapidity(),loc_phi);
            f->SetMarkerStyle(sty);
            f->SetMarkerColor(col);
            f->SetMarkerSize(saizu);
            f->Draw("sameP");
            if ( loc_phi < -TMath::Pi()/2 ) {
                TGraph *g = new TGraph(1);
                g->SetPoint(0,t.Rapidity(),loc_phi+2*TMath::Pi());
                g->SetMarkerStyle(sty);
                g->SetMarkerColor(col);
                g->SetMarkerSize(saizu);
                g->Draw("sameP");
            }
        }
        break;
   }

}
