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

void Append(TH2D* th, double eta, double phi, double pt) {
    for (int phi_ind = 1, Np = th->GetNbinsY(); phi_ind <= Np; ++phi_ind) {
        double phi_pos = th->GetYaxis()->GetBinCenter(phi_ind);
        double phi_diff = fabs( phi - phi_pos );
        if ( phi_diff > TMath::Pi() )
            phi_diff = 2*TMath::Pi() - phi_diff;
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
        double phi_sum = 0, pt_sum = 0;
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
            phi_sum += t.Phi()*t.Pt();
            pt_sum += t.Pt();
        }
        if ( /*stahp ||*/ (count > 1 && choice == -1) ) {
            continue;
        }
        cout << jentry << " " << count << endl;
        if ( !accept() )
            continue;
        phi_sum /= pt_sum;

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
        break;
   }

}
