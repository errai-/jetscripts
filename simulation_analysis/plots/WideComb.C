#define EtaPhi_cxx
#include "etaphi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::map;
using std::set;

const double radius2 = TMath::Power(0.5,2);
const bool lines = true;
const bool connect = false;
const bool useflav = true;
const bool printflav = true;
const bool printremn = true;
const bool skip = false;
const double pi2 = 2*TMath::Pi();
const double separation2 = TMath::Power(25,2);
const unsigned stop = 600;

bool accept() {
    cout << "OK? (y/n) ";
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

double dphi(double phi1, double phi2) {
    double diff = fabs( phi1 - phi2 );
    if ( diff > TMath::Pi() )
        diff = pi2 - diff;
    return diff;
}

double dR2(double phi1, double phi2, double eta1, double eta2) {
    double phi_diff = dphi(phi1,phi2);
    return TMath::Power(phi_diff,2)+TMath::Power(eta1-eta2,2);
}

void Append(TH2D* th, double eta, double phi, double pt) {
    for (int phi_ind = 1, Np = th->GetNbinsY(); phi_ind <= Np; ++phi_ind) {
        double phi_pos = th->GetYaxis()->GetBinCenter(phi_ind);
        double phi_diff = dphi( phi, phif(phi_pos,0) );
        for (int eta_ind = 1, Ne = th->GetNbinsX(); eta_ind <= Ne; ++eta_ind) {
            double eta_pos = th->GetXaxis()->GetBinCenter(eta_ind);
            double dist2 = TMath::Power(phi_diff, 2)+TMath::Power(eta - eta_pos, 2);

            th->Fill(eta_pos,phi_pos,pt*TMath::Exp(-dist2/radius2));
        }
    }
}

void EtaPhi::QuarkPairs( vector<pair<int,int> >& smaller, vector<pair<int,int> >& larger ) {
    map< double, pair<int,int> > dists;
    for (unsigned i = 0, Nx = smaller.size(); i < Nx; ++i) {
        int idx = smaller[i].first;
        int iflav = smaller[i].second;
        TLorentzVector s(fX[idx],fY[idx],fZ[idx],fT[idx]);
        for (unsigned j = 0, Ny = larger.size(); j < Ny; ++j) {
            int jdx = larger[j].first;
            int jflav = larger[j].second;
            TLorentzVector t(fX[jdx],fY[jdx],fZ[jdx],fT[jdx]);
            if ( !useflav || (iflav == jflav) ) {
                double phi_diff = dphi( s.Phi(), t.Phi() );
                double dist2 = dR2(s.Phi(),t.Phi(),s.Eta(),t.Eta());
                if ( dist2 > separation2 )
                    continue;
                //dist2 *= min( t.Perp2(), s.Perp2() );
                dists.insert( std::make_pair( dist2, std::make_pair(idx,jdx) ) );
            }
        }
    }

    set<int> used;
    for ( auto it = dists.begin(); it != dists.end(); ++it ) {
        int idx = it->second.first;
        int jdx = it->second.second;
        if ( used.find(idx) == used.end() && used.find(jdx) == used.end() ) {
            TLorentzVector s(fX[idx],fY[idx],fZ[idx],fT[idx]);
            TLorentzVector t(fX[jdx],fY[jdx],fZ[jdx],fT[jdx]);
            if ( fabs(s.Eta()) > 5 && fabs(t.Eta()) > 5 )
                continue;
            double phi1 = phif(s.Phi(),phi_sum), phi2 = phif(t.Phi(),phi_sum);
            if ( phi1 < -TMath::Pi()/2 && phi2 > TMath::Pi()/2 )
                phi1 += pi2;
            else if ( phi2 < -TMath::Pi()/2 && phi1 > TMath::Pi()/2 )
                phi2 += pi2;
            TLine *line = new TLine(s.Rapidity(),phi1,t.Rapidity(),phi2);
            line->SetLineWidth(2);
            line->Draw("same");

            used.insert( it->second.first );
            used.insert( it->second.second );
        }
    }
}

void EtaPhi::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    std::ofstream ulos("sus.txt");
    Long64_t nbytes = 0, nb = 0;
    cout << "Choose event, -1 for no choice" << endl;
    int choice;
    cin >> choice;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);
        if ( stop > 0 && jentry > stop )
            break;

        bool stahp = false;
        int count = 0;
        phi_sum = 0;
        double phi_diff = 0, pt_sum = 0;
        vector<int> jets;

        if ( choice >= 0 && choice != jentry )
            continue;

        if (printflav)
            cout << jentry << endl;
        for ( int j = 0; j < fJets; ++j ) {
            int fl = fFlav[j];
            if (fl < 0 || (fl >= 10 && fl < 20) )
                continue;
            TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
            jets.push_back(j);

            if (printflav)
                cout << " " << fl << " " << t.E() << endl;
            if ( fl != 0 )
                ++count;

            if ( choice==-1 && (fabs(t.Rapidity()) > 3 && jets.size() < 3) ) {
                stahp = true;
                break;
            }
            if ( jets.size() < 3 ) {
                ulos << jentry << " " << fFlav[j] << endl;
                phi_sum += t.Phi()*t.Pt();
                pt_sum += t.Pt();
                phi_diff += t.Phi()*TMath::Power(-1.0,int(jets.size()));
            }
        }
        if ( choice == -1 && ( stahp ) ) // || count > 1 ) )
            continue;
        if ( skip || !accept() ) {
            choice = -1;
            continue;
        }

        phi_sum /= pt_sum;
        if ( fabs(phi_diff) < TMath::Pi() )
            phi_sum += TMath::Pi();

        TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1200);
        c1->SetLogz(1);
        TH2D* grid = new TH2D("","",400,-5,5,400,-TMath::Pi(),(3/2.0)*TMath::Pi());
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
            int group = fConstituents[j];

            if ( group == 2 ) { // Final state
                TLorentzVector t(fX[j],fY[j],fZ[j],fT[j]);
                Append( grid, t.Rapidity(), phif(t.Phi(),phi_sum), t.Pt() );
            }
        }
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        grid->SetXTitle("y");
        grid->SetYTitle("#phi (rad)");
        grid->Draw("COLZ");

        for ( unsigned j = 0; j < jets.size(); ++j ) {
            TLorentzVector t(fX[jets[j]],fY[jets[j]],fZ[jets[j]],fT[jets[j]]);
            double loc_phi = phif(t.Phi(),phi_sum);
            if (loc_phi < -3*TMath::Pi()/4)
                loc_phi += pi2;
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

        vector<pair<int,int> > quarks;
        vector<pair<int,int> > antiquarks;
        vector<TLorentzVector> x_quark, x_aquark;

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
                    x_quark.push_back(t);
                    col = kBlue+2;
                    sty = kCircle;
                } else {
                    col = kGreen+3;
                    sty = kOpenSquare;
                }
            } else if ( group == 4 ) {
                x_aquark.push_back(t);
                col = kCyan;
                sty = kCircle;
            } else if ( group == 5 ) {
                col = kMagenta-3;
                sty = kDiamond;
            } else if ( group == 7 ) {
                if ( fl < 6 ) {
                    col = kBlue;
                    sty = kCircle;
                    quarks.push_back( std::make_pair( j, fl ) );
                } else if ( fl > 2000 ) {
                    col = kBlue;
                    sty = kOpenSquare;
                    if (printremn)
                        cout << "pos " << t.Eta() << endl;
                } else {
                    col = kRed;
                    sty = kOpenSquare;
                }
            } else if ( group == 8 ) {
                col = kMagenta;
                if ( fl > 2000 ) {
                    sty = kOpenSquare;
                    if (printremn)
                        cout << "anti " << t.Eta() << endl;
                } else {
                    sty = kCircle;
                    antiquarks.push_back( std::make_pair( j, fl ) );
                }
            } else if ( group == 9 ) {
                col = kCyan+2;
                sty = kOpenSquare;
            } else {
                continue;
            }

            if (saizu < 0)
                continue;

            if ( (group > 2 && group < 6) || group == 9 )
                saizu *= 2;

            TGraph *f = new TGraph(1);
            double loc_phi = phif(t.Phi(),phi_sum);
            f->SetPoint(0,t.Rapidity(),loc_phi);
            f->SetMarkerStyle(sty);
            f->SetMarkerColor(col);
            f->SetMarkerSize(saizu);
            f->Draw("sameP");
            if ( loc_phi < -TMath::Pi()/2 ) {
                TGraph *g = new TGraph(1);
                g->SetPoint(0,t.Rapidity(),loc_phi+pi2);
                g->SetMarkerStyle(sty);
                g->SetMarkerColor(col);
                g->SetMarkerSize(saizu);
                g->Draw("sameP");
            }
        }

        if (lines) {
            TLine *line1 = new TLine(-5,0,5,0);
            line1->SetLineStyle(8);
            line1->SetLineColor(17);
            line1->Draw("same");

            TLine *line2 = new TLine(-5,-TMath::Pi()/2,5,-TMath::Pi()/2);
            line2->SetLineStyle(2);
            line2->SetLineColor(8);
            line2->Draw("same");

            TLine *line3 = new TLine(-5,TMath::Pi()/2,5,TMath::Pi()/2);
            line3->SetLineStyle(2);
            line3->SetLineColor(8);
            line3->Draw("same");

            TLine *line4 = new TLine(-5,TMath::Pi(),5,TMath::Pi());
            line4->SetLineStyle(8);
            line4->SetLineColor(17);
            line4->Draw("same");
        }

        if ( connect ) {
            for ( unsigned i = 0, Ni = x_quark.size(); i < Ni; ++i ) {
                int best_match = -1;
                double min_dist = 100;
                for ( unsigned j = 0, Nj = quarks.size(); j < Nj; ++j ) {
                    int id = quarks[j].first;
                    TLorentzVector t(fX[id],fY[id],fZ[id],fT[id]);
                    double dist = dR2(t.Phi(),x_quark[i].Phi(),t.Eta(),x_quark[i].Eta());
                    if ( dist < min_dist ) {
                        best_match = j;
                        min_dist = dist;
                    }
                }
                quarks.erase( quarks.begin()+best_match );
            }
            for ( unsigned i = 0, Ni = x_aquark.size(); i < Ni; ++i ) {
                int best_match = -1;
                double min_dist = 100;
                for ( unsigned j = 0, Nj = antiquarks.size(); j < Nj; ++j ) {
                    int id = antiquarks[j].first;
                    TLorentzVector t(fX[id],fY[id],fZ[id],fT[id]);
                    double dist = dR2(t.Phi(),x_aquark[i].Phi(),t.Eta(),x_aquark[i].Eta());
                    if ( dist < min_dist ) {
                        best_match = j;
                        min_dist = dist;
                    }
                }
                antiquarks.erase( antiquarks.begin()+best_match );
            }
            if (printremn)
                cout << "quarks: " << quarks.size() << " antiquarks: " << antiquarks.size() << endl;


            if ( antiquarks.size() < quarks.size() )
                QuarkPairs( antiquarks, quarks );
            else
                QuarkPairs( quarks, antiquarks );
        }

        break;
   }

}
