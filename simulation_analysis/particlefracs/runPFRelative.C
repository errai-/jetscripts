// plots jet energy fractions, see README_ScriptInfo

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// ROOT, for histogramming.
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TColor.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TLatex.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod14_pythia8.C"


void stackModify(TH1D *setter){
  setter->GetXaxis()->SetTitle("p_{T} (GeV)");
  setter->GetYaxis()->SetTitle("Energy fraction difference (%)");
  setter->SetStats(0);
  setter->GetXaxis()->SetMoreLogLabels();
  setter->GetXaxis()->SetNoExponent();
  setter->GetYaxis()->SetTitleOffset(1.25);
  setter->GetXaxis()->SetTitleOffset(1.0);
  setter->GetXaxis()->SetLabelSize(0.045);
  setter->GetYaxis()->SetLabelSize(0.045);
  setter->GetXaxis()->SetTitleSize(0.045);
  setter->GetYaxis()->SetTitleSize(0.045);
}

void WorkHorse(vector<TH1D*>& PFfracs, std::string fileName) {
    TFile *inFile = new TFile(fileName.c_str(), "READ");
    // Read histograms.
    vector<TProfile*> fractionProfiles;
    for (int idx = 0; idx != 16; ++idx){
        std::stringstream tmpString("");
        tmpString << "a" << idx; 
        // a = all, g = gluonjets, q = quarkjets, lq = low pt quarks, hq = high pt quarks
        fractionProfiles.push_back( (TProfile*) inFile->Get( tmpString.str().c_str()) );
    }
    cout << fractionProfiles.size();

    vector<TH1D*> fractionHists;
    for (int idx = 0; idx != 16; ++idx){
        fractionHists.push_back( fractionProfiles[idx]->ProjectionX("","") );
    }
    PFfracs.push_back( (TH1D*) fractionHists[0]->Clone("chf") );
    PFfracs[0]->Add( fractionHists[1] );
    PFfracs[0]->Add( fractionHists[3] );
    PFfracs[0]->Add( fractionHists[4] );
    PFfracs[0]->Add( fractionHists[7] );
    PFfracs[0]->Add( fractionHists[8] );
    PFfracs[0]->Add( fractionHists[13] );
    PFfracs.push_back( (TH1D*) fractionHists[2]->Clone("phf") );
    PFfracs[1]->Add( fractionHists[11] );
    PFfracs.push_back( (TH1D*) fractionHists[5]->Clone("nhf") );
    PFfracs[2]->Add( fractionHists[6] );
    PFfracs[2]->Add( fractionHists[9] );
    PFfracs[2]->Add( fractionHists[10] );
    PFfracs[2]->Add( fractionHists[12] );
    PFfracs.push_back( (TH1D*) fractionHists[14]->Clone("elf+muf") );
    PFfracs.push_back( (TH1D*) fractionHists[15]->Clone("others") );

    PFfracs[0]->SetFillColor(kRed);
    PFfracs[1]->SetFillColor(kBlue);
    PFfracs[2]->SetFillColor(kGreen);
    PFfracs[3]->SetFillColor(kYellow);
    PFfracs[4]->SetFillColor(kCyan);
}

int runPFRelative(std::string fileName1, std::string fileName2) {

    // Create file on which histogram(s) can be saved.
    // THStack* partFracs = new THStack("particlestacks", "");
    vector<TH1D*> fracs1, fracs2;
    WorkHorse(fracs1,fileName1);
    WorkHorse(fracs2,fileName2);

    TH1D* dummy = (TH1D*) fracs1[0]->Clone();
    dummy->Divide( fracs1[0] );

    stackModify(fracs1[0]);
    TCanvas *canv = tdrCanvas("c1",fracs1[0],12,0,1);//new TCanvas("c1","c1",600,600);
    fracs1[0]->SetMarkerColor(kRed);
    fracs1[1]->SetMarkerColor(kBlue);
    fracs1[2]->SetMarkerColor(kGreen);
    fracs1[3]->SetMarkerColor(kYellow);
    fracs1[4]->SetMarkerColor(kCyan);
    for (unsigned i = 0; i < fracs1.size(); ++i) {
        fracs1[i]->Add( fracs2[i], -1 );
        //fracs1[i]->Divide( fracs2[i] );
        //fracs1[i]->Add( dummy, -1 );
        fracs1[i]->Scale(100);
        fracs1[i]->SetMarkerStyle(20);
        fracs1[i]->SetMaximum( 2.5 );
        fracs1[i]->SetMinimum( -1.5 );
        fracs1[i]->Draw("sameP");
    }


    TLegend *leg = tdrLeg(0.83,0.1,0.98,0.93);
    leg->AddEntry( fracs1[4], "rest", "f" );
    leg->AddEntry( fracs1[3], "elf+#muf", "f" );
    leg->AddEntry( fracs1[2], "nhf", "f" );
    leg->AddEntry( fracs1[1], "phf", "f" );
    leg->AddEntry( fracs1[0], "chf", "f" );
    leg->SetTextSize(0.045);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);

    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.6*0.07);

    latex.DrawLatex(0.7,0.88,"Pythia 6 Z2* - Pythia8 CUETP8S1");
//
//    fixOverlay();
//    canv->Print("efracs.pdf"); 
//    // Done.
    return 0;
}
