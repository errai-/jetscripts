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
  setter->GetYaxis()->SetTitle("Jet energy fraction");
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

int runPFSortedJets(std::string fileName) {

    // Create file on which histogram(s) can be saved.
    TFile *inFile = new TFile(fileName.c_str(), "READ");
    THStack* partFracs = new THStack("particlestacks", "");

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
    vector<TH1D*> PFfracs;
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

    for (unsigned i = 0; i != PFfracs.size(); ++i){
        partFracs->Add( PFfracs[i] );
    }

    TCanvas *canv = tdrCanvas("c1",fractionHists[0],12,0,1);//new TCanvas("c1","c1",600,600);

    stackModify(fractionHists[0]);
    partFracs->SetHistogram( fractionHists[0] );
    partFracs->GetHistogram()->SetMaximum(1.01);
    partFracs->GetHistogram()->SetMinimum(-0.01);
    //partFracs->SetFillStyle(4000);

    TLegend *leg = tdrLeg(0.83,0.1,0.98,0.93);
    leg->AddEntry( PFfracs[4], "rest", "f" );
    leg->AddEntry( PFfracs[3], "elf+muf", "f" );
    leg->AddEntry( PFfracs[2], "nhf", "f" );
    leg->AddEntry( PFfracs[1], "phf", "f" );
    leg->AddEntry( PFfracs[0], "chf", "f" );
    leg->SetTextSize(0.045);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);

    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.6*0.07);

    partFracs->Draw("same");
    latex.DrawLatex(0.77,0.85,"Pythia 6 Z2*");

    //fixOverlay();
    //canv->Print("efracs.pdf"); 
    // Done.
    return 0;
}
