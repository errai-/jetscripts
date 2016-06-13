#ifndef ProfileProjector_h
#define ProfileProjector_h

#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TKey.h>

#include "../analyzer/HistScripts.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::stringstream;

// A class for reading stack plots from a .root file
class ProfileProjector
{
private:
    TProfile2D *chfPu;
    TProfile2D *chf;
    TProfile2D *phf;
    TProfile2D *nhf;
    TProfile2D *elfMuf;
    TProfile2D *hfHf;
    TProfile2D *hfPhf;
    TH2D *ptBins;
    TH2D *ptBinsScaled;
    TFile *reading;

    vector<TH1D*> pileUp;

    TH1D *profileSlice;

    TH1D *chfPuSlice;
    TH1D *chfSlice;
    TH1D *phfSlice;
    TH1D *nhfSlice;
    TH1D *elfMufSlice;
    TH1D *hfHfSlice;
    TH1D *hfPhfSlice;
    THStack *eHistStack;

    TLegend *leg;

    string uniqueId;

public:
    ProfileProjector(string readFile, string id="default")
    {
        reading = new TFile(readFile.c_str());
        if ( ! reading->IsOpen() ) cout << "The file was not opened correctly" << endl;
        chfPu=(TProfile2D*)reading->Get("chfpu2D");
        chf=(TProfile2D*)reading->Get("chf2D");
        phf=(TProfile2D*)reading->Get("phf2D");
        nhf=(TProfile2D*)reading->Get("nhf2D");
        elfMuf=(TProfile2D*)reading->Get("elfmuf2D");
        hfHf=(TProfile2D*)reading->Get("hfhf2D");
        hfPhf=(TProfile2D*)reading->Get("hfphf2D");
        ptBins=(TH2D*)reading->Get("ptbins2D");
        ptBinsScaled=(TH2D*)reading->Get("ptbinsscaled2D");

        pileUp.push_back( (TH1D*)reading->Get("pileuphist1") );
        pileUp.push_back( (TH1D*)reading->Get("pileuphist2") );
        pileUp.push_back( (TH1D*)reading->Get("pileuphist3") );
        pileUp.push_back( (TH1D*)reading->Get("pileuphist4") );
        pileUp.push_back( (TH1D*)reading->Get("pileuphist5") );
        pileUp.push_back( (TH1D*)reading->Get("pileuphist6") );

        leg = new TLegend(0.4,0.2,0.6,0.35);
        uniqueId = id;
    }

    ~ProfileProjector(){
        delete chfPu;
        delete chf;
        delete phf;
        delete nhf;
        delete elfMuf;
        delete hfHf;
        delete hfPhf;
        delete ptBins;
        delete ptBinsScaled;
    }

    TProfile2D *GetChfPu(){ return chfPu; }
    TProfile2D *GetChf(){ return chf; }
    TProfile2D *GetPhf(){ return phf; }
    TProfile2D *GetNhf(){ return nhf; }
    TProfile2D *GetElfMuf(){ return elfMuf; }
    TProfile2D *GetHfHf(){ return hfHf; }
    TProfile2D *GetHfPhf(){ return hfPhf; }
    TLegend *GetLeg(){ return leg; }
    TH2D *GetPtBins(){ return ptBins; }
    TH2D *GetPtBinsScaled(){ return ptBinsScaled; }
    vector<TH1D*> GetPileUp(){return pileUp;}

    /* A function which gives a plot of pt or eta distribution of the jets for a 
     * given range of eta or pt */
    TH1D *PtProfile( int ptNotEta, double firstBin, double lastBin,
        string profileTitle="Energy profile"){
        
        int x=0,y=0,z=0;
        int lowerBin, upperBin;
        if (ptNotEta){
        if (firstBin==0 && lastBin==-1){
            lowerBin=0; upperBin=-1;
        } else {
            ptBins->GetBinXYZ( ptBins->FindBin(0,firstBin), x,y,z );
            lowerBin = y;
            ptBins->GetBinXYZ( ptBins->FindBin(0,lastBin), x,y,z );
            upperBin = y;
        }
        profileSlice = ptBins->ProjectionX( "ptnobins",lowerBin, upperBin );
        profileSlice->GetXaxis()->SetTitle( "p_{T} (GeV)" );
        } else {
        if (firstBin==0 && lastBin==-1){
            lowerBin=0; upperBin=-1;
        } else {
            ptBins->GetBinXYZ( ptBins->FindBin(firstBin), x,y,z );
            lowerBin = x;
            ptBins->GetBinXYZ( ptBins->FindBin(lastBin), x,y,z );
            upperBin = x;
        }
        profileSlice = ptBins->ProjectionY( "etanobins",lowerBin, upperBin );
        profileSlice->GetXaxis()->SetTitle( "#eta" );
        }
        profileSlice->SetTitle( profileTitle.c_str() );
        profileSlice->GetYaxis()->SetTitle( "Number of jets" );

        return (TH1D*) profileSlice->Clone();
    }

    /* A function which gives a plot of pt or eta distribution of the jets for a 
     * given range of eta or pt */
    TH1D *PtProfileScaled( int ptNotEta, double firstBin, double lastBin,
        string profileTitle="Energy profile"){
        
        int x=0,y=0,z=0;
        int lowerBin, upperBin;
        if (ptNotEta){
        if (firstBin==0 && lastBin==-1){
            lowerBin=0; upperBin=-1;
        } else {
            ptBinsScaled->GetBinXYZ( ptBinsScaled->FindBin(0,firstBin), x,y,z );
            lowerBin = y;
            ptBinsScaled->GetBinXYZ( ptBinsScaled->FindBin(0,lastBin), x,y,z );
            upperBin = y;
        }
        profileSlice = ptBinsScaled->ProjectionX( "ptnobins",lowerBin, upperBin );
        profileSlice->GetXaxis()->SetTitle( "p_{T} (GeV)" );
        } else {
        if (firstBin==0 && lastBin==-1){
            lowerBin=0; upperBin=-1;
        } else {
            ptBinsScaled->GetBinXYZ( ptBinsScaled->FindBin(firstBin), x,y,z );
            lowerBin = x;
            ptBinsScaled->GetBinXYZ( ptBinsScaled->FindBin(lastBin), x,y,z );
            upperBin = x;
        }
        profileSlice = ptBinsScaled->ProjectionY( "etanobins",lowerBin, upperBin );
        profileSlice->GetXaxis()->SetTitle( "#eta" );
        }
        profileSlice->SetTitle( profileTitle.c_str() );
        profileSlice->GetYaxis()->SetTitle( "Number of jets" );

        return (TH1D*) profileSlice->Clone();
    }

    void LegendFill( int forward ){
        if (forward){
            leg->AddEntry( hfPhfSlice, "Forward photons", "pf" );
            leg->AddEntry( hfHfSlice, "Forward hadrons", "pf" );
        }
        leg->AddEntry( elfMufSlice, "Electrons & Muons", "pf" );
        leg->AddEntry( nhfSlice, "Neutral hadrons", "pf" );
        leg->AddEntry( phfSlice, "Photons", "pf" );
        leg->AddEntry( chfSlice, "Charged hadrons", "pf" );
        leg->AddEntry( chfPuSlice, "Charged pile-up", "pf" );
    }

    /* Does the desired projection; a bit clumsy in order to avoid repetitive names */
    void Projector( int lowerBin, int upperBin, int ptNotEta, int errors)
    {
        stringstream chfp1, chfp2, chf1, chf2, phf1, phf2, nhf1, nhf2, elfMuf1,
        elfMuf2, hfHf1, hfHf2, hfPhf1, hfPhf2;
        chfp1 << uniqueId << "chfp"; chfp2 << uniqueId << "chfpuslice";
        chf1 << uniqueId << "chf"; chf2 << uniqueId << "chfslice";
        phf1 << uniqueId << "phf"; phf2 << uniqueId << "phfslice";
        nhf1 << uniqueId << "nhf"; nhf2 << uniqueId << "nhfslice";
        elfMuf1 << uniqueId << "elfmuf"; nhf2 << uniqueId << "elfmufslice";
        hfHf1 << uniqueId << "hfhf"; hfHf2 << uniqueId << "hfhfslice";
        hfPhf1 << uniqueId << "hfphf"; hfPhf2 << uniqueId << "hfphfslice";

        chfPuSlice = ptNotEta ? chfPu->ProfileX(chfp1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(chfp2.str().c_str(), errors ?"E":"") : chfPu->ProfileY(chfp1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(chfp2.str().c_str(), errors ?"E":"");
        chfSlice = ptNotEta ? chf->ProfileX(chf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(chf2.str().c_str(), errors ?"E":"") : chf->ProfileY(chf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(chf2.str().c_str(), errors ?"E":"");
        phfSlice = ptNotEta ? phf->ProfileX(phf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(phf2.str().c_str(), errors ?"E":"") : phf->ProfileY(phf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(phf2.str().c_str(), errors ?"E":"");
        nhfSlice = ptNotEta ? nhf->ProfileX(nhf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(nhf2.str().c_str(), errors ?"E":"") : nhf->ProfileY(nhf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(nhf2.str().c_str(), errors ?"E":"");
        elfMufSlice = ptNotEta ? elfMuf->ProfileX(elfMuf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(elfMuf2.str().c_str(), errors ?"E":"") : elfMuf->ProfileY(elfMuf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(elfMuf2.str().c_str(), errors ?"E":"");
        hfHfSlice = ptNotEta ? hfHf->ProfileX(hfHf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(hfHf2.str().c_str(), errors ?"E":"") : hfHf->ProfileY(hfHf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(hfHf2.str().c_str(), errors ?"E":"");
        hfPhfSlice = ptNotEta ? hfPhf->ProfileX(hfPhf1.str().c_str(),lowerBin,upperBin)->
        ProjectionX(hfPhf2.str().c_str(), errors ?"E":"") : hfPhf->ProfileY(hfPhf1.str().c_str(),
        lowerBin,upperBin)->ProjectionX(hfPhf2.str().c_str(), errors ?"E":"");
    }

  // A function that gives a stack plot of energy fractions
    THStack *EnergyFractions( int ptNotEta, double firstBin, double lastBin,
        int errors, int forwardInclude, string stackName="Energy stacks")
    {
        eHistStack = new THStack("energystacks", stackName.c_str() );
        int x=0,y=0,z=0;
        int lowerBin, upperBin;

        if (firstBin==0 && lastBin==-1){
            lowerBin=0; upperBin=-1;
        } else {
            ptBins->GetBinXYZ( ptBins->FindBin(ptNotEta ? 0 : firstBin, 
                ptNotEta ? firstBin : 0), x,y,z );
            lowerBin = ptNotEta ? y : x;
            ptBins->GetBinXYZ( ptBins->FindBin(ptNotEta ? 0 : lastBin,
                ptNotEta ? lastBin : 0), x,y,z );
            upperBin = ptNotEta ? y : x;
        }

        Projector( lowerBin, upperBin, ptNotEta, errors);

        if (ptNotEta){
            chfPuSlice->GetXaxis()->SetTitle( "p_{T} (GeV)" );
        } else {
            chfPuSlice->GetXaxis()->SetTitle( "#eta" );
        }
        
        chfPuSlice->SetTitle( stackName.c_str() );
        chfPuSlice->GetXaxis()->SetNoExponent();
        chfPuSlice->GetXaxis()->SetMoreLogLabels();
        chfPuSlice->GetYaxis()->SetTitle( "Energy fractions" );

        TH1D* h = (TH1D*) chfPuSlice;
        h->SetMaximum(0.95);
        eHistStack->SetHistogram( h );
        chfPuSlice->SetFillColor(kYellow);
        chfSlice->SetFillColor(kRed-7);
        phfSlice->SetFillColor(kBlue-7);
        nhfSlice->SetFillColor(kGreen-7);
        elfMufSlice->SetFillColor(kYellow+1);
        if (forwardInclude){
            hfHfSlice->SetFillColor(kGreen+7);
            hfPhfSlice->SetFillColor(kBlue+7);
            hfHfSlice->SetMarkerStyle(24);
            hfPhfSlice->SetMarkerStyle(25);
        }
        chfPuSlice->SetMarkerStyle(26);
        chfSlice->SetMarkerStyle(27);
        phfSlice->SetMarkerStyle(28);
        nhfSlice->SetMarkerStyle(31);
        elfMufSlice->SetMarkerStyle(32);
        
        eHistStack->Add( (TH1D*) chfPuSlice );
        eHistStack->Add( (TH1D*) chfSlice );
        eHistStack->Add( (TH1D*) phfSlice );
        eHistStack->Add( (TH1D*) nhfSlice );
        eHistStack->Add( (TH1D*) elfMufSlice );

        LegendFill( forwardInclude );
        if (forwardInclude){
            eHistStack->Add( (TH1D*) hfHfSlice );
            eHistStack->Add( (TH1D*) hfPhfSlice );
        }

        return (THStack*) eHistStack->Clone();
    }
};

#endif
