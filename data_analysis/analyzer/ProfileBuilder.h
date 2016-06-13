#ifndef ProfileBuilder_h
#define ProfileBuilder_h

#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>
#include <cassert>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <sstream>
#include "HistScripts.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

// A class for creating stack plots and storing them to a .root file
class ProfileBuilder
{
    private:

    TProfile2D *chfpu;
    TProfile2D *chf;
    TProfile2D *phf;
    TProfile2D *nhf;
    TProfile2D *elfMuf;
    TProfile2D *hfHf;
    TProfile2D *hfPhf;
    TH2D *ptBins;
    TH2D *ptBinsScaled;

    Long64_t ptBinAmnt;
    Long64_t etaBinAmnt;

    double minPt;
    double maxPt;
    double minEta;
    double maxEta;

    public:

    ProfileBuilder(Long64_t ptBinAmount, Long64_t etaBinAmount, double miPt, 
        double maPt, double miEta, double maEta){
        
        // Initializing some useful constants
        ptBinAmnt = ptBinAmount;
        etaBinAmnt = etaBinAmount;
        minPt = miPt;
        maxPt = maPt;
        minEta = miEta;
        maxEta = maEta;

        ptBinAmnt = 40;
        const double ptLimits[] = 
        //{1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 
            {56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 
            300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 
            790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 
            1588, 1684, 1784, 1890, 2000};
        //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037};
        // Slightly different boundaries:
        //{56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430, 
        //507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, 2000};
        //8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 2000, 2238, 2500, 2787, 3103, 3450};
        double logLimits[ptBinAmnt+1];
        double lattice = TMath::Log10( maPt/miPt )/ptBinAmnt;
        for (int idx = 0; idx != ptBinAmnt+1; ++idx){
            logLimits[ idx ] = miPt*TMath::Power( 10, idx*lattice );
        }

        chfpu = new TProfile2D("chfpu2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        chf = new TProfile2D("chf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        phf = new TProfile2D("phf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        nhf = new TProfile2D("nhf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        elfMuf = new TProfile2D("elfmuf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        hfHf = new TProfile2D("hfhf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        hfPhf = new TProfile2D("hfphf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        ptBins = new TH2D("ptbins2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
        ptBinsScaled = new TH2D("ptbinsscaled2D","",ptBinAmnt,logLimits,etaBinAmnt,miEta,maEta);
    }

    ~ProfileBuilder(){
        delete chfpu;
        delete chf;
        delete phf;
        delete nhf;
        delete elfMuf;
        delete hfHf;
        delete hfPhf;
    }

    TProfile2D *GetChfPu(){ return chfpu; }
    TProfile2D *GetChf(){ return chf; }
    TProfile2D *GetPhf(){ return phf; }
    TProfile2D *GetNhf(){ return nhf; }
    TProfile2D *GetElfMuf(){ return elfMuf; }
    TProfile2D *GetHfHf(){ return hfHf; }
    TProfile2D *GetHfPhf(){ return hfPhf; }
    TH2D *GetPtBins(){ return ptBins; }
    TH2D *GetPtBinsScaled(){ return ptBinsScaled; }

    void FillHelper( double _pt, double _eta, double _chfPu, double _chf, 
        double _phf, double _nhf, double _elfMuf, double _hfHf, double _hfPhf,
        double weight = 1., double prescale = 1.) 
    {
        chfpu->Fill( _pt, _eta, _chfPu, weight );
        chf->Fill( _pt, _eta, _chf, weight );
        phf->Fill( _pt, _eta, _phf, weight );
        nhf->Fill( _pt, _eta, _nhf, weight );
        elfMuf->Fill( _pt, _eta, _elfMuf, weight );
        hfHf->Fill( _pt, _eta, _hfHf, weight );
        hfPhf->Fill( _pt, _eta, _hfPhf, weight );
        ptBins->Fill( _pt, _eta, weight );
        ptBinsScaled->Fill( _pt, _eta, prescale*weight );
    }

    void WriteToFile(vector<TH1D*>* hists, string fileName) {
        TFile *fileToWrite = new TFile(fileName.c_str(),"RECREATE");
        
        chfpu->Write();
        chf->Write();
        phf->Write();
        nhf->Write();
        elfMuf->Write();
        hfHf->Write();
        hfPhf->Write();
        ptBins->Write();
        ptBinsScaled->Write();
        for (int i=0; i<6; i++){
            (*hists)[i]->Write();
        }
        fileToWrite->Close();
    }
};

#endif
