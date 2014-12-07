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
#include <assert.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <JetScripts.h>
#include <TFile.h>

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
  TH2D *ptBins;

  Double_t *ptLimits;

  Long64_t ptBinAmnt;
  Long64_t etaBinAmnt;

  Double_t minPt;
  Double_t maxPt;
  Double_t minEta;
  Double_t maxEta;

  Double_t ptLattice;

  public:
  ProfileBuilder(Long64_t ptBinAmount, Long64_t etaBinAmount, Double_t miPt, 
    Double_t maPt, Double_t miEta, Double_t maEta){
    
    // Initializing some useful constants
    ptBinAmnt = ptBinAmount;
    etaBinAmnt = etaBinAmount;
    minPt = miPt;
    maxPt = maPt;
    minEta = miEta;
    maxEta = maEta;
    ptLattice = TMath::Log10( maxPt/minPt )/ptBinAmnt;

    ptLimits = new Double_t[ptBinAmnt];
    for (Long64_t idx = 0; idx <= ptBinAmnt; idx++){
      ptLimits[idx] = minPt*TMath::Power(10, idx*ptLattice);
    }
    chfpu = new TProfile2D("chfpu2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
    chf = new TProfile2D("chf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
    phf = new TProfile2D("phf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
    nhf = new TProfile2D("nhf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
    elfMuf = new TProfile2D("elfmuf2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
    ptBins = new TH2D("ptbins2D","",ptBinAmnt,ptLimits,etaBinAmnt,miEta,maEta);
  }

  ~ProfileBuilder(){
    delete chfpu;
    delete chf;
    delete phf;
    delete nhf;
    delete elfMuf;
    delete[] ptLimits;
  }

  TProfile2D *GetChfPu(){ return chfpu; }
  TProfile2D *GetChf(){ return chf; }
  TProfile2D *GetPhf(){ return phf; }
  TProfile2D *GetNhf(){ return nhf; }
  TProfile2D *GetElfMuf(){ return elfMuf; }
  TH2D *GetPtBins(){ return ptBins; }

  void AddValue( Int_t idx, Double_t pt, Double_t eta, Double_t fraction ){
    switch( idx ){
    case 0:
      chfpu->Fill( pt, eta, fraction );
      break;
    case 1:
      chf->Fill( pt, eta, fraction );
      break;
    case 2:
      phf->Fill( pt, eta, fraction );
      break;
    case 3:
      nhf->Fill( pt, eta, fraction );
      break;
    case 4:
      elfMuf->Fill( pt, eta, fraction );
      break;
    case 5:
      ptBins->Fill( pt, eta, fraction );
      break;
    }
  }

  void WriteToFile(string fileName="testfile.root") {
    TFile *fileToWrite = new TFile(fileName.c_str(),"NEW");
    if ( fileToWrite->IsOpen() ) {
      cout << "The file " << fileName << " was opened successfully." << endl;
    } else {
      cout << "The file " << fileName << " already exists, please try another "
        << "name." << endl;
      return;
    }
    chfpu->Write();
    chf->Write();
    phf->Write();
    nhf->Write();
    elfMuf->Write();
    ptBins->Write();
    fileToWrite->Close();
  }
};

#endif
