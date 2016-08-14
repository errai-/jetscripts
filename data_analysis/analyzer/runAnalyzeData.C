#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <string>
#include <vector>
#include "AnalyzeData.C"

using std::string;
using std::vector;
using std::cout;
using std::endl;

/* Indices:
    0   - pythia8 generated data
    1   - cmssw MC (juska)
    2   - cmssw MC (old)
    3-10 - real detector data */
const int fileAmount = 11;
//string fileArray[fileAmount] = {
//    "pythia8.root",
//    "ProcessedTree_RDMC_START53_V7N_try2.root",
//    "27thJune_MC_NTuple.root",
//    "29thAugust_Run2012A.root",
//    "29thAugust_Run2012B_JetHT.root",
//    "29thAugust_Run2012B_JetMon.root",
//    "29thAugust_Run2012C_JetHT.root",
//    "29thAugust_Run2012C_JetMon.root",
//    "29thAugust_Run2012D_JetHT.root",
//    "29thAugust_Run2012D_JetMon.root",
//    "19thNovember2014_Run2012B_JetHT.root"};

void runAnalyzeData(Long64_t events=0, Int_t firstIdx=0, Int_t lastIdx = 0,
    string writeFile="test.root", Int_t testing = 0, string dataPath="data/")
{
    if (events==0) {
        cout << "Usage: runAnalyzeData( [Number of events], [First file index], ";
        cout << "[Last file index], [File to write], [Testing flag], [Data path]";
        cout << endl << "File indices:" << endl;
        //for ( int i = 0; i < fileAmount; ++i ) {
        //    cout << i << " - " << fileArray[i] << endl;
        //}
        return;
    }

    if (firstIdx<0 || lastIdx>=fileAmount || firstIdx>lastIdx) {
        cout << "Bad file indices" << endl;
        return;
    }
    /* Initialization (no need to change) */
    gROOT->Reset();
    gROOT->Clear();
    gROOT->ProcessLine(".L QCDModules/QCDEvent.cc+");
    gROOT->ProcessLine(".L QCDModules/QCDJet.cc+");
    gROOT->ProcessLine(".L QCDModules/QCDMET.cc+");
    gROOT->ProcessLine(".L QCDModules/QCDCaloJet.cc+");
    gROOT->ProcessLine(".L QCDModules/QCDPFJet.cc+");
    gROOT->ProcessLine(".L QCDModules/QCDEventHdr.cc+");
    
    TChain *forest = new TChain("ak5/ProcessedTree");
    //TChain *forest = new TChain("JetTree");

    /* isDT is complementary to isMC, but isMC has many options depending on
       the data file in use */
    int isMC = 0;
    if (firstIdx == 1 ) { isMC = 1; } 
    else if ( firstIdx == 2) { isMC = 2; } 
    else if (firstIdx == 0) { isMC = 3; }

    //cout << "First file read: " << fileArray[firstIdx] << endl;
    //cout << "Last file read: " << fileArray[lastIdx] << endl;

    //string tmpPath;
    ///* Read consequent files from the list */
    //for (Int_t fileChoice = firstIdx; fileChoice<=lastIdx; ++fileChoice) {
    //    tmpPath = dataPath;
    //    tmpPath += fileArray[fileChoice];
    //    forest->AddFile(tmpPath.c_str());
    forest->AddFile("data/19thNovember2014_Run2012B_JetHT.root");
    //}
    
    AnalyzeData scriptObject(forest, events, isMC, testing);
    scriptObject.Loop(writeFile);
}

