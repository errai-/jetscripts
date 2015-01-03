#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void runAnalyzeData(Long64_t loopLimit=10000, Int_t firstIdx=0, Int_t lastIdx = 0,
  Int_t testing = 0, string writeFile="test.root", 
  string dataPath="../../root_jetdata/")
{
  cout << "Processing runAnalyzeData.C ..." << endl;

  // These are not usually changed
  char scriptName[] = ".L AnalyzeData.C+";
  char treePath[] = "ak5/ProcessedTree";
  // This indicates the file to be opened
  Int_t isMC = 0;
  Int_t fileChoice = 1;
  gROOT->Reset();
  //gSystem->AddIncludePath("/home/hannu/Cern/jetscripts/data_analysis/analyzer");
  gROOT->ProcessLine(".L QCDModules/QCDEvent.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDMET.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDCaloJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDPFJet.cc+");
  gROOT->ProcessLine(".L QCDModules/QCDEventHdr.cc+");
  gROOT->ProcessLine(".L QCDModules/LorentzVector.h+");

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(scriptName);
  TChain *forest = new TChain(treePath);

  string tmpPath;

  /* 
  Indices:
  0   - pythia8 generated data
  1   - cmssw data (juska)
  2   - cmssw data (old)
  3-9 - real detector data 
  */
  string fileArray[] = {"fjpythia.root","ProcessedTree_RDMC_START53_V7N_try2.root",
    "27thJune_MC_NTuple.root",
    "29thAugust_Run2012A.root","29thAugust_Run2012B_JetHT.root",
    "29thAugust_Run2012B_JetMon.root","29thAugust_Run2012C_JetHT.root",
    "29thAugust_Run2012C_JetMon.root","29thAugust_Run2012D_JetHT.root",
    "29thAugust_Run2012D_JetMon.root"};
  Long64_t fileAmount = 9;

  /* 
  isDT is complementary to isMC, but isMC has many versions depending on
  the data file in use 
  */
  if (firstIdx == 1 ) { isMC = 1; } 
  else if ( firstIdx == 2) { isMC = 2; } 
  else if (firstIdx == 0) { isMC = 3; }

  if (firstIdx < 0) firstIdx = 0;
  if (lastIdx >= fileAmount) lastIdx = 8;
  if (firstIdx > lastIdx) firstIdx = lastIdx;

  cout << "First file read: " << fileArray[firstIdx] << endl;
  cout << "Last file read: " << fileArray[lastIdx] << endl;

  // For now, it is only possible to read some consequent files of the following list.
  for (Int_t fileChoice = firstIdx; fileChoice<=lastIdx; fileChoice++){
    tmpPath = dataPath;
    tmpPath += fileArray[fileChoice];
    /* 
    The AddFile routine adds only one tree at a time to a chain.
    Thus it is necessary to do some special handling if there are multiple
    trees in a single root file 
    */
    TFile *probe = new TFile(tmpPath.c_str());
    Int_t filesInDir = probe->GetDirectory("ak5")->GetNkeys();
    if (filesInDir > 3){
      TIter iter(probe->GetDirectory("ak5")->GetListOfKeys());
      TKey *tmpKey;
      while ( tmpKey = (TKey *) iter() ){
        if ( strcmp( tmpKey->GetName(), "ProcessedTree" ) !=  0) continue;
        string alterName(treePath);
        alterName += ";";
        std::stringstream numberConverter;
        numberConverter << tmpKey->GetCycle();
        alterName += numberConverter.str();
        forest->AddFile(tmpPath.c_str(),forest->kBigNumber,alterName.c_str());
      }
      delete tmpKey;
    } else {
      forest->AddFile(tmpPath.c_str());
    }
    delete probe;
  }
  AnalyzeData scriptObject(forest, loopLimit, isMC, testing);
  scriptObject.Loop(writeFile);
  delete forest;
}

