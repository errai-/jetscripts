#include "JetEvent.h"

ClassImp(JetData)
ClassImp(JetEvent)

TClonesArray *JetEvent::fgJets = 0;

using std::cout;
using std::endl;


void JetData::SetPxPyPzE(double px, double py, double pz, double e)
{
    fP4.SetPxPyPzE(px,py,pz,e);
}


void JetData::SetParams(JetVariables& vars, int flav)
{
    fChf = vars.chf;
    fNhf = vars.nhf;
    fPhf = vars.phf;
    fElf = vars.elf;
    fMuf = vars.muf;
    fChm = vars.chm;
    fNhm = vars.nhm;
    fPhm = vars.phm;
    fElm = vars.elm;
    fMum = vars.mum;
    fFlav = flav;
    fConstituents = vars.constituents;
    fPTD = vars.PTD;
    fSigma2 = vars.Sigma2;
    fDR = vars.DR;
    fAlpha = vars.Alpha;
    fDPhi = vars.DPhi;
    
    fPartonPT = vars.partonPT;
    fMatchPT = vars.matchPT;
    fMaxDR = vars.maxDR;
}


JetEvent::JetEvent(size_t tmpStore)
{
    Class()->IgnoreTObjectStreamer();
    JetData::Class()->IgnoreTObjectStreamer();
    fN_Jet = 0;
    if (!fgJets) fgJets = new TClonesArray("JetData",tmpStore);
    fSizeLim = tmpStore;
    fJets = fgJets;
}


JetEvent::~JetEvent()
{
    Reset();
}


void JetEvent::AddJet( double Px, double Py, double Pz, double E, JetVariables& jetVars, double weight, int flav)
{
    Int_t ObjectNumber = TProcessID::GetObjectCount();
    
    JetData *jet;
    
    jet = InitJet();
    fWeight = weight;
    jet->SetPxPyPzE(Px,Py,Pz,E);
    jet->SetParams(jetVars,flav);
    
    TProcessID::SetObjectCount(ObjectNumber);
}


JetData* JetEvent::InitJet()
{
    assert(fSizeLim>fN_Jet);
    JetData *jet = (JetData*) fJets->ConstructedAt(fN_Jet++);
    return jet;
}


void JetEvent::Clear(Option_t* /* option */)
{
    fJets->Clear("C");
    fN_Jet = 0;
}


void JetEvent::Reset(Option_t* option)
{
    delete fgJets; fgJets = 0;
    fN_Jet = 0;
}