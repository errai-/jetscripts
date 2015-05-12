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


void JetData::SetParams(double chf, double nhf, double phf, double elf, double 
   muf, double chm, double nhm, double phm, double elm, double mum, double weight, int flav)
{
    fChf = chf;
    fNhf = nhf;
    fPhf = phf;
    fElf = elf;
    fMuf = muf;
    fChm = chm;
    fNhm = nhm;
    fPhm = phm;
    fElm = elm;
    fMum = mum;
    fFlav = flav;
}


JetEvent::JetEvent(size_t tmpStore)
{
    Class()->IgnoreTObjectStreamer();
    fN_Jet = 0;
    if (!fgJets) fgJets = new TClonesArray("JetData",tmpStore);
    fSizeLim = tmpStore;
    fJets = fgJets;
}


JetEvent::~JetEvent()
{
    Reset();
}


void JetEvent::AddJet( double Px, double Py, double Pz, double E, double Chf, 
    double Nhf, double Phf, double Elf, double Muf, double Chm, double Nhm, 
    double Phm, double Elm, double Mum, double weight, int flav)
{
    Int_t ObjectNumber = TProcessID::GetObjectCount();
    
    JetData *jet;
    
    jet = InitJet();
    jet->SetPxPyPzE(Px,Py,Pz,E);
    jet->SetParams(Chf,Nhf,Phf,Elf,Muf,Chm,Nhm,Phm,Elm,Mum,weight,flav);
    
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