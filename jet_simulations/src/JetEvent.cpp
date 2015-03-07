#include "JetEvent.h"

ClassImp(JetData)
ClassImp(JetEvent)

TClonesArray *JetEvent::fgJets = 0;


JetEvent::JetEvent(size_t tmpStore)
{
  Class()->IgnoreTObjectStreamer();
  SetNjet(0);
  if (!fgJets) fgJets = new TClonesArray("JetData",tmpStore);
  fJets = fgJets;
}

JetEvent::~JetEvent()
{
  Reset();
}

void JetEvent::Build( double Px, double Py, double Pz, double E, double Chf, 
  double Nhf, double Phf, double Elf, double Muf, char flav){
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  
  JetData *jet;
  
  jet = AddJet();
  jet->fPx = Px;
  jet->fPy = Py;
  jet->fPz = Pz;
  jet->fE = E;

  jet->fChf = Chf;
  jet->fNhf = Nhf;
  jet->fPhf = Phf;
  jet->fElf = Elf;
  jet->fMuf = Muf;
  
  jet->fFlavour = flav;
  
  TProcessID::SetObjectCount(ObjectNumber);
}

JetData* JetEvent::AddJet()
{
  if (fNjet>=10000) {
    fNjet = 0; // As a last resort, set the indexing to loop from beginning.
    cout << "The container size for SimParticle is too small, set a larger one in the constructor" << endl;
  }
  JetData *jet = (JetData*) fJets->ConstructedAt(fNjet++);
  return jet;
}

void JetEvent::Clear(Option_t* /* option */)
{
  fJets->Clear("C");
  SetNjet(0);
}

void JetEvent::Reset(Option_t* option)
{
  delete fgJets; fgJets = 0;
  SetNjet(0);
}

