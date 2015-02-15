#include "SimEvent.h"

ClassImp(SimParticle)
ClassImp(SimEvent)

TClonesArray *SimEvent::fgParts = 0;


SimEvent::SimEvent(size_t tmpStore)
{
  Class()->IgnoreTObjectStreamer();
  SetNpart(0);
  if (!fgParts) fgParts = new TClonesArray("SimParticle",tmpStore);
  fParts = fgParts;
}

SimEvent::~SimEvent()
{
  Reset();
}

void SimEvent::Build( double Px, double Py, double Pz, double E, int Id, double Charge, int pi0Gamma, int jetFlavor, int isExcitedState){
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  
  SimParticle *part;
  
  part = AddParticle();
  part->fPx = Px;
  part->fPy = Py;
  part->fPz = Pz;
  part->fE = E;

  part->fPDGCode = Id;
  part->fChargeTimes3 = TMath::Nint(3*Charge);
  part->IsPi0Photon = pi0Gamma ? true : false;
  part->IsJetFlavor = jetFlavor ? true : false;
  part->IsExcitedState = isExcitedState ? true : false;
  TProcessID::SetObjectCount(ObjectNumber);
}

SimParticle* SimEvent::AddParticle()
{
  if (fNpart>=10000) {
    fNpart = 0; // As a last resort, set the indexing to loop from beginning.
    cout << "The container size for SimParticle is too small, set a larger one in the constructor" << endl;
  }
  SimParticle *part = (SimParticle*) fParts->ConstructedAt(fNpart++);
  return part;
}

void SimEvent::Clear(Option_t* /* option */)
{
  fParts->Clear("C");
  SetNpart(0);
}

void SimEvent::Reset(Option_t* option)
{
  delete fgParts; fgParts = 0;
  SetNpart(0);
}

