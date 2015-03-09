#include "PrtclEvent.h"

ClassImp(PrtclData)
ClassImp(PrtclEvent)

TClonesArray *PrtclEvent::fgParts = 0;

PrtclEvent::PrtclEvent(size_t tmpStore)
{
  Class()->IgnoreTObjectStreamer();
  SetNpart(0);
  if (!fgParts) fgParts = new TClonesArray("PrtclData",tmpStore);
  fSizeLim = tmpStore;
  fParts = fgParts;
}

void PrtclEvent::Build( double Px, double Py, double Pz, double E, int Id, double Charge, int pi0Gamma, int jetFlavor, int isExcitedState){
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  
  PrtclData *part;
  
  part = AddParticle();
  part->fP4.SetPxPyPzE(Px,Py,Pz,E);

  part->fPDGCode = Id;
  part->fChargeTimes3 = TMath::Nint(3*Charge);
  part->IsPi0Photon = (pi0Gamma ? true : false);
  part->IsJetFlavor = (jetFlavor ? true : false);
  part->IsExcitedState = (isExcitedState ? true : false);
  TProcessID::SetObjectCount(ObjectNumber);
}

PrtclData* PrtclEvent::AddParticle()
{
  if (fNpart>fSizeLim) {
    fNpart = 0; // As a last resort, set the indexing to loop from beginning.
    cout << "The container size for PrtclData is too small, set a larger one in the constructor" << endl;
  }
  PrtclData *part = (PrtclData*) fParts->ConstructedAt(fNpart++);
  return part;
}

void PrtclEvent::Clear(Option_t* /* option */)
{
  fParts->Clear("C");
  SetNpart(0);
}

void PrtclEvent::Reset(Option_t* option)
{
  delete fgParts; fgParts = 0;
  SetNpart(0);
}
