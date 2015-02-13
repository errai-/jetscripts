#include "SimEvent.h"

ClassImp(SimParticle);
ClassImp(SimEvent);

TClonesArray *SimEvent::fgParts = 0;

Double_t SimParticle::Pt() const
{
  return sqrt(pow(fPx,2)+pow(fPy,2));
}

Double_t SimParticle::Eta() const
{
  Double_t P = sqrt(pow(fPx,2) + pow(fPy,2) + pow(fPz,2));
  return 0.5*log( (P+fPz)/(P-fPz));
}

Double_t SimParticle::Phi() const
{
  Double_t rawAngle = TMath::ATan2( fPy, fPx );
  if ( fPx >= 0 ) {
    return rawAngle;
  } else if ( fPy >= 0 ) {
    return TMath::Pi() - rawAngle;
  } else {
    return -TMath::Pi() - rawAngle;
  }
}

SimEvent::SimEvent(size_t tmpStore)
{
  SetNpart(0);
  if (!fgParts) fgParts = new TClonesArray("SimParticle",tmpStore);
  fParts = fgParts;
}

SimEvent::~SimEvent()
{
  Reset();
}

void SimEvent::Build( double Px, double Py, double Pz, double E, int Id){
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  
  SimParticle *part;
  
  part = AddParticle();
  part->fPx = Px;
  part->fPy = Py;
  part->fPz = Pz;
  part->fE  = E;
  part->fPDGCode = Id;
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

