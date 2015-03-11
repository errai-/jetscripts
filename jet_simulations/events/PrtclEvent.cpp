#include "events/PrtclEvent.h"

ClassImp(PrtclData)
ClassImp(PrtclEvent)


TClonesArray *PrtclEvent::fgPrtcls = 0;

using std::cout;
using std::endl;


void PrtclData::SetPxPyPzE(double px, double py, double pz, double e)
{
   fP4.SetPxPyPzE(px,py,pz,e);
}


void PrtclData::SetParams(int id, double charge, int pi0Gamma, int jetFlavor, 
int    isExcitedState )
{
   fPDGCode = id;
   fChargeTimes3 = TMath::Nint(3*charge);
   fIsPi0Photon = (pi0Gamma ? true : false);
   fIsJetFlavor = (jetFlavor ? true : false);
   fIsExcitedState = (isExcitedState ? true : false);
}


PrtclEvent::PrtclEvent(size_t tmpStore)
{
   Class()->IgnoreTObjectStreamer();
   fN_Prtcl = 0;
   if (!fgPrtcls) fgPrtcls = new TClonesArray("PrtclData",tmpStore);
   fSizeLim = tmpStore;
   fPrtcls = fgPrtcls;
}


void PrtclEvent::AddPrtcl( double px, double py, double pz, double e, int id, 
   double charge, int pi0Gamma, int jetFlavor, int isExcitedState)
{
   int ObjectNumber = TProcessID::GetObjectCount();
  
   PrtclData *part;
   part = InitPrtcl();
   part->SetPxPyPzE(px,py,pz,e);
   part->SetParams(id,charge,pi0Gamma,jetFlavor,isExcitedState);

   TProcessID::SetObjectCount(ObjectNumber);
}


PrtclData* PrtclEvent::InitPrtcl()
{
   if (fN_Prtcl>fSizeLim) {
      // As a last resort, set the indexing to loop from beginning.
      fN_Prtcl = 0; 
      cout << "The container size for PrtclData is too small, set a larger one" 
         "in the constructor" << endl;
   }
   PrtclData *part = (PrtclData*) fPrtcls->ConstructedAt(fN_Prtcl++);
   return part;
}


void PrtclEvent::Clear(Option_t* /* option */)
{
   fPrtcls->Clear("C");
   fN_Prtcl = 0;
}


void PrtclEvent::Reset(Option_t* option)
{
   delete fgPrtcls; fgPrtcls = 0;
   fN_Prtcl = 0;
}
