#include "PrtclEvent.h"

ClassImp(PrtclData)
ClassImp(PrtclEvent)


TClonesArray *PrtclEvent::fgPrtcls = 0;

using std::cout;
using std::endl;


void PrtclData::SetPxPyPzE(double px, double py, double pz, double e)
{
    fP4.SetPxPyPzE(px,py,pz,e);
}


void PrtclData::SetParams(int id, int status)
{
    fPDGCode = id;
    fAnalysisStatus = status;
}


PrtclEvent::PrtclEvent()
{
    Class()->IgnoreTObjectStreamer();
    PrtclData::Class()->IgnoreTObjectStreamer();
    fN_Prtcl = 0;
    /* TClonesArray is dynamic in size, use 1000 for an approximation of the memory needed */
    if (!fgPrtcls) fgPrtcls = new TClonesArray("PrtclData",5000);
    fPrtcls = fgPrtcls;
}


void PrtclEvent::AddPrtcl(double px, double py, double pz, double e,
                          int id, int status)
{
    int ObjectNumber = TProcessID::GetObjectCount();

    PrtclData *part;
    part = InitPrtcl();
    part->SetPxPyPzE(px,py,pz,e);
    part->SetParams(id,status);
    part->SetBit(kCanDelete);
    part->SetBit(kMustCleanup);

    TProcessID::SetObjectCount(ObjectNumber);
}


PrtclData* PrtclEvent::InitPrtcl()
{
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
