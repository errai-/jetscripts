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
//     fChm = vars.chm;
//     fNhm = vars.nhm;
//     fPhm = vars.phm;
//     fElm = vars.elm;
//     fMum = vars.mum;

    fFlav = flav;                       /* Jet flavor, depends on the definition in use */
    fConstituents = vars.constituents;  /* Number of constituents, exact definition varies (e.g. pT threshold) */
    fPTD = vars.PTD;                    /* pT partition function */
    fSigma2 = vars.Sigma2;              /* The inner jet radius */

    fDR = vars.DR;                      /* The distance between the jet axis and the parton */
    fAlpha = vars.Alpha;                /* The alpha parameter, depends on the event type */
    fDPhi = vars.DPhi;                  /* The phi parameter, depends on the event type */

    fPartonPT = vars.partonPT;          /* pT of the parton matched with the jet */
    fMatchPT = vars.matchPT;            /* pT of the opposing jet or such */
    //fDR_Next = vars.nextDR;             /* distance to the next parton to match */
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
