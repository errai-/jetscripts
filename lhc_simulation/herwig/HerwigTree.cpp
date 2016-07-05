#include "HerwigTree.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


/* STRUCTURES FOR STORING NECESSARY PARTICLE DATA ACCORDING TO THE EVENT TYPE: */

/* The danish number system is very complicated, so better just divide with GeV */
void HerwigTree::particleAdd(const tPPtr& part, int saveStatus)
{
    mPrtclEvent->AddPrtcl(part->momentum().x()/GeV,
                          part->momentum().y()/GeV,
                          part->momentum().z()/GeV,
                          part->momentum().t()/GeV,
                          part->id(),
                          saveStatus);
    cout << part->number() << " " << saveStatus << endl;
    return;
}

void HerwigTree::print_parents(const tPPtr& part) {
    cout << " " << part->number() << " " << part->id() << "  ";
    const tParticleVector parents = part->parents();
    for (unsigned i = 0u; i < parents.size(); ++i)
        cout << parents[i]->number() << "-" << getStatusCode(parents[i]) << " ";
    cout << endl;
    for (unsigned i = 0u; i < parents.size(); ++i)
        print_parents(parents[i]);
}

void smMomenta(Lorentz5Momentum & psum, tPPtr parent) {
  if(!parent->children().empty()) {
    for(unsigned int ix=0;ix<parent->children().size();++ix)
      smMomenta(psum,parent->children()[ix]);
  }
  else {
    cout << parent->id() << endl;
    psum += parent->momentum();
  }
}


// Fancy plots: event->printGraphviz()
// All particles:
//  tPVector all;
//  event->select(std::back_inserter(all),SelectAll());
//  all.begin(), all.end() etc.
// Hard process (method 1 - brings also the id 82 collimations)
//  const ParticleVector hordProk = event->primarySubProcess()->outgoing(); (const_iterator)
// Hard process (method 2 - preferred)
//  tPVector hardProc = event->primaryCollision()->step(0)->getFinalState(); (const_iterator)
void HerwigTree::analyze(tEventPtr event, long ieve, int loop, int status) {

    if (ieve%mTimerStep==0&&ieve>0) mTimer.printTime();
    cout << ieve << endl;
    if ( loop > 0 || status != 0 || !event ) return;

    try {
        mPrtclEvent->Clear();
        mSpecialIndices.clear();
        int hardProcCount = 0;
        event->printGraphviz();

        mPrtclEvent->fWeight = event->weight();
        mHard = event->primaryCollision()->handler();

        /* The hardest subprocess */
        tPVector hardProc = event->primaryCollision()->step(0)->getFinalState();
        int leptons = 0;
        for (tPVector::const_iterator part = hardProc.begin(); part != hardProc.end(); ++part) {
            int absId = abs((*part)->id());
            bool gammaCase = (mMode==2 && absId==ParticleID::gamma );
            bool ZCase = (mMode==3 && absId==ParticleID::muminus );

            if (gammaCase) {
                gammaAdd(*part);
                ++hardProcCount;
            } else if (ZCase) {
                muonAdd(*part);
                ++hardProcCount;
            } else if ( absId < 10 || absId==21 ) {
                particleAdd(*part,3);
            } else {
                if (absId > 10 && absId < 20 ) {
                    leptonAdd(*part);
                    ++leptons;
                }
                continue;
            }
            ++hardProcCount;
        }

        /* ttbar events: seek lepton+jets (one w to leptons, one to quarks) */
        if ( mMode==4 && mSpecialIndices.size()!=2 ) {
            if (!generator()->repeatEvent()) throw runtime_error("Repeating events does not work.");
            return;
        }

        if (mMode==3) --hardProcCount;

        /* Sanity checks */
        if (   (mMode<4 && hardProcCount !=2)
            || (mMode==2 && mSpecialIndices.size()!=1)
            || (mMode==3 && mSpecialIndices.size()!=2)
            || (mMode==4 && hardProcCount != 4) )
        {
            throw std::logic_error("Unexpected hard process structure");
        }

        /* Final state particles */
        tPVector finals = event->getFinalState();
        for (tPVector::const_iterator part = finals.begin(); part != finals.end(); ++part) {
            int finalIdx = (*part)->number();

            if ( std::count( mSpecialIndices.begin(), mSpecialIndices.end(), finalIdx)>0 ) continue;

            int absId = abs( (*part)->id() );

            /* pi0 photons in a generic event have the status 2 */
            int saveStatus = 1;
            if ( (mMode==0) && absId==ParticleID::gamma && gammaChecker(*part) ) saveStatus = 2;
            particleAdd( *part, saveStatus );
        }

        mTree->Fill();

    } catch (std::exception& e) {
        cout << "An error occurred: " << e.what() << endl;
    }
}


bool HerwigTree::gammaAdd(tPPtr gamma)
{
    while (gamma->decayed()) {
        const ParticleVector children = gamma->children();
        /* No pair production */
        if (children.size()!=1) {
            if (!generator()->repeatEvent()) throw runtime_error("Repeating events does not work.");
            return false;
        }
        gamma = children[0];
    }
    mSpecialIndices.push_back( gamma->number() );
    particleAdd(gamma,2);
    return true;
}


bool HerwigTree::muonAdd(tPPtr muon)
{
    while (muon->decayed() > 0) {
        const ParticleVector children = muon->children();
        for (ParticleVector::const_iterator child = children.begin(); child != children.end(); ++child) {
            if ( abs((*child)->id())==ParticleID::muminus ) { muon = *child; break; }
        }
    }
    mSpecialIndices.push_back( muon->number() );
    particleAdd(muon,2);
    return true;
}


bool HerwigTree::leptonAdd(tPPtr lepton)
{
    int type = abs( lepton->id() )%2;

    while (lepton->decayed()) {
        const ParticleVector children = lepton->children();
        bool stuck = true;
        for (ParticleVector::const_iterator child = children.begin(); child != children.end(); ++child) {
            int absId = abs((*child)->id());
            if ( absId<20 && absId>10 && absId%2==type ) {
                lepton = *child; stuck = false; break;
            }
        }
        /* Check if stuck in a loop (for instance if lepton goes to hadrons) */
        if ( stuck ) return false;
    }
    mSpecialIndices.push_back( lepton->number() );
    particleAdd(lepton,2);
    return true;
}


/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
int HerwigTree::gammaChecker(const tPPtr& photon) {

    /* One mother, which is pi0 */
    const tParticleVector parents = photon->parents();
    if ( parents.size()!=1 || parents[0]->id()!=ParticleID::pi0 ) return 0;

    const ParticleVector children = parents[0]->children();
    if ( children.size()!=2 ) return 0;

    double eDifference = fabs( parents[0]->momentum().t() -
    children[0]->momentum().t() - children[1]->momentum().t() );
    if ( eDifference > 0.001 ) return 0;

    return 1;
}


/* See HadronAndPartonSelector.cc in cmssw, indicates whether a hadron (used for
 * flavour inspection) is in an excited state or not. Checks whether a hadron has
 * a daughter of the same flavour. The parameter quarkId is a PDG quark flavour. */
int HerwigTree::isExcitedHadronState(const tPPtr& part, int quarkId) {
    assert( quarkId>=0 && quarkId<=6 );

    ParticleVector children = part->children();
    for (ParticleVector::const_iterator child = children.begin();
         child != children.end(); ++child )
         {
             if ( HadrFuncs::StatusCheck( quarkId, (*child)->id() ) ) return 1;
         }
         return 0;
}


/* NOT IN CURRENT PRODUCTION, KEPT FOR REFERENCE: */


/* Implemented similarly as in cmssw. */
int HerwigTree::getStatusCode(const tPPtr& part) const
{
    int status = 1;
    if ( !part->children().empty() || part->next() ) {
        tStepPtr step = part->birthStep();
        if ((!step || (step && (!step->handler() || step->handler() == mHard)))
            && part->id() != 82)
        {
            status = 3;
        } else {
            status = 2;
        }
    }
    return status;
}
