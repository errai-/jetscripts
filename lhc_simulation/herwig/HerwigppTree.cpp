#include "HerwigppTree.h"

using namespace jetanalysis;
using namespace ThePEG;

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

#include <iostream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::string;

/* STRUCTURES FOR STORING NECESSARY PARTICLE DATA ACCORDING TO THE EVENT TYPE: */

/* The danish number system is very complicated, so better just divide with GeV */
void HerwigppTree::particleAdd(const tPPtr& part, int saveStatus) 
{
    pEvent->AddPrtcl(part->momentum().x()/GeV,part->momentum().y()/GeV,
                     part->momentum().z()/GeV,part->momentum().t()/GeV,
                     part->id(), saveStatus);
}

void HerwigppTree::analyze(tEventPtr event, long ieve, int loop, int status)
{
    if (ieve%timerStep==0&&ieve>0) timer->printTime();

    if ( loop > 0 || status != 0 || !event ) return;

    pEvent->fWeight = event->weight();
    eh = event->primaryCollision()->handler();
    //event->printGraphviz();

    vector<std::size_t> specialIndices;
    int hardProcCount = 0;

    /* The hardest subprocess */
    if (mode==4) {
        tPVector out = event->primaryCollision()->step(0)->getFinalState();
        int leptons = 0;
        for (tPVector::const_iterator part = out.begin(); part != out.end(); ++part) {
            int absId = abs( (*part)->id() );
            
            if ( absId < 10 ) {
                particleAdd(*part,3);
                ++hardProcCount;
            } else if ( absId > 10 && absId < 20 ) {
                // TODO: chase the final leptons
                particleAdd(*part,2);
                ++leptons;
            }
        }
        if (leptons==2) {
            pEvent->Clear();
            return;
        }
    } else {
        const ParticleVector hardProc = event->primarySubProcess()->outgoing();
        
        for (ParticleVector::const_iterator part = hardProc.begin(); part != hardProc.end(); ++part) {
            bool gammaCase = (mode==2 && abs((*part)->id())==ParticleID::gamma );
            bool ZCase = (mode==3 && abs((*part)->id())==ParticleID::muminus );
            
            if (gammaCase) {
                PPtr gamma = *part;
                while (gamma->decayed()) {
                    const ParticleVector children = gamma->children();
                    /* No pair production */
                    if (children.size()!=1) {
                        pEvent->Clear();
                        return;
                    }
                    gamma = children[0];
                }
                specialIndices.push_back( gamma->number() );
                particleAdd(gamma,2);
            } else if (ZCase) {
                tPPtr muon = (*part);
                while (muon->decayed() > 0) {
                    const ParticleVector children = muon->children();
                    for (ParticleVector::const_iterator child = children.begin(); child != children.end(); ++child) {
                        if ( abs((*child)->id())==ParticleID::muminus ) { muon = *child; break; }
                    }
                }
                specialIndices.push_back( muon->number() );
                particleAdd(muon,2);
            } else {
                particleAdd(*part,3);
            }
        }
        
        if (  (mode==2&&specialIndices.size()!=1) 
           || (mode==3&&specialIndices.size()!=2)
           || (hardProc.size() != 2 && hardProc.size() != 3) )
        {
            cout << "Unexpected hard process structure" << endl;
        }
    }
    
    /* Final state particles */
    tPVector finals = event->getFinalState();
    for (tPVector::const_iterator part = finals.begin(); part != finals.end(); ++part) {
        int finalIdx = (*part)->number();
        
        if ( std::count( specialIndices.begin(), specialIndices.end(), finalIdx)>0 ) continue;
        
        int absId = abs( (*part)->id() );
        
        /* pi0 photons in a generic event have the status 2 */
        int saveStatus = 1;
        if ( (mode==0) && absId==ParticleID::gamma && gammaChecker(*part) ) saveStatus = 2;
        particleAdd( *part, saveStatus );
    }
    
    tree->Fill();
    pEvent->Clear();
}


/* NOT IN CURRENT PRODUCTION, KEPT FOR REFERENCE: */

/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
int HerwigppTree::gammaChecker(const tPPtr& photon) {
    
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
int HerwigppTree::isExcitedHadronState(const tPPtr& part, int quarkId) {
    assert( quarkId>=0 && quarkId<=6 );

    ParticleVector children = part->children();
    for (ParticleVector::const_iterator child = children.begin();
         child != children.end(); ++child ) 
    {
        if ( HadrFuncs::StatusCheck( quarkId, (*child)->id() ) ) return 1;
    }
    return 0;
}


/* Implemented similarly as in cmssw. */
int HerwigppTree::getStatusCode(const tPPtr& part) const
{
    int status = 1;
    if ( !part->children().empty() || part->next() ) {
        tStepPtr step = part->birthStep();
        if ((!step || (step && (!step->handler() || step->handler() == eh))) 
            && part->id() != 82)
        {
            status = 3;
        } else {
            status = 2;
        }
    }
    return status;
}

