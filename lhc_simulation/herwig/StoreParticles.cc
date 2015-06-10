#include "StoreParticles.h"

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

/* A function that checks whether a photon is originated from a pi0 and that
 * the energy of the photon-pair corresponds to the pion. returns 0 if
 * the origin is not a pion with good energy and 1 if it is */
int StoreParticles::gammaChecker(const tPPtr& photon) {
    
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
int StoreParticles::isExcitedHadronState(const tPPtr& part, int quarkId) {
    assert( quarkId>=0 && quarkId<=6 );

    ParticleVector children = part->children();
    for (ParticleVector::const_iterator child = children.begin();
            child != children.end(); ++child ) 
        {
        if ( HadrFuncs::StatusCheck( quarkId, (*child)->id() ) ) return 1;
    }
    return 0;
}

void StoreParticles::particleAdd(const tPPtr& part, int saveStatus) 
{
    pEvent->AddPrtcl(part->momentum().x(),part->momentum().y(),part->momentum().z(),
                     part->momentum().t(),part->id(),part->data().charge(), saveStatus);
}

/* Implemented similarly as in cmssw. Not in active use, kept for reference. */
int StoreParticles::getStatusCode(const tPPtr& part) const
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
  
void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int status) 
{
    /* Rotate to CMS, extract final state particles and call analyze(particles). */
    AnalysisHandler::analyze(event, ieve, loop, status);
    
    if ( loop > 0 || status != 0 || !event ) return;
    
    pEvent->fWeight = event->weight();
    
    /* The hardest subprocess */
    eh = event->primaryCollision()->handler();
    const ParticleVector hardProc = event->primarySubProcess()->outgoing();
    vector<int> outgoingHard;
    for (std::size_t i = 0; i < hardProc.size(); ++i) {
        outgoingHard.push_back( hardProc[i]->number() );
        particleAdd(hardProc[i],3);
    }
    
    /* Final state particles */
    tPVector finals = event->getFinalState();
    for (tPVector::const_iterator pit = finals.begin(); pit != finals.end(); ++pit) {
        int absId = abs( (*pit)->id() );
        int saveStatus = 1;
        if ( (mode==0) && absId==ParticleID::gamma && gammaChecker(*pit) ) {
            saveStatus = 2;
        }
        
        particleAdd( *pit, saveStatus );
    }

//     tPVector parts;
//     event->select(std::back_inserter(parts),SelectAll());
//     /* Loop over all particles. */ 
//     for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit) {
//         int absId = abs( (*pit)->id() );
//         if ( absId==2101 || absId ==2203 || absId==82 ) continue; // uu, ud, p+rem
//         
//         int pStatus = getStatusCode( *pit );
//     }
    
    herwigTree->Fill();
    pEvent->Clear();
}

void StoreParticles::dofinish() 
{
    AnalysisHandler::dofinish();
    
    herwigTree->GetCurrentFile();
    herwigTree->AutoSave("Overwrite");
    herwigFile->Close();
    cout << "StoreParticles: a root tree has been written to a file" << endl;  
}

void StoreParticles::doinitrun() 
{
    AnalysisHandler::doinitrun();
    string fileName = "particles_herwig";
    mode = 0;
    cout << "halp" << endl;
    
    /* In a general multithread-case, generate a thread-unique root file name */
    string runInfo = generator()->runName();
    size_t pos = runInfo.find("++")+2;
    if (runInfo.size() > 0) {
        runInfo = runInfo.substr(pos);
        pos = runInfo.find("_");
        int multiplier = std::stoi(runInfo.substr(0,pos));
        
        runInfo = runInfo.substr(pos+1);
        pos = runInfo.find("_");
        int runIdx = std::stoi(runInfo.substr(0,pos));
        
        runInfo = runInfo.substr(pos+1);
        pos = runInfo.find("_");
        mode = std::stoi(runInfo.substr(0,pos));
        
        fileName += "_";
        if (mode == 0) {
            fileName += "generic_";
        } else if (mode == 1) {
            fileName += "dijet_";
        } else if (mode == 2) {
            fileName += "gammajet_";
        } else if (mode == 3) {
            fileName += "Zjet_";
        } else {
            cout << "Bad input file" << endl;
            fileName += "dijet_";
        }
        
        string fileNameFinal = fileName;
        fileName += std::to_string( generator()->N() );
        if (multiplier > 1) {
            if (runIdx==1) {
                fileNameFinal += std::to_string( multiplier*generator()->N() );
                fileNameFinal += ".root";
                
                TFile *outFile = new TFile(fileNameFinal.c_str(), "RECREATE");
                outFile->Close();
            }
            fileName += "_";
            fileName += std::to_string( runIdx );
        }
    }
    fileName += ".root";
    
    /* Setup a root tree in the selected file */
    herwigFile = new TFile (fileName.c_str(),"RECREATE");
    herwigFile->SetCompressionLevel(1); // by default file is compressed 

    if (!herwigFile) {
        cout << "StoreParticles: root file has not been created." << endl;
        return;
    }
    
    herwigTree = new TTree ("HerwigTree","Tree filled with herwig data.");
    if (!herwigTree) {
        cout << "StoreParticles: root tree has not been created." << endl;
        return;
    }
    herwigTree->SetAutoSave(100000000); /* autosave when 0.1 Gbyte written */
    herwigTree->SetCacheSize(10000000);  /* set a 10 MBytes cache (useless when writing local files) */

    TTree::SetBranchStyle(1); /* new style by default */
    pEvent = new PrtclEvent;
    TBranch *branch = herwigTree->Branch("event", &pEvent, 32000,4);
    branch->SetAutoDelete(kFALSE);
    herwigTree->BranchRef();
}

/* *** Attention *** This class-description must be correct. */
DescribeNoPIOClass<StoreParticles,AnalysisHandler>
  describejetanalysisStoreParticles("jetanalysis::StoreParticles", "../lib/libStoreParticles.so");

