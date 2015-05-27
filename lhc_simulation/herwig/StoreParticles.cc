#include "StoreParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace jetanalysis;
using namespace ThePEG;

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

#include <iostream>
#include <fstream>
#include <string>

#include "../generic/help_functions.h"
#include "herwig_functions.h"

using std::cout;
using std::endl;
using std::string;

int StoreParticles::getStatusCode(tPPtr part) const
{
   int status = 1;
   size_t nChildren = part->children().size();
   if ( nChildren > 0 || part->next() ) {
      tStepPtr step = part->birthStep();
      if ((!step || (step && (!step->handler() || step->handler() == eh))) && part->id() != 82)
         status = 3;
      else
         status = 2;
   }
   /* This is the official ThePEG method for obtaining status, but it does not
    * work with the CMSSW interpretation of jet flavours. */
   // if ( nChildren > 0 || part->next() ) status = 11; 
   // if ( nChildren > 1 ) { 
   //    long id = part->data().id();
   //    if ( BaryonMatcher::Check(id) || MesonMatcher::Check(id) ||
   //       id == ParticleID::muminus || id == ParticleID::muplus ||
   //       id == ParticleID::tauminus || id == ParticleID::tauplus )
   //          if ( part->mass() <= part->data().massMax() &&
   //             part->mass() >= part->data().massMin() ) status = 2;
   // }
   return status;
}
  
void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int status) 
{
   /* Rotate to CMS, extract final state particles and call analyze(particles). */
   AnalysisHandler::analyze(event, ieve, loop, status);
   if ( loop > 0 || !event ) return;
   // if ( loop > 0 || status != 0 || !event ) return;

   //tPVector parts=event->getFinalState();
   tPVector parts;
   event->select(std::back_inserter(parts),SelectAll());
   eh = event->primaryCollision()->handler();
   
   /* Loop over all particles. */ 
   for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit) {
      int absId = abs( (*pit)->id() );
      if ( absId==2101 || absId ==2203 || absId==82 ) continue; // uu, ud, p+rem
      int pStatus = getStatusCode( *pit );
      if ( pStatus==3 ) continue; // Beam particles and partons, not of interest
      
      /* Normal end state particles */
      if (pStatus == 1) {
         int pi0Gamma = 0;
         if ( absId==22 && gammaChecker(*pit) ) pi0Gamma = 1;
         pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
            (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(), pi0Gamma ? 10 : 1);
      }
      
      /* Ghost particles */
      int ghostStatus = 0;
      
      if ( (pStatus==2) && ( (absId<=6 && absId>0) || absId==21 ) ) {
            ghostStatus = 11; /* Partons */
      } else if (absId >= 100) { /* Status codes below this are not conventional hadrons */
         if (HadrFuncs::hasBottom(absId) && !isExcitedHadronState(*pit,5) ) {
            ghostStatus = 12; /* b Hadrons */
         }
         if (HadrFuncs::hasCharm(absId) && !isExcitedHadronState(*pit,4) ) {
            ghostStatus = 13; /* c Hadrons */
         }
      } /* Add Leptons here if needed */
      
      if (ghostStatus) {
         pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
            (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(),ghostStatus);
      }
   }
   // TODO ghost particles/partons
   
   herwigTree->Fill();
   pEvent->Clear();
}

void StoreParticles::dofinish() 
{
  AnalysisHandler::dofinish();
  
  herwigTree->GetCurrentFile();
  herwigTree->AutoSave("Overwrite");
  herwigFile->Close();
  cout << "StoreParticles: root tree has been written to a file" << endl;  
}

void StoreParticles::doinitrun() 
{
    AnalysisHandler::doinitrun();
    string fileName = "particles_herwig";
    
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
        int mode = std::stoi(runInfo.substr(0,pos));
        
        fileName += "_";
        if (mode == 0) {
            fileName += "generic_";
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
    
    herwigFile = new TFile (fileName.c_str(),"RECREATE");
    herwigFile->SetCompressionLevel(1); // by default file is compressed 

    if (!herwigFile) {
        cout << "StoreParticles: root file has not been created..." << endl;
        return;
    }
    
    // create ROOT Tree
    herwigTree = new TTree ("HerwigTree","Tree filled with herwig data.");
    if (!herwigTree) {
        cout << "StoreParticles: root tree has not been created..." << endl;
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

/* *** Attention *** The following static variable is needed for the type
 * description system in ThePEG. Please check that the template arguments
 * are correct (the class and its base class), and that the constructor
 * arguments are correct (the class name and the name of the dynamically
 * loadable library where the class implementation can be found). */
DescribeNoPIOClass<StoreParticles,AnalysisHandler>
  describejetanalysisStoreParticles("jetanalysis::StoreParticles", "../lib/libStoreParticles.so");

