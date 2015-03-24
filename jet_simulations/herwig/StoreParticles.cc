#include "StoreParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


using namespace jetanalysis;
using namespace ThePEG;

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

#include <iostream>

using std::cout;
using std::endl;

void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int status) 
{
   /* Rotate to CMS, extract final state particles and call analyze(particles). */
   AnalysisHandler::analyze(event, ieve, loop, status);
   if ( loop > 0 || !event ) return;
   cout << status << endl;
   // if ( loop > 0 || status != 0 || !event ) return;

   //tPVector parts=event->getFinalState();
   tPVector parts;
   event->select(std::back_inserter(parts),SelectAll());
   
   /* Loop over all particles. */ 
   for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit){
      pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
         (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(),1);
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

  // create ROOT File
  herwigFile = new TFile ("herwig_particles.root","RECREATE");
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
  herwigTree->SetAutoSave(1000000000); /* autosave when 1 Gbyte written */
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

