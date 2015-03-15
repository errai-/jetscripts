/* -*- C++ -*-
 *
 * This is the implementation of the non-inlined, non-templated member
 * functions of the StoreParticles class. */

#include "herwig/StoreParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"



using namespace jetanalysis;

StoreParticles::StoreParticles() {}

StoreParticles::~StoreParticles() {}


#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int status) {
   AnalysisHandler::analyze(event, ieve, loop, status);
   
   /* Rotate to CMS, extract final state particles and call analyze(particles). */
   if ( loop > 0 || status != 0 || !event ) return;

   /* Get the final-state particles */
   tPVector parts=event->getFinalState();
  
   /* Loop over all particles. 
    * This should be cleaned up. Herwig is still a bit of a work in progress,
    * so some of these commented lines could be useful.
    * pEvent is used only partially, since the "more advanced" parts are at
    * the moment working only with pythia8 */
   for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit){
      pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
         (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(),1);
   }
   // TODO ghost particles/partons
   
   herwigTree->Fill();
   pEvent->Clear();
}

LorentzRotation StoreParticles::transform(tcEventPtr event) const {
   return LorentzRotation();
   /* Return the Rotation to the frame in which you want to perform the analysis. */
}

void StoreParticles::analyze(const tPVector & parts, double weight) {
   AnalysisHandler::analyze(parts);
   /* Calls analyze() for each particle. */
}

void StoreParticles::analyze(tPPtr, double weight) {}


IBPtr StoreParticles::clone() const {
   return new_ptr(*this);
}

IBPtr StoreParticles::fullclone() const {
   return new_ptr(*this);
}


/* If needed, insert default implementations of virtual function defined
 * in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs). */

void StoreParticles::doupdate() {
   AnalysisHandler::doupdate();
   // First update base class.
   bool redo = touched();
   // redo if touched.
   //  UpdateChecker::check(aDependentMember, redo);
   // Update referenced objects on which this depends redo is set to true
   // if the dependent object is touched.
   //  for_each(ContainerOfDependencies, UpdateChecker(redo));
   // Update a container of references.
   //  for_each(MapOfDependencies, UpdateMapChecker(redo));
   // Update a map of references.
   if ( !redo ) return;
   // return if nothing has been touched. Otherwise do the actual update.
   //  touch()
  // Touch if anything has changed.
}

void StoreParticles::doinit() {
  AnalysisHandler::doinit();
}

void StoreParticles::dofinish() {
  AnalysisHandler::dofinish();
  
  herwigTree->GetCurrentFile();
  herwigTree->AutoSave("Overwrite");
  herwigFile->Close();
  cout << "StoreParticles: root tree has been written to a file" << endl;  
}

void StoreParticles::doinitrun() {
  AnalysisHandler::doinitrun();

  // create ROOT File
  herwigFile = new TFile ("herwig_particles.root","RECREATE");
  Int_t comp   = 1;       // by default file is compressed
  herwigFile->SetCompressionLevel(comp);

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

void StoreParticles::rebind(const TranslationMap & trans) {
  AnalysisHandler::rebind(trans);
}

IVector StoreParticles::getReferences() {
  IVector ret = AnalysisHandler::getReferences();
  return ret;
}



// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<StoreParticles,AnalysisHandler>
  describejetanalysisStoreParticles("jetanalysis::StoreParticles", "StoreParticles.so");

void StoreParticles::Init() {

  static ClassDocumentation<StoreParticles> documentation
    ("There is no documentation for the StoreParticles class");

}

