// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoreParticles class.
//

#include "StoreParticles.h"
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

void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  if ( loop > 0 || state != 0 || !event ) return;

  /** get the final-state particles */
  tPVector parts=event->getFinalState();
  
  particles = parts.size();
  
  int counter = 0; // index of the current particle
  /** loop over all particles */
  for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit){
    //if( ChargedSelector::Check(**pit) )
    px[counter] = (*pit)->momentum().x();
    py[counter] = (*pit)->momentum().y();
    pz[counter] = (*pit)->momentum().z();
    e[counter] = (*pit)->momentum().t();
    id[counter] = (*pit)->id();
    counter++;
  }
  // Fill TTree record
  herwigTree->Fill();
}

LorentzRotation StoreParticles::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void StoreParticles::analyze(const tPVector & parts, double weight) {
  AnalysisHandler::analyze(parts);
  // Calls analyze() for each particle.
}

void StoreParticles::analyze(tPPtr, double weight) {}


IBPtr StoreParticles::clone() const {
  return new_ptr(*this);
}

IBPtr StoreParticles::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

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
  herwigTree->Write();
  herwigFile->Close();
  cout << "StoreParticles: root tree has been written to a file" << endl;  
}

void StoreParticles::doinitrun() {
  AnalysisHandler::doinitrun();

  // create ROOT Tree
  herwigTree = new TTree ("h76","myAnalysis root tree", 1);
  if (!herwigTree) {
    cout << "StoreParticles: root tree has not been created..." << endl;
    return;
  }
  // create ROOT File
  herwigFile = new TFile ("herwigData.root","RECREATE");
  if (!herwigFile) {
    cout << "StoreParticles: root file has not been created..." << endl;
    return;
  }
  herwigTree->SetDirectory (herwigFile);

  // define ROOT Tree branches/leaves  
  herwigTree->Branch ("particles",  &particles, "particles/I");
  herwigTree->Branch ("px",      px,      "px[particles]/D");
  herwigTree->Branch ("py",      py,      "py[particles]/D");
  herwigTree->Branch ("pz",      pz,      "pz[particles]/D");
  herwigTree->Branch ("e",       e,        "e[particles]/D");
  herwigTree->Branch ("id",      id,      "id[particles]/I");
//   herwigTree->Branch ("Nqurk",   &Nqurk,  "Nqurk/I");
//   herwigTree->Branch ("Nhdrn",   &Nhdrn,  "Nhdrn/I");
//   herwigTree->Branch ("Kp",      Kp,      "Kp[Nentry]/I");
//   herwigTree->Branch ("Kn",      Kn,      "Kn[Nentry]/I");
//   herwigTree->Branch ("Stat",    Stat,    "Stat[Nentry]/I");
//   herwigTree->Branch ("Wgt",     &Wgt,    "Wgt/D");
//   herwigTree->Branch ("Qscl",    Qscl,    "Qscl[5]/D");
//   herwigTree->Branch ("Pm",      Pm,      "Pm[Nentry]/D");
}

void StoreParticles::rebind(const TranslationMap & trans) {
  // dummy = trans.translate(dummy);
  AnalysisHandler::rebind(trans);
}

IVector StoreParticles::getReferences() {
  IVector ret = AnalysisHandler::getReferences();
  // ret.push_back(dummy);
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

