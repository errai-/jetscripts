///////////////////////////////////////////////////////////////////
// A class for Herwig to store simulation data into a .root file //
// Based on a template: contains some peculiarities              //
// Hannu Siikonen 7.3.2015                                       //
///////////////////////////////////////////////////////////////////

#ifndef STOREPARTICLES_H
#define STOREPARTICLES_H

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"

#include <iostream>
using std::cout;
using std::endl;

// ROOT:
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "../events/PrtclEvent.h"

namespace jetanalysis {

using namespace ThePEG;

/**
 * Here is the documentation of the StoreParticles class.
 *
 * @see \ref StoreParticlesInterfaces "The interfaces"
 * defined for StoreParticles.
 */
class StoreParticles: public AnalysisHandler {

public:

   /** @name Standard constructors and destructors. */
   //@{
   StoreParticles() {}
   
   ~StoreParticles() {}
   //@}
   
   /** @name Virtual functions required by the AnalysisHandler class. */
   //@{
   /**
    * Analyze a given Event. Note that a fully generated event
    * may be presented several times, if it has been manipulated in
    * between. The default version of this function will call transform
    * to make a lorentz transformation of the whole event, then extract
    * all final state particles and call analyze(tPVector) of this
    * analysis object and those of all associated analysis objects. The
    * default version will not, however, do anything on events which
    * have not been fully generated, or have been manipulated in any
    * way.
    * @param event pointer to the Event to be analyzed.
    * @param ieve the event number.
    * @param loop the number of times this event has been presented.
    * If negative the event is now fully generated.
    * @param state a number different from zero if the event has been
    * manipulated in some way since it was last presented.
    */
   virtual void analyze(tEventPtr event, long ieve, int loop, int state);

   /**
    * Return a LorentzTransform which would put the event in the
    * desired Lorentz frame.
    * @param event a pointer to the Event to be considered.
    * @return the LorentzRotation used in the transformation.
    */
   LorentzRotation transform(tcEventPtr) const { return LorentzRotation(); }

   /**
    * Analyze the given vector of particles. The default version calls
    * analyze(tPPtr) for each of the particles.
    * @param parts the vector of pointers to particles to be analyzed
    * @param weight the weight of the current event.
    */
   void analyze(const tPVector & parts, double) { AnalysisHandler::analyze(parts); }

   /**
    * Analyze the given particle.
    * @param particle pointer to the particle to be analyzed.
    * @param weight the weight of the current event.
    */
   void analyze(tPPtr, double) {}
   //@}


   /**
    * The standard Init function used to initialize the interfaces.
    * Called exactly once for each class by the class description system
    * before the main function starts or
    * when this class is dynamically loaded.
    */
   static void Init()
   {
      static ClassDocumentation<StoreParticles> documentation
         ("There is no documentation for the StoreParticles class");
   }

protected:

   /** @name Clone Methods. */
   //@{
   /**
    * Make a simple clone of this object.
    * @return a pointer to the new object.
    */
   IBPtr clone() const { return new_ptr(*this); }

   /** Make a clone of this object, possibly modifying the cloned object
    * to make it sane.
    * @return a pointer to the new object.
    */
   IBPtr fullclone() const { return new_ptr(*this); }
   //@}

   /** @name Standard Interfaced functions. */
   //@{
   /**
   * Check sanity of the object during the setup phase. Does nothing at the moment.
   */
   void doupdate()
   {
      AnalysisHandler::doupdate();
      bool redo = touched(); if ( !redo ) return;
   }
   
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  void doinit() { AnalysisHandler::doinit(); }

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  void rebind(const TranslationMap & trans) { AnalysisHandler::rebind(trans); }

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences()
  {
  IVector ret = AnalysisHandler::getReferences();
  return ret;
  }
  //@}

  
private:

   /**
    * The assignment operator is private and must never be called.
    * In fact, it should not even be implemented.
    */
   StoreParticles & operator=(const StoreParticles &);

   TTree * herwigTree;
   TFile * herwigFile;

   PrtclEvent *pEvent;
};

}

#endif /* STOREPARTICLES_H */
