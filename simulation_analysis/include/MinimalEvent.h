#ifndef MINIMALEVENT_H
#define MINIMALEVENT_H

//#include "TROOT.h"
#include <vector>
#include <iostream>
#include <fstream>
// #include <cstdint>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TArrayD.h"

using std::vector;

// This file format is designed for the use of ptcut, so that
// a minimal amount of storage space is consumed. For applications
// like jetsorter, it is adviced to use an indicator byte instead
// of saving daughter & mother codes etc. Charges and other discrete
// properties can be easily saved in a binary form.

// Data is saved in binary form without any spaces or other formatting
class MinimalEvent {
public:
  // The four momentum of a particle
  vector<double>px;
  vector<double>py;
  vector<double>pz;
  vector<double>e;
  // Id
  vector<int> id;
  // Number of particles
  size_t particles;
  
  MinimalEvent() : particles(0) {}
  virtual ~MinimalEvent(){ }

  void SetVals( double, double, double, double, int );

  void Nullify();

  TArrayD Pt() const;

  TArrayD Eta() const;

  ClassDef(MinimalEvent, 1);
};

#endif // MINIMALEVENT_H
