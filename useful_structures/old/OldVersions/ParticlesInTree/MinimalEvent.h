#ifndef MINIMALEVENT_H
#define MINIMALEVENT_H

//#include "TROOT.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using std::vector;

typedef union DoubleInBinary{
  double dValue;
  char binary[8];
} DoubleInBinary;

typedef union IntInBinary{
  int iValue;
  char binary[4];
} IntInBinary;

// This file format is designed for the use of ptcut, so that
// a minimal amount of storage space is consumed. For applications
// like jetsorter, it is adviced to use an indicator byte instead
// of saving daughter & mother codes etc. Charges and other discrete
// properties can be easily saved in a binary form.

// Data is saved in binary form without any spaces or other formatting
class MinimalEvent {
private:
  // Editor params
  DoubleInBinary transformer;
  IntInBinary iTransformer;
public:
  // The four momentum of a particle
  vector<double>px;
  vector<double>py;
  vector<double>pz;
  vector<double>e;
  // Id
  vector<int> id;
  size_t eventno;
  // Number of particles
  size_t particles;
  
  MinimalEvent() : particles(0) {}
  ~MinimalEvent(){ }

  void SetVals( double _px, double _py, double _pz, double _e, int _id, size_t _eno){
    px.push_back(_px); 
    py.push_back(_py);
    pz.push_back(_pz);
    e.push_back(_e);
    id.push_back(_id);
    // Pi0 photons - this is indicated
    // by changing the actual particle code to 20 (not in any other use).
    particles++;
    eventno = _eno;
  }

  void Nullify(){
    particles = 0;
    px.clear();
    py.clear();
    pz.clear();
    e.clear();
    id.clear();
  }


  void Write(std::ofstream *output){
    iTransformer.iValue = particles;
    output->write(iTransformer.binary,4);
    for (size_t i = 0; i!=particles; ++i){
      // Float valued parameters
      transformer.dValue = px[i];
      output->write(transformer.binary,8);
      transformer.dValue = py[i];
      output->write(transformer.binary,8);
      transformer.dValue = pz[i];
      output->write(transformer.binary,8);
      transformer.dValue = e[i];
      output->write(transformer.binary,8);
      // Integer valued parameters:
      iTransformer.iValue = id[i];
      output->write(iTransformer.binary,4);
    }
  }

  void Read(std::ifstream *input){
    Nullify();
    input->read(iTransformer.binary,4);
    particles = iTransformer.iValue;
    for (size_t i = 0; i!=particles; ++i){
      // Float valued params
      input->read(transformer.binary,8);
      px.push_back(transformer.dValue);
      input->read(transformer.binary,8);
      py.push_back(transformer.dValue);
      input->read(transformer.binary,8);
      pz.push_back(transformer.dValue);
      input->read(transformer.binary,8);
      e.push_back(transformer.dValue);
      // Integer valued params
      input->read(iTransformer.binary,4);
      id.push_back(iTransformer.iValue);
    }
  }
};

#endif // MINIMALEVENT_H
