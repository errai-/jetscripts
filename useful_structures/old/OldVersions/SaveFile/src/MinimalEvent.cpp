#include "../include/MinimalEvent.h"

void MinimalEvent::SetVals( double _px, double _py, double _pz, double _e, int _id){
  px.push_back(_px); 
  py.push_back(_py);
  pz.push_back(_pz);
  e.push_back(_e);
  id.push_back(_id);
  // Pi0 photons - this is indicated
  // by changing the actual particle code to 20 (not in any other use).
  particles++;
}

void MinimalEvent::Nullify(){
  particles = 0;
  px.clear();
  py.clear();
  pz.clear();
  e.clear();
  id.clear();
}

void MinimalEvent::Write(std::ofstream *output){
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

void MinimalEvent::Read(std::ifstream *input){
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