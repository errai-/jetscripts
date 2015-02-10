#include "MinimalEvent.h"

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

TArrayD MinimalEvent::Pt() const {
  TArrayD pt(particles);
  for (size_t i = 0; i != particles; ++i){
    pt.AddAt(sqrt(pow(px[i],2)+pow(py[i],2)), i);
  }
  return pt;
}

TArrayD MinimalEvent::Eta() const {
  TArrayD eta(particles);
  for (size_t i = 0; i != particles; ++i){
    double P = sqrt(pow(px[i],2) + pow(py[i],2) + pow(pz[i],2));
    eta.AddAt(0.5*log( (P+pz[i])/(P-pz[i])), i);
  }
  return eta;
}
