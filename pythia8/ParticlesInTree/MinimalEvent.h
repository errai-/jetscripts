#ifndef MINIMALEVENT_H
#define MINIMALEVENT_H

#include "TROOT.h"
#include <vector>
#include <iostream>

class MinimalEvent {
public:
  // The four momentum of a particle
  vector<double>px;
  vector<double>py;
  vector<double>pz;
  vector<double>e;
  vector<int> status;
  vector<int> id;
  vector<int> mother1;
  vector<int> mother2;
  vector<int> daughter1;
  vector<int> daughter2;
  size_t idx;
  
  MinimalEvent() : idx(0) {}
  ~MinimalEvent(){ }

  void SetVals( double _px, double _py, double _pz, double _e, int _status, 
    int _id, int _m1, int _m2, int _d1, int _d2 ){
    px.push_back(_px); 
    py.push_back(_py);
    pz.push_back(_pz);
    e.push_back(_e);
    status.push_back(_status);
    id.push_back(_id);
    mother1.push_back(_m1);
    mother2.push_back(_m2);
    daughter1.push_back(_d1);
    daughter2.push_back(_d2);
    idx++;
  }

  void Nullify(){
    idx = 0;
    px.clear();
    py.clear();
    pz.clear();
    e.clear();
    status.clear();
    id.clear();
    mother1.clear();
    mother2.clear();
    daughter1.clear();
    daughter2.clear();
  }
};

#endif // MINIMALEVENT_H
