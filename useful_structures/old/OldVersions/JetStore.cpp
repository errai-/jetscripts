
#include "../include/JetStore.h"


void JetStore::Insert(double px, double py, double pz, double e, int _status, int _id){
  fastjet::PseudoJet tmpJet(px,py,pz,e);
  Momentum.push_back(tmpJet);
  status.push_back(_status);
  id.push_back(_id);
}

void JetStore::Clear(){
  Momentum.clear();
  status.clear();
  id.clear();
}