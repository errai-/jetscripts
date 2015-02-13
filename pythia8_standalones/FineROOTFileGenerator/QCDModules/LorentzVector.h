#ifndef LorentzVector_h
#define LorentzVector_h
#include <math.h>
#include <TMath.h>

class Coordinates {
  public:
  
  double fX;
  double fY;
  double fZ;
  double fT;
};

class LorentzVector {
  private:

  Coordinates fCoordinates;

  public:

  double getX() const { return fCoordinates.fX; }
  double getY() const { return fCoordinates.fY; }
  double getZ() const { return fCoordinates.fZ; }
  double getT() const { return fCoordinates.fT; }

  void setX(double _x) { fCoordinates.fX = _x; }
  void setY(double _y) { fCoordinates.fY = _y; }
  void setZ(double _z) { fCoordinates.fZ = _z; }
  void setT(double _t) { fCoordinates.fT = _t; }

  void addX(double _x) { fCoordinates.fX += _x; }
  void addY(double _y) { fCoordinates.fY += _y; }
  void addZ(double _z) { fCoordinates.fZ += _z; }
  void addT(double _t) { fCoordinates.fT += _t; }

  void scaleX(double _x) { fCoordinates.fX *= _x; }
  void scaleY(double _y) { fCoordinates.fY *= _y; }
  void scaleZ(double _z) { fCoordinates.fZ *= _z; }
  void scaleT(double _t) { fCoordinates.fT *= _t; }

  double p() const { return pow( pow(fCoordinates.fX,2) + pow(fCoordinates.fY,2)
    + pow(fCoordinates.fZ,2), 0.5); }
  double pt() const { return pow( pow(fCoordinates.fX,2) + pow(fCoordinates.fY,2)
    , 0.5 ); }
  double eta() const { 
    if ( p() - fCoordinates.fZ == 0 ){ return 1000000000000;
    }else{ return 0.5*TMath::Log( ( p() + fCoordinates.fZ )/( p() - fCoordinates.fZ ) );}
  }
  double phi() const { return TMath::ATan2( fCoordinates.fY, 
    fCoordinates.fX ); }
  double Rapidity() const { return 0.5*log( ( fCoordinates.fT + 
    fCoordinates.fZ )/( fCoordinates.fT - fCoordinates.fZ ) ); }
  double energy() const { return fCoordinates.fT; }
  double mass() const { return pow( pow(fCoordinates.fT,2) - pow( p(), 2 ), 0.5 ); }

  LorentzVector(){ }

  ~LorentzVector(){ }

  int operator=(const LorentzVector& right){
    fCoordinates.fX = right.getX();
    fCoordinates.fY = right.getY();
    fCoordinates.fZ = right.getZ();
    fCoordinates.fT = right.getT();
  }

};

#endif
