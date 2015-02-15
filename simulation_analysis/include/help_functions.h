#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
// Nice libraries from C
#include <cmath>
#include <ctime>
#include <cstdint>
// ROOT, histogramming
#include "TROOT.h"
#include "TMath.h"
#include "TProfile.h"

using std::cout;
using std::endl;

// Timer class for showing the progress of a simulation

class Timer{
  private:
    std::clock_t start;
    int nEvent;
    int dEvent; // spacing of time/event lattice
    int cEvent = 0; // current event
    int hours; 
    int minutes; 
    int seconds; 
    
  public:
    Timer(){}
    
    Timer(int _nEvent, int _dEvent): nEvent(_nEvent),dEvent(_dEvent) {}
    
    ~Timer(){}

    void set_params(int _nEvent, int _dEvent){
      nEvent = _nEvent;
      dEvent = _dEvent;
    }
    
    void start_timing(){
      start = std::clock(); 
    }
    
    void count_time(){
      cEvent += dEvent;
      double time_processor = (std::clock() - start)/(( (double) CLOCKS_PER_SEC ) );
      time_processor = time_processor*( ((double) nEvent)/cEvent-1); 
      minutes =  time_processor/60; 
      hours = minutes/60;
      seconds = time_processor-60*minutes;
      minutes = minutes - hours*60;
    }

    void print_time(){
      count_time();
      cout << cEvent << " events created, ETA : " << hours << "h" <<
        minutes << "m" << seconds << "s." << endl;
    }
};


// Difference in R of two jets

static double deltaR( double phi1, double phi2, double eta1, double eta2 ){
  
  double dPhi = abs(phi1-phi2);
  while ( dPhi >  TMath::Pi() ) dPhi -= TMath::Pi();
  
  double dEta = eta1 - eta2;
  
  return pow( pow( dPhi, 2 ) + pow( dEta, 2 ) , 0.5 );
}

// Some functions for pdg particle code recognition:

static int isBottom( int id ) {
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

static int isCharm( int id ) {
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

static int isStrange( int id ) {
  int code1;
  int code2;
  bool tmpHasStrange = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 3 || code2 == 3) tmpHasStrange = true;
  return tmpHasStrange;
}

static int isDown( int id ) {
  int code1;
  int code2;
  bool tmpHasDown = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 2 || code2 == 2) tmpHasDown = true;
  return tmpHasDown;
}

static int isUp( int id ) {
  int code1;
  int code2;
  bool tmpHasUp = false;
  code1 = (int)( ( abs( id ) / 100)%10 );
  code2 = (int)( ( abs( id ) /1000)%10 );
  if ( code1 == 1 || code2 == 1) tmpHasUp = true;
  return tmpHasUp;
}

static int statusCheck( int id1, int id2 ){
  if ( id1 == 5 && isBottom( id2 ) ) return 1;
  if ( id1 == 4 && isCharm( id2 ) ) return 1;
  if ( id1 == 3 && isStrange( id2 ) ) return 1;
  if ( id1 == 2 && isDown( id2 ) ) return 1;
  if ( id1 == 1 && isUp( id2 ) ) return 1;
  return 0;
}


#endif // HELP_FUNCTIONS_H
