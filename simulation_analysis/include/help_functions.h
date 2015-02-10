#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
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
#include "TProfile.h"

using std::cout;
using std::endl;

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

static void histFiller( vector<TProfile*> &hists, double pt, double eTot, double piPlus,
  double piMinus, double pi0Gamma, double kaPlus, double kaMinus, double kSZero,
  double kLZero, double proton, double aproton, double neutron, double aneutron,
  double gamma, double lambda0, double sigma, double elecmuon, double others ){
  hists[0]->Fill( pt, piPlus/eTot ); hists[1]->Fill( pt, piMinus/eTot );
  hists[2]->Fill( pt, pi0Gamma/eTot ); hists[3]->Fill( pt, kaPlus/eTot );
  hists[4]->Fill( pt, kaMinus/eTot ); hists[5]->Fill( pt, kSZero/eTot );
  hists[6]->Fill( pt, kLZero/eTot ); hists[7]->Fill( pt, proton/eTot );
  hists[8]->Fill( pt, aproton/eTot ); hists[9]->Fill( pt, neutron/eTot );
  hists[10]->Fill( pt, aneutron/eTot ); hists[11]->Fill( pt, gamma/eTot );
  hists[12]->Fill( pt, lambda0/eTot ); hists[13]->Fill( pt, sigma/eTot );
  hists[14]->Fill( pt, elecmuon/eTot ); hists[15]->Fill( pt, others/eTot );
}


#endif // HELP_FUNCTIONS_H
