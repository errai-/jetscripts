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

#endif // HELP_FUNCTIONS_H