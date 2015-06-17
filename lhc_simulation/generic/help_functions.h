#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <chrono>
// Nice libraries from C
#include <cassert>
#include <cmath>
#include <cstdint>
// ROOT, histogramming
#include "TROOT.h"
#include "TMath.h"
#include "TProfile.h"

using std::cout;
using std::endl;
using std::string;
using namespace std::chrono;

// Timer class for showing the progress of a simulation

class Timer{
private:
    high_resolution_clock::time_point mStart;
    int mTotal;
    int mDelta; // spacing of time/event lattice
    int mCurr;  // current event
    int mHours; 
    int mMinutes; 
    int mSeconds;    

public:
    Timer(){}
        
    Timer(int eventNo, int dEvent): mTotal(eventNo), mDelta(dEvent) {}
        
    ~Timer() 
    {
        cout << "Time elapsed: "; 
        cout << duration_cast<seconds>(high_resolution_clock::now()-mStart).count();
        cout << " s." << endl;
    }

    void setParams(int eventNo, int dEvent) 
    {
        mTotal = eventNo;
        mDelta = dEvent;
    }
        
    void startTiming() 
    {
        mCurr = 0;
        mStart = high_resolution_clock::now();
    }
        
    void countTime() 
    {
        mCurr += mDelta;
        high_resolution_clock::time_point now = high_resolution_clock::now();
        double timeProcessor = duration_cast<seconds>(now-mStart).count();
        timeProcessor = timeProcessor*( ((double) mTotal)/mCurr-1); 
        mMinutes =  timeProcessor/60; 
        mHours = mMinutes/60;
        mSeconds = timeProcessor-60*mMinutes;
        mMinutes = mMinutes - mHours*60;
    }

    void printTime() {
        countTime();
        cout << mCurr << " events created, ETA : " << mHours << "h" <<
            mMinutes << "m" << mSeconds << "s." << endl;
    }
};


/* Turn a given string into a suitable root file name */
static string rootFileName(char *inputName) {
    string editString = inputName;
    assert(editString.size());
    if (editString.find(".root",editString.size()-5) != string::npos) {
        return editString;
    }
    editString += ".root";
    return editString;
}


/* Some functions for PDG particle code recognition for hadrons: 
 * (these are basically implementations from CMSSW) */
namespace HadrFuncs
{
    static int HasTop(int id) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 6 || code2 == 6) return true;
        return false;
    }

    static int HasBottom(int id) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 5 || code2 == 5) return true;
        return false;
    }

    static int HasCharm(int id) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 4 || code2 == 4) return true;
        return false;
    }

    static int HasStrange( int id ) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 3 || code2 == 3) return true;
        return false;
    }

    static int HasDown( int id ) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 2 || code2 == 2) return true;
        return false;
    }

    static int HasUp( int id ) 
    {
        int code1 = ( abs( id ) / 100)%10;
        int code2 = ( abs( id ) /1000)%10;
        if ( code1 == 1 || code2 == 1) return true;
        return false;
    }


    /* Checks whether a hadron has the given quark */
    static int StatusCheck( int quarkId, int hadronId )
    {
        switch (quarkId) {
            case 6: return HasTop(hadronId); break;
            case 5: return HasBottom(hadronId); break;
            case 4: return HasCharm(hadronId); break;
            case 3: return HasStrange(hadronId); break;
            case 2: return HasDown(hadronId); break;
            case 1: return HasUp(hadronId); break;
            default: return 0; break;
        }
    }
}

#endif // HELP_FUNCTIONS_H
