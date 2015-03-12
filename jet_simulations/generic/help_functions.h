#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <string>
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
public:
   Timer(){}
    
   Timer(int eventNo, int dEvent): mTotal(eventNo), mDelta(dEvent) {}
    
   ~Timer(){}

   void setParams(int eventNo, int dEvent) 
   {
      mTotal = eventNo;
      mDelta = dEvent;
   }
    
   void startTiming() 
   {
      mCurr = 0;
      mStart = std::clock(); 
   }
    
   void countTime() 
   {
      mCurr += mDelta;
      double timeProcessor = (std::clock() - mStart)/( (double) CLOCKS_PER_SEC);
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

private:
   std::clock_t mStart;
   int mTotal;
   int mDelta; // spacing of time/event lattice
   int mCurr;  // current event
   int mHours; 
   int mMinutes; 
   int mSeconds;    
};


/* Some functions for PDG particle code recognition: 
 * (these are basically implementations from CMSSW) */
static int isBottom(int id) 
{
   int code1;
   int code2;
   bool tmpHasBottom = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
   return tmpHasBottom;
}

static int isCharm(int id) 
{
   int code1;
   int code2;
   bool tmpHasCharm = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
   return tmpHasCharm;
}

static int isStrange( int id ) 
{
   int code1;
   int code2;
   bool tmpHasStrange = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 3 || code2 == 3) tmpHasStrange = true;
   return tmpHasStrange;
}

static int isDown( int id ) 
{
   int code1;
   int code2;
   bool tmpHasDown = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 2 || code2 == 2) tmpHasDown = true;
   return tmpHasDown;
}

static int isUp( int id ) 
{
   int code1;
   int code2;
   bool tmpHasUp = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 1 || code2 == 1) tmpHasUp = true;
   return tmpHasUp;
}


static int chargeSign( int id )
{
   if ( id == 1 ) return 1;
   if ( id == -2 ) return 1;
   if ( id == -3 ) return 1;
   if ( id == 4 ) return 1;
   if ( id == -5 ) return 1;
   if ( id == 6 ) return 1;
   if ( id == -1 ) return -1;
   if ( id == 2 ) return -1;
   if ( id == 3 ) return -1;
   if ( id == -4 ) return -1;
   if ( id == 5 ) return -1;
   if ( id == -6 ) return -1;
   return 1;
}

#endif // HELP_FUNCTIONS_H
