#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

// Stdlib header file for input and output.
#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <cassert>
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
using std::string;

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


/* Turn a given string into a suitable root file name */
string rootFileName(char *inputName) {
   string editString = inputName;
   assert(editString.size());
   if (editString.find(".root",editString.size()-5) != string::npos) {
      return editString;
   }
   editString += ".root";
   return editString;
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


/* Helper function for isExcitedState */
static int statusCheck( int quarkId, int hadronId )
{
   switch (quarkId) {
      case 6: return hasTop(hadronId);
      case 5: return hasBottom(hadronId);
      case 4: return hasCharm(hadronId);
      case 3: return hasStrange(hadronId);
      case 2: return hasDown(hadronId);
      case 1: return hasUp(hadronId);
      default: return 0;
   }
}


/* Some functions for PDG particle code recognition: 
 * (these are basically implementations from CMSSW) */
static int hasTop(int id) 
{
   int code1;
   int code2;
   bool tmpHasTop = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 6 || code2 == 6) tmpHasTop = true;
   return tmpHasTop;
}

static int hasBottom(int id) 
{
   int code1;
   int code2;
   bool tmpHasBottom = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
   return tmpHasBottom;
}

static int hasCharm(int id) 
{
   int code1;
   int code2;
   bool tmpHasCharm = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
   return tmpHasCharm;
}

static int hasStrange( int id ) 
{
   int code1;
   int code2;
   bool tmpHasStrange = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 3 || code2 == 3) tmpHasStrange = true;
   return tmpHasStrange;
}

static int hasDown( int id ) 
{
   int code1;
   int code2;
   bool tmpHasDown = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 2 || code2 == 2) tmpHasDown = true;
   return tmpHasDown;
}

static int hasUp( int id ) 
{
   int code1;
   int code2;
   bool tmpHasUp = false;
   code1 = (int)( ( abs( id ) / 100)%10 );
   code2 = (int)( ( abs( id ) /1000)%10 );
   if ( code1 == 1 || code2 == 1) tmpHasUp = true;
   return tmpHasUp;
}


#endif // HELP_FUNCTIONS_H
