#ifndef _HAD_CommonCalc_H_
#define _HAD_CommonCalc_H_

#include "EVGCore/EventRecord.h"

namespace genie {
namespace mc_vs_data {

class CommonCalc
{

   public:
   
      ~CommonCalc();
      
      static CommonCalc* GetInstance();
      
      void AnalyzeEvent( const EventRecord& );
      void ResetForNextEvent() { fEventProcessed=false; fCurrentW2=0.; fCurrentW2Bin=-1; return; }
      void Clear() { for ( int i=0; i<fNW; ++i ) { nv[i]=0; nv2[i]=0; } return; }

      
      int GetNW() const { return fNW; }
      double* GetNV() const { return nv; }
      double* GetNV2() const { return nv2; }
      double GetWeight() const { return weight; }
      double GetCurrentW2() const { return fCurrentW2; }
      int GetCurrentW2Bin() const { return fCurrentW2Bin; }
      
   private:
      
      CommonCalc();
   
      int fNW; // # of W2 bins (also called sometimes as "energy bins")
      // naming cipoed over from the original code
      double* W2lo; // = { 1, 2, 3, 4, 6, 8, 12, 16, 23, 32, 45, 63, 90, 125 };
      double* W2hi; // = { 2, 3, 4, 6, 8, 12, 16, 23, 32, 45, 63, 90, 125, 225 };
      double* nv; // = {0}; // # of neutrino interactions, per W2 bin
      double* nv2; //
      double weight;
      //
      double fCurrentW2;
      int fCurrentW2Bin;
      bool fEventProcessed;

      static CommonCalc* fInstance;

};


} //mc_vs_data
} //genie


#endif

