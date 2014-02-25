#ifndef _HAD_GenPlotsChHadW2_H_
#define _HAD_GenPlotsChHadW2_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class GenPlotsChHadW2 : public GenPlotsBase
{

   public:
    
      GenPlotsChHadW2();
      ~GenPlotsChHadW2();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
// These variables are meant for kChHad_W2 case
      int nw;
      double* W2lo; // = { 1, 2, 3, 4, 6,  8, 12, 16, 23, 32, 45, 63,  90, 125 };
      double* W2hi; // = { 2, 3, 4, 6, 8, 12, 16, 23, 32, 45, 63, 90, 125, 225 };
      double* aW2; // = {0.0};
      double* nch; //      = {0.0};   // no. of charge particles
      double* nv;  //       = {0};     // no. of neutrino interactions

};


} //mc_vs_data
} //genie


#endif
