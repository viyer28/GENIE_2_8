#ifndef _HAD_ChHadDoverMultW2_H_
#define _HAD_ChHadDoverMultW2_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class ChHadDoverMultW2 : public GenPlotsBase
{

   public:
    
      ChHadDoverMultW2();
      virtual ~ChHadDoverMultW2();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventLoopPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
      double* aW2; // = {0.0};
      double* nch; //      = {0.0};   // no. of charge particles
      double* D;
      double* D_nch;
      double* errnch;

};


} //mc_vs_data
} //genie


#endif
