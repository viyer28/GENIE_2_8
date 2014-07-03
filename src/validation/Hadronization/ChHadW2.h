#ifndef _HAD_ChHadW2_H_
#define _HAD_ChHadW2_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class ChHadW2 : public GenPlotsBase
{

   public:
    
      ChHadW2();
      virtual ~ChHadW2();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventLoopPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
// These variables are meant for kChHad_W2 case
      double* aW2; // = {0.0};
      double* nch; //      = {0.0};   // no. of charge particles

};


} //mc_vs_data
} //genie


#endif
