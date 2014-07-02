#ifndef _HAD_Pi0W2_H_
#define _HAD_Pi0W2_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class Pi0W2 : public GenPlotsBase
{

   public:
    
      Pi0W2();
      virtual ~Pi0W2();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventLoopPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
      double* aW2; // = {0.0};
      double* npizero; //      = {0.0};   // no. of charge particles

};


} //mc_vs_data
} //genie


#endif
