#ifndef _HAD_Pi0D_H_
#define _HAD_Pi0D_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class Pi0D : public GenPlotsBase
{

   public:
    
      Pi0D();
      virtual ~Pi0D();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventLoopPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
      double* npizero; //      = {0.0};   // no. of charge particles
      double* errnpi0;
      double* Dnpi0;

};


} //mc_vs_data
} //genie


#endif
