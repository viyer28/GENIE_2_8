#ifndef _HAD_NegChHadD_H_
#define _HAD_NegChHadD_H_

//#include <vector>
//#include <string>
//#include <map>

#include "EVGCore/EventRecord.h"

#include "GenPlotsBase.h"

namespace genie {
namespace mc_vs_data {

class NegChHadD : public GenPlotsBase
{

   public:
    
      NegChHadD();
      ~NegChHadD();
       
      virtual void Clear();
      virtual void AnalyzeEvent( const EventRecord& );
      virtual void EndOfEventLoopPlots();

   private:
         
// FIXME ???
// copied over from the original code (with small modifications)
// These variables are meant for kNegChHad_D case

      double* Dneg; // negative hadron dispersion D- vs average multiplicity <n-> (p. 2, 14)
      double* nneg;
      double* errnneg;
};


} //mc_vs_data
} //genie


#endif
