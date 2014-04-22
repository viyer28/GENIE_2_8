#ifndef _HAD_GenPlotsBase_H_
#define _HAD_GenPlotsBase_H_

#include <vector>
#include <string>
#include <map>

#include "EVGCore/EventRecord.h"

class TH1F;
// class TProfile;

namespace genie {
namespace mc_vs_data {

class GenPlotsBase
{

   public:
    
      GenPlotsBase() : fName(""), fResult(0) {}
      virtual ~GenPlotsBase();
       
      std::string         GetName() const { return fName; }
      virtual void        Clear();
      virtual void        AnalyzeEvent( const EventRecord& ) = 0;
      virtual void        EndOfEventPlots() = 0;
      
      TH1F*               GetPlot() { return fResult; }

   protected:
         
      std::string                             fName;
      TH1F*                                   fResult;              

};


} //mc_vs_data
} //genie


#endif
