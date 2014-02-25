#ifndef Val_HAD_RunConfig_H
#define Val_HAD_RunConfig_H

#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>
#include <string>

#include "Utils/GSimFiles.h"
#include "Utils/CmdLnArgParser.h"


namespace genie {
namespace mc_vs_data {

class RunConfig
{

   public:
    
      RunConfig( int argc, char ** argv );
      ~RunConfig();
         
      void                 Next();
      bool                 IsDone();
      int                  GetCurrentModelId()     const { return fCurrentModel; }
      int                  GetCurrentGSampleId()   const { return fCurrentGSample; }
      std::string          GetCurrentModelName()   const;  
      std::string&         GetCurrentGSampleName() const;
      std::string          GetOutputFormat()       const { return fOutputFormat; }

   private:
     
     // standard Genie utility classes
     //
     CmdLnArgParser* fLineParser;
     GSimFiles*      fGSimFiles;
     //
     int             fCurrentModel;
     int             fCurrentGSample;
     std::string     fOutputFormat;
     
};

} // mc_vs_data
} // genie

#endif
