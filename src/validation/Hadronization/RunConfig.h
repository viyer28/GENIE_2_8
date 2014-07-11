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

class ExpData;

class RunConfig
{

   public:
    
      RunConfig( int argc, char ** argv );
      ~RunConfig();
         
      void                 Next();
      bool                 IsDone();
      int                  GetNModels()            const { return fGSimFiles->NModels(); }
      int                  GetCurrentModelId()     const { return fCurrentModel; }
      int                  GetCurrentGSampleId()   const { return fCurrentGSample; }
      std::string          GetCurrentModelName()   const;  
      std::string&         GetCurrentGSampleName() const;
      std::string          GetOutputFormat()       const { return fOutputFormat; }
//      const ExpData*       GetExpData()            const { return fExpData; }
      ExpData*       GetExpData()            const { return fExpData; }

   private:
     
     // standard Genie utility classes
     //
     CmdLnArgParser* fLineParser;
     GSimFiles*      fGSimFiles;
     ExpData*        fExpData;
     
     //
     int             fCurrentModel;
     int             fCurrentGSample;
     std::string     fOutputFormat;
     
};

} // mc_vs_data
} // genie

#endif
