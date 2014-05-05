
#include "RunConfig.h"

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

#include "ExpData.h"

using namespace genie;
using namespace genie::mc_vs_data;
using namespace genie::utils;

RunConfig::RunConfig( int argc, char ** argv )
{

  fLineParser = new CmdLnArgParser( argc, argv );
  fGSimFiles  = new GSimFiles( false, 10 ); // no chain; max 10 models/versions
  fExpData    = new ExpData();

  fCurrentModel=-1;
  fCurrentGSample=-1; 
  
  // get GENIE inputs
  if( fLineParser->OptionExists('g') ) 
  {
     string inputs = fLineParser->ArgAsString('g');
     bool ok = fGSimFiles->LoadFromFile(inputs);
     if(!ok) { 
        LOG("gvldtest", pFATAL) 
          << "Could not read validation program inputs from: " << inputs;
//        PrintSyntax();
        exit(1);
     }
  }

  if ( fLineParser->OptionExists('d') )
  {
     string dsets = fLineParser->ArgAsString('d');
     bool ok = fExpData->LoadExpData( dsets );
     if(!ok) 
     { 
        LOG("gvldtest", pFATAL) 
          << "Could not read validation program datasets from: " << dsets;
//        PrintSyntax();
        exit(1);
     }
  }

  if ( fGSimFiles->NModels() > 0 ) 
  {
     fCurrentModel = 0;
     // --> if ( (fGSimFiles->EvtFileNames(fCurrentModel)).size() > 0 ) fCurrentGSample = 0;
  }
  
  fOutputFormat = "gif";
  if(fLineParser->OptionExists('f')) 
  {
    fOutputFormat = fLineParser->ArgAsString('f');
  }  
     
}

RunConfig::~RunConfig()
{

   delete fLineParser;
   delete fGSimFiles;

}

bool RunConfig::IsDone()
{
   
   if ( fCurrentModel == -1 ) return true;
   return false;

}

void RunConfig::Next()
{

   int EndSample = (fGSimFiles->EvtFileNames(fCurrentModel)).size()-1;
   if ( fCurrentGSample == EndSample )
   {
      // we're at the end of list of generated samples for a given model
      //
      fCurrentModel++;
//      fCurrentGSample = 0;
      fCurrentGSample = -1; // 0 is causing a skip of the 1st sample in the next model's list
                            // as an alternative, I can set it at 0 and make it return here;
   }

   if ( fCurrentModel == fGSimFiles->NModels() ) 
   {
     // we're at past the end of the config
     //
     fCurrentModel = -1;
     fCurrentGSample = -1;
     return;
   } 

   // increment sample counter within a given model
   //
   fCurrentGSample++;  

   return;

}

std::string RunConfig::GetCurrentModelName() const
{

// FIXME !!!
// need to assert if valid !!!

   return fGSimFiles->ModelTag( fCurrentModel ) ;

}
std::string& RunConfig::GetCurrentGSampleName() const
{

   // FIXME !!!
   // need to asset it's all valid !!!
   //
   return  (fGSimFiles->EvtFileNames(fCurrentModel))[fCurrentGSample];

}

