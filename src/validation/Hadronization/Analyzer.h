#ifndef Val_HAD_Analyzer_H
#define Val_HAD_Analyzer_H

#include <string>

#include "ExpData.h"

class TTree;

namespace genie {

// fwd declaration (namespace genie)
class NtpMCEventRecord;

namespace mc_vs_data {

// fwd declaration (namespace genie::mc_vs_data)
class GenPlotsBase;

class Analyzer
{

   public:
   
      Analyzer(); 
      ~Analyzer();
      
      void SetOutputFormat( const std::string& format ) { fOutputFormat=format; return; }
      void SetExpDataPtr( ExpData* ptr ) { fExpDataPtr = ptr; return; }
      void Analyze( const std::string&, const std::string& );
      void DrawResults( const int ); 
   
   private:
   
      // methods 
      //
      bool CheckInteractionType( const int );
      void AddPlots( std::vector<GenPlotsBase*>& );
      
      // data members
      //
      int                        fCurrentBeam;
      int                        fCurrentTarget;
      TTree*                     fEvtTree;
      NtpMCEventRecord*          fMCRec;
      ExpData::InteractionType   fInteractionType; //interaction type, 0:v+p, 1:v+n, 2:vbar+p, 3:vbar+n
      std::map< std::string, std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > > fGenPlots;
      std::string                fOutputFormat;
      ExpData*                   fExpDataPtr;

};

} // end namespace mc_vs_data
} // end namespace genie

#endif
