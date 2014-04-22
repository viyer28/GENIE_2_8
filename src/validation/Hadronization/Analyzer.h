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
      void Analyze( const std::string& model, const std::string& sample );
      void DrawResults( const ExpData* );
   
   private:
   
      // methods 
      //
      bool CheckInteractionType( const int );
      // void DrawResults();
      
      // data members
      //
      int                        fCurrentBeam;
      int                        fCurrentTarget;
      TTree*                     fEvtTree;
      NtpMCEventRecord*          fMCRec;
      ExpData::InteractionType   fInteractionType; //interaction type, 0:vp, 1:vn, 2:vbarp, 3:vbarn
      std::vector<GenPlotsBase*> fGenPlots;
      std::string                fOutputFormat;

};

} // end namespace mc_vs_data
} // end namespace genie

#endif
