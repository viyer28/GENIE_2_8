#ifndef Val_HAD_ExpData_H
#define Val_HAD_ExpData_H

#include <vector>
#include <string>
#include <iostream>

#include "libxml/xmlreader.h"

class TGraph;
class TGraphErrors;

namespace genie {
namespace mc_vs_data {

class ExpData
{
   
   public:
   
      enum InteractionType { kInvalid=-1, kNuP=0, kNuN=1, kNubarP=2, kNubarN=3 };
      
      ExpData();
      ~ExpData();
      
      bool LoadExpData( const std::string& );
            
      const std::map< std::string, std::vector<std::string> >*    GetExpDataNames( const InteractionType& ) const;
      const std::map< std::string, std::vector<TGraphErrors*> >*  GetExpDataGraphs( const InteractionType& ) const; 
   
   private:
   
      bool ProcessRecord( xmlTextReader* );
      InteractionType CheckInteractionType( const xmlChar* , const xmlChar* );
      void AddExpData( const InteractionType&, const std::string& );
      TGraphErrors* MakeGraph( const InteractionType&, const std::string& ); 
      // void MakeGraph( const InteractionType&, const Observable&, const std::string& );
      
      std::string     fExpDataDirPath;
      InteractionType fCurrentIntType;
      std::string     fCurrentDSLocation;
      std::string     fCurrentDSReference;
      TGraphErrors*   fCurrentGraph;
      
      // can this be a generic approach ???
      // so far looks like it works...
      //
      std::map< InteractionType, std::map< std::string, std::vector<TGraphErrors*> > > fGraphs;
      std::map< InteractionType, std::map< std::string, std::vector<std::string> > >   fExpDataHolder;
      
      // Note (by Marc P.)
      //
      // Nested maps and vectors are not completely insane.
      // In particular, the above cosntruct map < InteractionType, map<Observable,vector<string> > >
      // is OK if the access methid is like this:
      // get(InteractionType) --> map< ..., ... >
      // BUT !!! it's twice the binary search, 
      // and the whole thing is much wider distributed in memory.
      // 
      // It could be more efficient to do like this:
      // map< KEY, vector<string> >
      // where
      // KEY = struct{ InteractionType, Observable } (with "less" operand implemented)
      // if the access method is meant to be like this:
      // getstuff(InteractionType,Observable) ---> vector<string>
      //
      // In case of C++11, one can also go for an unordered_map (hash map)
      // which uses the hash function instead of binary search, as map does

};

} // end namespace mc_vs_data
} // end namespace genie


#endif
