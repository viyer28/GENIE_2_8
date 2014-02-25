#ifndef Val_HAD_ExpData_H
#define Val_HAD_ExpData_H

#include <vector>
#include <string>
#include <iostream>

class TGraph;
class TGraphErrors;

namespace genie {
namespace mc_vs_data {

class ExpData
{
   
   public:
   
      enum InteractionType { kInvalid=-1, kNuP=0, kNuN=1, kNubarP=2, kNubarN=3 };
      enum Observable      { kNone=-1, kChHad_W2=0, kChHad_D_n=1 };
      
      ExpData();
      ~ExpData() {}
      
      void AddExpData( const InteractionType&, const Observable&, const std::string& );
      const std::map< Observable, std::vector<std::string> >*    GetExpDataNames( const InteractionType& ) const;
      const std::map< Observable, std::vector<TGraphErrors*> >*  GetExpDataGraphs( const InteractionType& ) const; 
   
   private:
   
      void MakeGraph( const InteractionType&, const Observable&, const std::string& );
      //void          DrawGraph( TGraph* gr, int mstyle, int mcol, double msize=0.8, string opt = "def" );
      
      std::string   fExpDataDirPath;
      
      // can this be a generic approach ???
      // so far looks like it works...
      //
      std::map< InteractionType, std::map< Observable, std::vector<TGraphErrors*> > > fGraphs;
      std::map< InteractionType, std::map< Observable, std::vector<std::string> > >   fExpDataHolder;
      
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
