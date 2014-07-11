#ifndef Val_HAD_ExpData_H
#define Val_HAD_ExpData_H

#include <vector>
#include <string>
#include <iostream>

#include "libxml/xmlreader.h"

// class TGraph;
// class TGraphErrors;
#include "TGraphErrors.h"

namespace genie {
namespace mc_vs_data {

class ExpGraph 
{

   public:

      ExpGraph() : fGraph(0), fXLog(false), fYLog(false) {}
      ExpGraph(Int_t n, const Float_t *x, const Float_t *y, const Float_t *ex=0, const Float_t *ey=0) { fGraph=new TGraphErrors(n,x,y,ex,ey); 
                                                                                                        fXLog=false;
												        fYLog=false; }
      ExpGraph(Int_t n, const Double_t *x, const Double_t *y, const Double_t *ex=0, const Double_t *ey=0) { fGraph=new TGraphErrors(n,x,y,ex,ey); 
                                                                                                            fXLog=false;
												            fYLog=false; }
      ExpGraph(const ExpGraph& gr ) { fGraph=new TGraphErrors(*(gr.fGraph)); fXLog=gr.fXLog; fYLog=gr.fYLog; }
      ~ExpGraph() { if ( fGraph ) delete fGraph; }

      void SetXLog( const std::string& xlog ) { if ( xlog == "ON" || xlog == "YES" ) fXLog = true; return; }
      void SetYLog( const std::string& ylog ) { if ( ylog == "ON" || ylog == "YES" ) fYLog = true; return; }
   
      TGraphErrors* GetGraph() { return fGraph; }
      bool          GetXLog()  const { return fXLog; }
      bool          GetYLog()  const { return fYLog; } 
   
   private:
   
      TGraphErrors* fGraph;
      bool          fXLog;
      bool          fYLog; 

};

/* JVY: just some thoughts ahead... maybe will update the storage structure later on... although it can mess up on the MC side...

class Observable
{

   
   public:
      
      Observable() : fName(""), fXLog(false), fYLog(false), fXMin(0.), fXMax(0.), fYMin(0.), fYMax(0.) {}
      Observable( const Observable& var ) { fName=var.fName; fXLog=var.fXLog; fYLog=var.fYLog; 
                                            fXMin=var.fXMin; fXMax=var.fXMax; 
					    fYMin=var.fYMin; fYMax=var.fYMax; } 
      bool operator<( const Observable& var ) { if ( fName < var.fName ) return true; return false; }
   
   private:
   
      std::string fName;
      bool        fXLog;
      bool        fYLog;
      double      fXMin;
      double      fXMax;
      double      fYMin;
      double      fYMax;

};
*/

class ExpData
{
   
   public:
   
      // FIXME !!!
      //
      // At present, we treat CHORUS target (composite/mix) as a special case.
      // We assume that it's an exception, and there wont be any other dataset
      // of that kind in the near future.
      // This is not an ideal solution.However, implementing a machinery 
      // to identify and store a composite target info is a complex task,
      // and it'll require a justification before investing any manpower.
      // Besides, GENIE output does NOT record this kind of mata-data either.
      //
      enum InteractionType { kInvalid=-1, kNuP=0,    kNuN=1,    kNuIon=2,    kNuCHORUS=3,
                                          kNubarP=4, kNubarN=5, kNubarIon=6, kNubarCHORUS=7 };
      
      ExpData();
      ~ExpData();
      
      bool LoadExpData( const std::string& );
            
      const std::map< std::string, std::vector<std::string> >*    GetExpDataNames( const InteractionType& ) const;
      const std::map< std::string, std::vector<ExpGraph*> >*      GetExpDataGraphs( const InteractionType& ) const;
      
      bool Exists( const InteractionType&, const std::string& ) const ;
   
   private:
   
      bool ProcessRecord( xmlTextReader* );
      InteractionType CheckInteractionType( const xmlChar* , const xmlChar* );
      void AddExpData( const InteractionType&, const std::string& );
      ExpGraph*     MakeGraph( const InteractionType&, const std::string& );
      
      std::string     fExpDataDirPath;
      InteractionType fCurrentIntType;
      std::string     fCurrentDSLocation;
      std::string     fCurrentDSReference;
      ExpGraph*       fCurrentGraph;
      
      // can this be a generic approach ???
      // so far looks like it works...
      //
      std::map< InteractionType, std::map< std::string, std::vector<ExpGraph*> > > fGraphs;
//      std::map< InteractionType, std::map< std::string, std::vector<TGraphErrors*> > > fGraphs;
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
