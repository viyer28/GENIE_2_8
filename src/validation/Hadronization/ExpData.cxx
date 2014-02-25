
#include <TSystem.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TLegend.h>
#include <TLatex.h>
#include <Riostream.h>

#include "Utils/SystemUtils.h"
#include "Messenger/Messenger.h"

#include "ExpData.h"

using namespace genie;
using namespace genie::mc_vs_data;

ExpData::ExpData()
{

   std::string GDir = gSystem->Getenv("GENIE"); // (it's a const char* return type ) 
   fExpDataDirPath = GDir + "/data/validation/vA/hadronics/bubble_chambers/";
      
   // default datasets
   // charged hadrons mult vs W^2
   //

   fExpDataHolder[kNuP][kChHad_W2].push_back("charged_hadron_multiplicity/prd27_47_fig4_vp.dat");
   fExpDataHolder[kNuP][kChHad_W2].push_back("charged_hadron_multiplicity/npb181_385_fig1_nch_W.dat");
   MakeGraph( kNuP, kChHad_W2, "charged_hadron_multiplicity/prd27_47_fig4_vp.dat" );
   MakeGraph( kNuP, kChHad_W2, "charged_hadron_multiplicity/npb181_385_fig1_nch_W.dat" );
   
   fExpDataHolder[kNuN][kChHad_W2].push_back("charged_hadron_multiplicity/prd27_47_fig4_vn.dat");
   fExpDataHolder[kNuN][kChHad_W2].push_back("charged_hadron_multiplicity/zpc24_119_vn.dat");
   MakeGraph( kNuN, kChHad_W2, "charged_hadron_multiplicity/prd27_47_fig4_vp.dat" );
   MakeGraph( kNuN, kChHad_W2, "charged_hadron_multiplicity/npb181_385_fig1_nch_W.dat" );

}

void ExpData::AddExpData( const InteractionType& type, const Observable& var, const std::string& dfile_location )
{

   std::map< InteractionType, std::map< Observable, std::vector<std::string> > >::iterator i1 = fExpDataHolder.find(type);   
   
   // check if data set(s) exist(s) for this interaction type
   // if not, just add a new element
   //
   if ( i1 == fExpDataHolder.end() )
   {
      fExpDataHolder[type][var].push_back( dfile_location );
      MakeGraph( type, var, dfile_location );
      return;
   }

   // data for this interaction type do exist
   //
   std::map< Observable, std::vector<std::string> >::iterator i2 = (i1->second).find(var);

   // check if data exist for this observable
   // if not, just make a new entry
   //
   if ( i2 == (i1->second).end() )
   {
      ((i1->second)[var]).push_back( dfile_location );
      MakeGraph( type, var, dfile_location );
      return;
   } 
   
   int NEntries = (i2->second).size();
   for ( int ie=0; ie<NEntries; ++ie )
   {
      if ( (i2->second)[ie] == dfile_location ) return; // don't add this dataset as it's already in !!!
   }
   
   (i2->second).push_back( dfile_location );
   MakeGraph( type, var, dfile_location );

   return;

}

const std::map< ExpData::Observable, std::vector<std::string> >* ExpData::GetExpDataNames( const InteractionType& type ) const
{

   if ( type == kInvalid ) return NULL;
   
   std::map< InteractionType, std::map< Observable, std::vector<std::string> > >::const_iterator i = fExpDataHolder.find(type);

   if ( i != fExpDataHolder.end() )
   {
      return &(i->second);
   }

   return NULL;

}

const std::map< ExpData::Observable, std::vector<TGraphErrors*> >* ExpData::GetExpDataGraphs( const InteractionType& type ) const
{

   if ( type == kInvalid ) return NULL;
   
   std::map< InteractionType, std::map< Observable, std::vector<TGraphErrors*> > >::const_iterator i = fGraphs.find(type);
   
   if ( i != fGraphs.end() )
   {
      return &(i->second);
   }
   
   return NULL;

}

// TGraphErrors* ExpData::MakeGraph( const InteractionType& type, const Observable& var, std::string& file )
void ExpData::MakeGraph( const InteractionType& type, const Observable& var, const std::string& file )
{

  // fragment copied over from the original code, 
  // (with just a couple of small modifications)
  //
  string filename = fExpDataDirPath + file;

  bool file_exists = !(gSystem->AccessPathName(filename.c_str()));
  if(!file_exists) {
    LOG("gvldtest", pERROR) << "Can not find file: " << filename;
    return ;
  }

  ifstream in;
  in.open(filename.c_str());
  double x, y, ex, ey;
  vector<double> vx;
  vector<double> vy;
  vector<double> vex;
  vector<double> vey;
  while (1) {
    in >>x >> y >>ey;
    ex = 0;
    if (!in.good()) break;
    vx.push_back(x);
    vy.push_back(y);
    vex.push_back(ex);
    vey.push_back(ey);
  }

  fGraphs[type][var].push_back( new TGraphErrors(vx.size(),&vx[0],&vy[0],&vex[0],&vey[0]) ); 
    
  return;
   
}
