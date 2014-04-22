#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>

#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

#include "Analyzer.h"

#include "ExpData.h"

// here there should also be others, for different observables...
#include "GenPlotsChHadW2.h"

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

using namespace genie;
using namespace genie::mc_vs_data;

Analyzer::Analyzer() 
   : fCurrentBeam(-1), fCurrentTarget(-1), fEvtTree(0), fMCRec(0), 
     fInteractionType(ExpData::kInvalid) 
{ 

   fGenPlots.push_back( new GenPlotsChHadW2() );
   fOutputFormat = "gif"; 

}

Analyzer::~Analyzer()
{

   int NCounts = fGenPlots.size();
   for ( int nc=0; nc<NCounts; ++nc )
   {
      delete fGenPlots[nc];
   }
   fGenPlots.clear();

}

void Analyzer::Analyze( const std::string& model, const std::string& sample )
{
   
// FIXME !!!
// Right now model(genie version) is NOT used anywhere.
// Need to add it up, like a regression test.

   // reset interaction type to default
   //
   fInteractionType = ExpData::kInvalid;

   LOG("gvldtest", pNOTICE) 
     << "Analyzing input file " << sample << " for model: " << model;
   
   TFile* fin = new TFile( sample.c_str(), "READ" );
   if ( !fin )
   {
      LOG("gvldtest", pERROR) 
         << "Invalid Genie sample " << sample;
      fInteractionType = ExpData::kInvalid;
      return;
   }
   
   fEvtTree = dynamic_cast <TTree *>( fin->Get("gtree")  );
   
   if (!fEvtTree) 
   {
      LOG("gvldtest", pERROR) 
         << "Null input GHEP event tree found in file: " << sample;
      fInteractionType = ExpData::kInvalid;
      return;
   }

   fEvtTree->SetBranchAddress("gmcrec", &fMCRec);   

   if (!fMCRec) 
   {
      LOG("gvldtest", pERROR) << "Null MC record";
      fInteractionType = ExpData::kInvalid;
      return;
   }
   
   Long64_t NEvents = fEvtTree->GetEntries();
   if ( NEvents <= 0 ) 
   {
      LOG("gvldtest", pERROR) << "Number of events = 0";
      fInteractionType = ExpData::kInvalid;
      return;
   }
   
   LOG("gvldtest", pNOTICE) 
     << "*** Analyzing: " << NEvents << " events found in file " << sample;
   

   if ( !CheckInteractionType(NEvents) ) 
   {
      LOG("gvldtest", pERROR) << " INVALID interaction type " ;
      fInteractionType = ExpData::kInvalid;
      return; // in fact, it needs to do more than just a return...
   }
      
   int NCounts = 0; // number of plotters for specific observables (size of fGenPlots container)
   
   // clear the plots (gen sample) from previous analysis round (if any)
   //  
   if ( !fGenPlots.empty() )
   {
      NCounts = fGenPlots.size();
      for ( int nc=0; nc<NCounts; ++nc )
      {
         fGenPlots[nc]->Clear();
      }
   }
   
   // loop over Genie event records
   //
   for(Long64_t iev = 0; iev < NEvents; ++iev) 
   {
      
      fEvtTree->GetEntry(iev);

      EventRecord&   event      = *(fMCRec->event);
      
      LOG("gvldtest", pDEBUG) << event;
      
      // go further only if the event is physical
      //
      bool EvtIsUnphysical = event.IsUnphysical();
      if( EvtIsUnphysical ) {
	LOG("gvldtest", pNOTICE) << "Skipping unphysical event";
	fMCRec->Clear();
	continue;
      }

      NCounts = fGenPlots.size();
      for ( int nc=0; nc<NCounts; ++nc )
      {
         fGenPlots[nc]->AnalyzeEvent( event );
      }
      
   } // end loop over Genie event records
      
   // some of the plots are actually produced as graphs
   //
   NCounts = fGenPlots.size();
   for ( int nc=0; nc<NCounts; ++nc )
   {
      fGenPlots[nc]->EndOfEventPlots();
   }
      
   fin->Close();
   delete fin;

   return;

}

bool Analyzer::CheckInteractionType( const int NEvt ) 
{

   ExpData::InteractionType type = ExpData::kInvalid;
   
   assert (fEvtTree);
   assert(fMCRec);
   
   for ( Long64_t i=0; i<NEvt; ++i )
   {

      fEvtTree->GetEntry(i);

      EventRecord&   event      = *(fMCRec->event);
            
      // go further only if the event is physical
      //
      bool EvtIsUnphysical = event.IsUnphysical();
      if( EvtIsUnphysical ) {
	LOG("gvldtest", pNOTICE) << "Skipping unphysical event";
	fMCRec->Clear();
	continue;
      }

      // input particles and primary (???) lepton
      //
      GHepParticle* probe = event.Probe();
      assert(probe);
      GHepParticle* target = event.Particle(1); 
      assert(target);

      if (probe->Pdg() == kPdgNuMu)
      {
	    if      (target->Pdg() == kPdgProton ) { type = ExpData::kNuP; }
	    else if (target->Pdg() == kPdgNeutron) { type = ExpData::kNuN; }
      }
      else if (probe->Pdg() == kPdgAntiNuMu)
      {
	    if      (target->Pdg() == kPdgProton ) { type = ExpData::kNubarP; }
	    else if (target->Pdg() == kPdgNeutron) { type = ExpData::kNubarN; }
      }
      if ( type < 0 )
      {
	 LOG("gvldtest", pERROR) 
           << "Unexpected interaction: neutrino = " << probe->Pdg()
           << " target = " << target->Pdg();
	 return false;
      }
      if ( fInteractionType < 0 )
      {
         
	 fInteractionType = type; 
	      
	 // store probe+target combo 
         //
         fCurrentBeam   = probe->Pdg();
         fCurrentTarget = target->Pdg();
      }
      else
      {
         if ( type != fInteractionType )
	 {
	    LOG("gvldtest", pERROR) 
               << "Interaction type is inconsistent with earlier settings. Abort this sample." ;
	    return false;
	 }               
      }
   }
   
   return true;

}

void Analyzer::DrawResults( const ExpData* dsets )
{
      
   const std::map< std::string, std::vector<TGraphErrors*> >* ExpGraphs = dsets->GetExpDataGraphs( fInteractionType );
   std::map< std::string, std::vector<TGraphErrors*> >::const_iterator itr;
   
   std::stringstream ss;
   ss << fInteractionType;

   int NCounts = fGenPlots.size();
   
   for ( int nc=0; nc<NCounts; ++nc )
   {
      
      TCanvas cnv("cnv","",500,450); 
      cnv.SetLogx(); 
          
      TH1F* plot = fGenPlots[nc]->GetPlot();
      plot->SetStats(0);
      plot->Draw("p");
     
      std::string         name = fGenPlots[nc]->GetName();      

      std::string output = name + "-" + ss.str() + "." + fOutputFormat.c_str();
      
      if ( !ExpGraphs )
      {
         LOG("gvldtest", pNOTICE) << " NO EXP.DATA FOUND IN STORAGE FOR INTERACTION TYPE " << fInteractionType;
         cnv.Print( output.c_str() );
	 continue;
      }

      itr = ExpGraphs->find(name);
      if ( itr == ExpGraphs->end() )
      {
         LOG("gvldtest", pNOTICE) << " NO EXP.DATA FOUND IN STORAGE FOR THIS OBSERVALE " << name; 
         cnv.Print( output.c_str() );
	 continue;
      }
      
      int NDataSets = (itr->second).size();
      for ( int nds=0; nds<NDataSets; ++nds )
      {
         TGraphErrors* gr = (itr->second)[nds];
	 gr->SetMarkerStyle(21+nds);
	 gr->SetMarkerColor(kBlack);
	 gr->SetMarkerSize(0.8);
	 gr->Draw("psame");
      }
            
      cnv.Print(output.c_str()); 
      
   }
   
   return;

}
