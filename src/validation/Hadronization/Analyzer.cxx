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

#include "Interaction/Interaction.h"
#include "Interaction/Target.h"

#include "Analyzer.h"

#include "ExpData.h"

// here there should also be others, for different observables...
#include "ChHadW2.h"
#include "NegChHadD.h"
#include "ChHadDoverMultW2.h"
#include "Pi0W2.h"
#include "Pi0D.h"

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "CommonCalc.h"

using namespace genie;
using namespace genie::mc_vs_data;

Analyzer::Analyzer() 
   : fCurrentBeam(-1), fCurrentTarget(-1), fEvtTree(0), fMCRec(0), 
     fInteractionType(ExpData::kInvalid) 
{ 

   fOutputFormat = "gif"; 
   fExpDataPtr = 0 ;

}

Analyzer::~Analyzer()
{

   std::map< std::string, std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > >::iterator itr;
   for ( itr=fGenPlots.begin(); itr!=fGenPlots.end(); ++itr )
   {
      std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >::iterator itr1;
      for ( itr1=(itr->second).begin(); itr1!=(itr->second).end(); ++itr1 )
      {      
         int NCounts1 = (itr1->second).size();
         for ( int nc1=0; nc1<NCounts1; ++nc1 )
         {
            delete (itr1->second)[nc1];
         }
	 (itr1->second).clear();
      }
   } 
   
   // do NOT delete fExpDataPtr !!!
   
}

void Analyzer::AddPlots( std::vector<GenPlotsBase*>& plots )
{

   assert(fExpDataPtr);
   
   if ( fExpDataPtr->Exists( fInteractionType, "ChHad_W2" ) )   plots.push_back( new ChHadW2() );
   if ( fExpDataPtr->Exists( fInteractionType, "NegChHad_D" ) ) plots.push_back( new NegChHadD() );
   if ( fExpDataPtr->Exists( fInteractionType, "ChHad_D_W2ChHad_D_W2" ) ) plots.push_back( new ChHadDoverMultW2() );
   if ( fExpDataPtr->Exists( fInteractionType, "Pi0_W2" ) )     plots.push_back( new Pi0W2() );
   if ( fExpDataPtr->Exists( fInteractionType, "Pi0_D" ) )      plots.push_back( new Pi0D() );

   return;

}

void Analyzer::Analyze( const std::string& model, const std::string& sample )
{
   
   // reset interaction type to default
   //
   fInteractionType = ExpData::kInvalid;
   fCurrentBeam = -1;
   fCurrentTarget = -1;
   
   // clear out common-use arrays
   //
   CommonCalc::GetInstance()->Clear();

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
      
   std::map< std::string, std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > >::iterator itr = fGenPlots.find(model);
   std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >::iterator itr1;
   if ( itr == fGenPlots.end() )
   {
//      // nothing found for this "model"
//      // add plots
      std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > tmp;
      tmp.insert( std::pair< ExpData::InteractionType, std::vector<GenPlotsBase*> >( fInteractionType, std::vector<GenPlotsBase*>() ) );
      itr1 = tmp.find(fInteractionType);
      if ( itr1 != tmp.end() ) // alternatively it needs to abort !
      {
         AddPlots( (itr1->second) );
      } 
      fGenPlots.insert( std::pair< std::string, std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > >( model, tmp ) );
      // reset the iterators properly once more
      itr = fGenPlots.find(model);
      if ( itr == fGenPlots.end() )
      {
         LOG("gvldtest", pERROR) << " INVALID attempt to add analyzer(s) for model " << model << " and interaction type " << fInteractionType;
	 return;
      }
      itr1 = (itr->second).find(fInteractionType);
      if ( itr1 == (itr->second).end() )
      {
         LOG("gvldtest", pERROR) << " INVALID attempt to add analyzer(s) for model " << model << " and interaction type " << fInteractionType;
	 return;
      }
   }
   else 
   {
      // OK, entry for this model (Genie version) is already there;
      // need to expand
      std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >& tmp = itr->second;
      itr1 = tmp.find(fInteractionType);
      if ( itr1 == tmp.end() )
      {
         //... but no entry for this particular interaction type yet
	 // now add a new (sub)entry
	 tmp.insert( std::pair< ExpData::InteractionType, std::vector<GenPlotsBase*> >( fInteractionType, std::vector<GenPlotsBase*>() ) );
         std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >::iterator itr2 = tmp.find(fInteractionType);
	 if ( itr2 != tmp.end() )
	 {
	    AddPlots( (itr2->second) );
	 }
      }      
      itr1 = (itr->second).find(fInteractionType);
      if ( itr1 == (itr->second).end() )
      {
         LOG("gvldtest", pERROR) << " INVALID attempt to add analyzer(s) for model " << model << " and interaction type " << fInteractionType;
	 return;
      }
   }
   
// FIXME !!!
//   if ( !itr1 )
//   {
//      // need to abort !!!
//   }

   int NCounts = 0; // number of plotter per model/inttype
   
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
      
      CommonCalc::GetInstance()->ResetForNextEvent();

      NCounts = (itr1->second).size();
      for ( int nc=0; nc<NCounts; ++nc )
      {
         (itr1->second)[nc]->AnalyzeEvent( event );
      }
      
   } // end loop over Genie event records
      
   // some of the plots are actually produced as graphs
   //
   NCounts = (itr1->second).size();
   for ( int nc=0; nc<NCounts; ++nc )
   {
      (itr1->second)[nc]->EndOfEventLoopPlots();
   }
         
   fin->Close();
   delete fin;

   return;

}

/*
This method relies on the standard PDG coding, including for ions:
PDG2006 convention: 10LZZZAAAI

In Genie, PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z)
There're also convenient methods: 
genie::pdg::IonPdgCodeToZ(int ion_pdgc)
genie::pdg::IonPdgCodeToA(int ion_pdgc)

Names/codes for some commonly used nuclear PDG codes:
const int kPdgTgtFreeP     = 1000010010;
const int kPdgTgtFreeN     = 1000000010;
const int kPdgTgtDeuterium = 1000010020;
const int kPdgTgtC12       = 1000060120;
const int kPdgTgtO16       = 1000080160;
const int kPdgTgtFe56      = 1000260560;

*/

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

      // incident particle ("beam/projectile" neutrino)
      //
      GHepParticle* probe = event.Probe();
      assert(probe);
      
      // Original implementation: in fact, this is a crude assumption that the target is at ipos=1
      //
      // GHepParticle* target = event.Particle(1);
      //
      // A cleaner way to do is as follows.
      // However, one needs to remember that Particle(1) is likely return the PDG 
      // of the nucleon (kPdgProton or kPdgNeutron), even if the target is kPdgTgtFreeP 
      // or kPdgTgtFreeN, while Target will keep it as originally specified. 
      //
      Interaction* interaction = event.Summary();
      Target* target = (interaction->InitState()).TgtPtr();       
      assert(target);
      
      if (probe->Pdg() == kPdgNuMu)
      {
//	    if      ( target->Pdg() == kPdgProton  ) { type = ExpData::kNuP; }
//	    else if ( target->Pdg() == kPdgNeutron ) { type = ExpData::kNuN; }
//	    else if ( target->Pdg() > 1000000000 && target->Pdg() < 2000000000 ) { type = ExpData::kNuIon; }
      
         if      ( target->Z() == 1 && target->A() == 1 ) { type = ExpData::kNuP; }
	 else if ( target->Z() == 0 && target->A() == 1 ) { type = ExpData::kNuN; }
	 else if ( target->Z() > 1 && target->A() > 1 )   { type = ExpData::kNuIon; }
      }
      else if (probe->Pdg() == kPdgAntiNuMu)
      {
//	    if      (target->Pdg() == kPdgProton ) { type = ExpData::kNubarP; }
//	    else if (target->Pdg() == kPdgNeutron) { type = ExpData::kNubarN; }
//	    else if ( target->Pdg() > 1000000000 && target->Pdg() < 2000000000 ) { type = ExpData::kNuIon; }
         if      ( target->Z() == 1 && target->A() == 1 ) { type = ExpData::kNubarP; }
	 else if ( target->Z() == 0 && target->A() == 1 ) { type = ExpData::kNubarN; }
	 else if ( target->Z() > 1 && target->A() > 1 )   { type = ExpData::kNubarIon; }
      }
      if ( type < 0 )
      {
	 LOG("gvldtest", pERROR) 
           << "Unexpected interaction: neutrino = " << probe->Pdg()
           << " target = " << target->Pdg();
	 return false;
      }
      
//      if ( fCurrentTarget > -1 && fCurrentTarget != target->Pdg() )
//      {
//         std::cout << " Different type of TARGET " << fCurrentTarget << " --> " << target->Pdg() << std::endl;
//      }
      
      if ( fInteractionType < 0 )
      {
         
	 fInteractionType = type; 
	      
	 // store probe+target combo 
         //
         fCurrentBeam   = probe->Pdg();
         fCurrentTarget = target->Pdg();
	 continue;
      }
      else
      {
         
	 // NOTE (JVY): We assume that CHORUS case (composite/mix target) is an exception.
	 //             The solution isn't ideal or scalable, but we'll leave with it for now.
	 //             Also, the implementation is somewhat affected by the fact that GENIE
	 //             does NOT record to the output any meta-data that would contain info
	 //             on a composite target. GENIE's output event only stores the actual
	 //             nucleus that was hit (and also the nucleon from it).
	 
	 if ( fInteractionType == ExpData::kNuCHORUS || fInteractionType == ExpData::kNubarCHORUS )
	 {
	    // return true;
	    continue;
	 }
	 
	 
	 if ( type == ExpData::kNuIon || type == ExpData::kNubarIon )
	 {
	    if ( fCurrentTarget != target->Pdg() )
	    {
               if ( probe->Pdg() == kPdgNuMu )
	       {
	          fInteractionType = ExpData::kNuCHORUS;
	       }
	       else if ( probe->Pdg() == kPdgAntiNuMu )
	       {
	          fInteractionType = ExpData::kNubarCHORUS;
	       }
	       // return true;
	       continue;
	    }
	 }
	 	 
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

void Analyzer::DrawResults( const int nmodels ) 
{
      
// FIXME !!!
// We need a security mechanism to ensure that the number of samples/variables 
// specified per model match vs what is specified for another model s!!!


   const std::map< std::string, std::vector<ExpGraph*> >* ExpGraphs; 

   std::map< std::string, std::vector<ExpGraph*> >::const_iterator ditr;

   std::map< std::string, std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> > >::const_iterator gitr;
   std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >::const_iterator gpitr;

   std::map< std::string, std::vector<ExpGraph*> >::const_iterator itr; // iterator over exp.data/graphs
      
   int NModels = fGenPlots.size();
   
   std::map< ExpData::InteractionType, std::vector<TCanvas*> > cnv;
   std::map< ExpData::InteractionType, std::vector<TCanvas*> >::iterator icnv;
   
   int im = 0;
   
   for ( gitr=fGenPlots.begin(); gitr!=fGenPlots.end(); ++gitr ) // loop over models (genie versions)
   {
      ExpData::InteractionType itype = ExpData::kInvalid;
      const std::map< ExpData::InteractionType, std::vector<GenPlotsBase*> >& tmp = (gitr->second);
      for ( gpitr=tmp.begin(); gpitr!=tmp.end(); ++gpitr )
      {
         itype = gpitr->first;
	 //
	 // for the very first model, instiate storage for TCanvas
	 // while other results for models will be overlaid
	 //
         if ( im == 0 ) cnv.insert( std::pair< ExpData::InteractionType, std::vector<TCanvas*> >( itype, std::vector<TCanvas*>() ) );
         icnv = cnv.find(itype); 
	 const std::vector<GenPlotsBase*>& plots = gpitr->second;
	 int NPlots = plots.size();
	 for ( int npl=0; npl<NPlots; ++npl )
	 { 
	   std::string name = plots[npl]->GetName();
	   if ( im == 0 ) // for the very first model, prepare and store TCanvas
	   {
              //
	      // NOTE: TCanvas name must be unique !!!!!
	      //
	      std::stringstream sss;
              sss << itype;
	      std::string cname = name + "-" + sss.str();
	      (icnv->second).push_back( new TCanvas( cname.c_str(), "", 500, 400 ) );
	   }
	   TH1F* plot = plots[npl]->GetPlot();
	   plot->SetStats(0);
           (icnv->second)[npl]->cd();
	   if ( im == 0 ) 
	   {
	      
	      ExpGraphs =  fExpDataPtr->GetExpDataGraphs( itype );
              if ( !ExpGraphs )
              {
                 LOG("gvldtest", pNOTICE) << " NO EXP.DATA FOUND IN STORAGE FOR INTERACTION TYPE " << itype;
	         continue;
              }
              itr = ExpGraphs->find(name);
              if ( itr == ExpGraphs->end() )
              {
                 LOG("gvldtest", pNOTICE) << " NO EXP.DATA FOUND IN STORAGE FOR THIS OBSERVABLE " << name; 
	         continue;
              }      	      
              
	      ExpGraph* egr = (itr->second)[0];
	      if ( egr->GetXLog() ) (icnv->second)[npl]->SetLogx();
	      if ( egr->GetYLog() ) (icnv->second)[npl]->SetLogy();
	      	      
              TGraphErrors* gr = egr->GetGraph();
	      
	      plot->GetXaxis()->SetRangeUser( gr->GetXaxis()->GetXmin(), gr->GetXaxis()->GetXmax() );
	      plot->GetYaxis()->SetRangeUser( gr->GetYaxis()->GetXmin(), gr->GetYaxis()->GetXmax() );
	      plot->GetXaxis()->SetTitle( gr->GetXaxis()->GetTitle() );
	      plot->GetYaxis()->SetTitle( gr->GetYaxis()->GetTitle() );
	      plot->Draw("p");
	   }
	   else
	   {
// Here I'm scaling the same plots from Genie 2.8.0
// just to emulate a "different model" (version)
// 
//	      double scale = 1. + 0.1*im;
//	      plot->Scale( scale );
//	      plot->SetMarkerColor(kGreen);
	      plot->Draw("psame");
	   }
	   if ( im == NModels-1 )
	   {
              std::stringstream ss;
              ss << itype;
	      std::string output = name + "-" + ss.str() + "." + fOutputFormat.c_str();
              int NDataSets = (itr->second).size();
              for ( int nds=0; nds<NDataSets; ++nds )
              {
                 TGraphErrors* gr = ((itr->second)[nds])->GetGraph();
	         gr->SetMarkerStyle(21+nds);
	         gr->SetMarkerColor(kBlack);
	         gr->SetMarkerSize(0.8);
	         gr->Draw("psame");
              }            
              (icnv->second)[npl]->Print(output.c_str());
	   }
	 }
      } // end loop over interaction types within a given model      
      im++;
   } // end loop over models
      
   return;

}
