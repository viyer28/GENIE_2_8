#include "ChHadW2.h"

#include <string>

// #include <TGraphErrors.h>
#include <TH1F.h>
#include <TLorentzVector.h>

#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

#include "CommonCalc.h"

using namespace genie;
using namespace genie::mc_vs_data;

ChHadW2::ChHadW2() 
   : GenPlotsBase()
{

   fName = "ChHad_W2";
   
   int nw = CommonCalc::GetInstance()->GetNW();
   aW2  = new double[nw];
   nch  = new double[nw];
   
   for ( int i=0; i<nw; ++i )
   {
      aW2[i] = 0.;
      nch[i] = 0.;
   }

}

ChHadW2::~ChHadW2()
{

   Clear();
   delete [] aW2;
   delete [] nch;

}

void ChHadW2::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();
   
   int nw = CommonCalc::GetInstance()->GetNW();
   for ( int i=0; i<nw; ++i )
   {
      aW2[i] = 0.;
      nch[i] = 0.;
   }
   
   return;

}

void ChHadW2::AnalyzeEvent( const EventRecord& event )
{
      
// For efficiency, it needs to be done only once per event record,
// and there's a flag in CommonCalc that makes this method
// simply return if the event record has been processed already.
// The flag gets reset at the end of running the last histo module
// (also internal machinery, nothing for a user to worry about).
//
   CommonCalc* calc = CommonCalc::GetInstance();
   calc->AnalyzeEvent( event );
   
   int ncharged = 0;
   double weight = CommonCalc::GetInstance()->GetWeight(); // if it's just 1., why would one need it ???
      
   GHepParticle * p = 0;
   TIter event_iter(&event);   
   while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
   {
      if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) 
      {
         if (p->Charge()) ncharged++;
      }     
   }

   int ipos = calc->GetCurrentW2Bin();
   if (ipos==-1) return; //out of the energy range

   aW2[ipos] += weight*calc->GetCurrentW2();
   nch[ipos] += weight*(ncharged);
   
   return;

}

void ChHadW2::EndOfEventLoopPlots()
{
    
   int nw = CommonCalc::GetInstance()->GetNW();
   double* nv = CommonCalc::GetInstance()->GetNV();
   
   for (int i = 0; i<nw; i++)
   {
       if (nv[i]>0.)
       {
	  nch[i]   /= nv[i];
	  aW2[i]   /= nv[i];
//	  std::cout << " aW2[i] = " << aW2[i] << std::endl;
       }
   }
    
   // why did the original authour wanted a graph w/errors if the errors are all zeros ???
   // TGraphErrors* curgraph = new TGraphErrors(nw-1,aW2,nch,gerrx,gerrx);   
   // 
   // All in all, for the uniformity it's better to make it a histo with variable bin size
   //
   // Note: this one below is a valid ctor:  
   //       TH1F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   //
   fResult = new TH1F( "ChHadW2", "ch.hadrons mult vs W2", nw-1, aW2 );
   fResult->SetDirectory(0); // in order to NOT add it to the already open Events TFile (Genie's) 
   
   for ( int i=0; i<nw; ++i )
   {
      fResult->Fill( aW2[i], nch[i] );
   }
   
   fResult->SetLineColor(kRed);
   fResult->SetMarkerStyle(20);
   fResult->SetMarkerColor(kRed);
   fResult->SetMarkerSize(0.8);
      
   return;
}
