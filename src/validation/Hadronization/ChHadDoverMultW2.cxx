#include "ChHadDoverMultW2.h"

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

ChHadDoverMultW2::ChHadDoverMultW2() 
   : GenPlotsBase()
{

   fName = "ChHad_D_W2";
   
   int nw = CommonCalc::GetInstance()->GetNW();
   aW2    = new double[nw];
   nch    = new double[nw];
   D      = new double[nw];
   D_nch  = new double[nw];
   errnch = new double[nw];
   
   for ( int i=0; i<nw; ++i )
   {
      aW2[i]    = 0.;
      nch[i]    = 0.;
      D[i]      = 0.;
      D_nch[i]  = 0.;
      errnch[i] = 0.;
   }

}

ChHadDoverMultW2::~ChHadDoverMultW2()
{

   Clear();
   delete [] aW2;
   delete [] nch;
   delete [] D;
   delete [] D_nch;
   delete [] errnch;

}

void ChHadDoverMultW2::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();
   
   int nw = CommonCalc::GetInstance()->GetNW();
   for ( int i=0; i<nw; ++i )
   {
      aW2[i]    = 0.;
      nch[i]    = 0.;
      D[i]      = 0.;
      D_nch[i]  = 0.;
      errnch[i] = 0.;
   }
   
   return;

}

void ChHadDoverMultW2::AnalyzeEvent( const EventRecord& event )
{
      
// For efficiency, it needs to be done only once per event record,
// and there's a flag in CommonCalc that makes this method
// simply return if the event record has been processed already.
// The flag gets reset at the end of running the last histo module
// (also internal machinery, nothing for a user to worry about).
//
   CommonCalc::GetInstance()->AnalyzeEvent( event );
   
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

   int ipos = CommonCalc::GetInstance()->GetCurrentW2Bin();
   if (ipos==-1) return; //out of the energy range

   aW2[ipos] += weight*CommonCalc::GetInstance()->GetCurrentW2();
   nch[ipos] += weight*(ncharged);
   errnch[ipos]  += weight*pow((double)ncharged,2);
   
   return;

}

void ChHadDoverMultW2::EndOfEventLoopPlots()
{
    
   int nw = CommonCalc::GetInstance()->GetNW();
   double* nv = CommonCalc::GetInstance()->GetNV();
   
   for (int i = 0; i<nw; i++)
   {
       if (nv[i]>0.)
       {
	  nch[i]   /= nv[i];
	  aW2[i]   /= nv[i];
	  D[i]      = sqrt(errnch[i]/nv[i]-pow(nch[i],2));
	  D_nch[i]  = D[i]/nch[i];
       }
   }
    
   // original graph:
   // TGraph* curgraph = new TGraph(nw-1,aW2,D_nch) 
   // 
   // All in all, for the uniformity it's better to make it a histo with variable bin size
   //
   // Note: this one below is a valid ctor:  
   //       TH1F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   //
   fResult = new TH1F( "ChHad_D_W2", "ch.hadrons D/<mult> vs W2", nw-1, aW2 );
   fResult->SetDirectory(0); // in order to NOT add it to the already open Events TFile (Genie's) 
   
   for ( int i=0; i<nw; ++i )
   {
      fResult->Fill( aW2[i], D_nch[i] );
   }
   
   fResult->SetLineColor(kRed);
   fResult->SetMarkerStyle(20);
   fResult->SetMarkerColor(kRed);
   fResult->SetMarkerSize(0.8);
      
   return;
}
