#include "Pi0W2.h"

#include <string>

// #include <TGraphErrors.h>
#include <TH1F.h>
#include <TLorentzVector.h>
// #include <TProfile.h>
#include "TCanvas.h"

#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

#include "CommonCalc.h"

using namespace genie;
using namespace genie::mc_vs_data;

Pi0W2::Pi0W2() 
   : GenPlotsBase()
{

   fName = "Pi0_W2";
   
   int nw   = CommonCalc::GetInstance()->GetNW();
   aW2      = new double[nw];
   npizero  = new double[nw];
   
   for ( int i=0; i<nw; ++i )
   {
      aW2[i]     = 0.;
      npizero[i] = 0.;
   }

}

Pi0W2::~Pi0W2()
{

   Clear();
   delete [] aW2;
   delete [] npizero;

}

void Pi0W2::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();
   
   int nw = CommonCalc::GetInstance()->GetNW();
   for ( int i=0; i<nw; ++i )
   {
      aW2[i] = 0.;
      npizero[i] = 0.;
   }
   
   return;

}

void Pi0W2::AnalyzeEvent( const EventRecord& event )
{
      
// For efficiency, it needs to be done only once per event record,
// and there's a flag in CommonCalc that makes this method
// simply return if the event record has been processed already.
// The flag gets reset at the end of running the last histo module
// (also internal machinery, nothing for a user to worry about).
//
   CommonCalc::GetInstance()->AnalyzeEvent( event );
   
   int npi0 = 0;
   double weight = CommonCalc::GetInstance()->GetWeight(); // if it's just 1., why would one need it ???
      
   GHepParticle * p = 0;
   TIter event_iter(&event);   
   while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
   {
      if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) 
      {
         if (p->Pdg() == kPdgPi0) npi0++;
      }     
   }

   int ipos = CommonCalc::GetInstance()->GetCurrentW2Bin();
   if (ipos==-1) return; //out of the energy range

   aW2[ipos] += weight*CommonCalc::GetInstance()->GetCurrentW2();
   npizero[ipos] += weight*npi0;
   
   return;

}

void Pi0W2::EndOfEventLoopPlots()
{
    
   int nw = CommonCalc::GetInstance()->GetNW();
   double* nv = CommonCalc::GetInstance()->GetNV();
   
   for (int i = 0; i<nw; i++)
   {
       if (nv[i]>0.)
       {
	  npizero[i] /= nv[i];
	  aW2[i]     /= nv[i];
//	  std::cout << " aW2[i] = " << aW2[i] << std::endl;
       }
   }
    
   // 
   // TGraphErrors* curgraph = new TGraphErrors(nw-1, aW2,npizero,gerrx,errnpi0); 
   // 
   // All in all, for the uniformity it's better to make it a histo with variable bin size
   //
   // Note: this one below is a valid ctor:  
   //       TH1F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   //
   fResult = new TH1F( "Pi0W2", "pi0 mult vs W2", nw-1, aW2 );
   fResult->SetDirectory(0); // in order to NOT add it to the already open Events TFile (Genie's) 
   
   for ( int i=0; i<nw; ++i )
   {
      fResult->Fill( aW2[i], npizero[i] );
   }
   
   fResult->SetLineColor(kRed);
   fResult->SetMarkerStyle(20);
   fResult->SetMarkerColor(kRed);
   fResult->SetMarkerSize(0.8);
      
   return;
}
