#include "Pi0D.h"

#include <string>

#include <TGraphErrors.h>
#include <TH1F.h>
#include "TCanvas.h"

#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

#include "CommonCalc.h"

using namespace genie;
using namespace genie::mc_vs_data;

Pi0D::Pi0D() 
   : GenPlotsBase()
{

   fName = "Pi0_D";
   
   int nw   = CommonCalc::GetInstance()->GetNW();
   npizero  = new double[nw];
   errnpi0  = new double[nw];
   Dnpi0    = new double[nw];
   
   for ( int i=0; i<nw; ++i )
   {
      npizero[i] = 0.;
      errnpi0[i] = 0.;
      Dnpi0[i]   = 0.;
   }

}

Pi0D::~Pi0D()
{

   Clear();
   delete [] npizero;
   delete [] errnpi0;
   delete [] Dnpi0;

}

void Pi0D::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();
   
   int nw = CommonCalc::GetInstance()->GetNW();
   for ( int i=0; i<nw; ++i )
   {
      npizero[i] = 0.;
      errnpi0[i] = 0.;
      Dnpi0[i] = 0.;
   }
   
   return;

}

void Pi0D::AnalyzeEvent( const EventRecord& event )
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

   npizero[ipos] += weight*npi0;
   errnpi0[ipos] += weight*pow((double)npi0,2);
   
   return;

}

void Pi0D::EndOfEventLoopPlots()
{
    
   int nw = CommonCalc::GetInstance()->GetNW();
   double* nv = CommonCalc::GetInstance()->GetNV();   
   double* nv2 = CommonCalc::GetInstance()->GetNV2();
   
   for (int i = 0; i<nw; i++)
   {
       if (nv[i]>0.)
       {
	  npizero[i]   /= nv[i];
	  Dnpi0[i]  = sqrt(errnpi0[i]/nv[i]-pow(npizero[i],2));
	  errnpi0[i]= Dnpi0[i]/sqrt(nv2[i]);
       }
   }
    
   // why did the original authour wanted a graph w/errors if the errors are all zeros ???
   // TGraph* curgraph = new  TGraph(nw-1,npizero,Dnpi0);
   // BTW, why isn't it TGraphErrors if we do have errors calculated ???!!! 
   // 
   // Note: I can't use the ctor: 
   //       TH1F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   //       because the content of npizero array in NOT guaranteed to be in increasing order.
   //       Thus I invent this trick.
   //
   TGraph* gr = new TGraph(nw-1,npizero,Dnpi0);
   fResult = new TH1F( *(gr->GetHistogram()) );
   fResult->SetNameTitle( "Pi0D", "pi0 dispersion vs mult" );
//    fResult = new TH1F( "Pi0D", "pi0 dispersion D- vs mult", nw-1, npizero );
   fResult->SetDirectory(0); // in order to NOT add it to the already open Events TFile (Genie's) 
   
   for ( int i=0; i<nw; ++i )
   {
      fResult->Fill( npizero[i], Dnpi0[i] );
//       fResult->SetBinError( i, errnpi0[i] );
   }
   
   fResult->SetLineColor(kRed);
   fResult->SetMarkerStyle(20);
   fResult->SetMarkerColor(kRed);
   fResult->SetMarkerSize(0.8);
   
   delete gr;
   
   return;
}
