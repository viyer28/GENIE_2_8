#include "NegChHadD.h"

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

NegChHadD::NegChHadD() 
   : GenPlotsBase()
{

   fName = "NegChHad_D";
   
   int nw = CommonCalc::GetInstance()->GetNW();   
   
   Dneg = new double[nw];
   nneg = new double[nw];
   errnneg = new double[nw];
   

   for ( int i=0; i<nw; ++i )
   {
      Dneg[i] = 0.;
      nneg[i] = 0.;
      errnneg[i] = 0.;
   }

}

NegChHadD::~NegChHadD()
{

   Clear();
   delete [] Dneg;
   delete [] nneg;
   delete [] errnneg;

}

void NegChHadD::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();

   int nw = CommonCalc::GetInstance()->GetNW();   
   for ( int i=0; i<nw; ++i )
   {
      Dneg[i] = 0.;
      nneg[i] = 0.;
      errnneg[i] = 0.; 
   }
   
   return;

}

void NegChHadD::AnalyzeEvent( const EventRecord& event )
{
     
   int ncharged = 0;
   int nnegative = 0;
   double weight = CommonCalc::GetInstance()->GetWeight(); // if it's just 1., why would one need it ???
      
   // copied over from the original code
   // (with a couple of minor modifications)
   //
   GHepParticle * p = 0;
   TIter event_iter(&event);   
   while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
   {
      if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) 
      {
         if (p->Charge()) 
	 {
	    ncharged++;
	    if (p->Charge()<0)
	    {
	       nnegative++;
	    }
	 }
      }     
   }

   int ipos = CommonCalc::GetInstance()->GetCurrentW2Bin();      
   if (ipos==-1) return; //out of the energy range

   nneg[ipos]    += weight*(nnegative);
   errnneg[ipos] += weight*pow((double)nnegative,2);
   
   return;

}

void NegChHadD::EndOfEventLoopPlots()
{
    
   static CommonCalc* service = CommonCalc::GetInstance();
   int nw = service->GetNW();
   double* nv = service->GetNV();
   double* nv2 = service->GetNV2();
   
   for (int i = 0; i<nw; i++)
   {
       if (nv[i]>0.)
       {
          nneg[i] /= nv[i];
	  Dneg[i]   = sqrt(errnneg[i]/nv[i]-pow(nneg[i],2));
	  errnneg[i]= Dneg[i]/sqrt(nv2[i]);
       }
   }
    
   // the original graf was no-errors:
   // D_nneg[imc] = new TGraph(nw-1,nneg,Dneg);  
   // 
   // All in all, for the uniformity it's better to make it a histo with variable bin size
   //
   // However, there's at least one flaw here: there's no guarantee if the array elements
   // will go in the increasing order, which is a requirement for the axis definition.
   // FIXME:: need to invent something to ensure it...
   //
   fResult = new TH1F( "NegChHadD", "negative hadron dispersion D- vs average multiplicity <n->", nw-1, nneg );
   fResult->SetDirectory(0); // in order to NOT add it to the already open Events TFile (Genie's) 
   
   for ( int i=0; i<nw; ++i )
   {
      fResult->Fill( nneg[i], Dneg[i] );
//      fResult->SetBinError( i, errnneg[i] );
   }
   
   fResult->SetLineColor(kRed);
   fResult->SetMarkerStyle(20);
   fResult->SetMarkerColor(kRed);
   fResult->SetMarkerSize(0.8);
   
   return;
}
