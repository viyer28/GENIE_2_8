#include "GenPlotsChHadW2.h"

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

using namespace genie;
using namespace genie::mc_vs_data;

GenPlotsChHadW2::GenPlotsChHadW2() 
   : GenPlotsBase()
{

   fName = "ChHad_W2";
   
   nw = 14;
   W2lo = new double[nw];
   W2hi = new double[nw];
   aW2  = new double[nw];
   nch  = new double[nw];
   nv   = new double[nw];
   
   W2lo[0] = 1.;
   W2lo[1] = 2.;
   W2lo[2] = 3.;
   W2lo[3] = 4.;
   W2lo[4] = 6.;
   W2lo[5] = 8.;
   W2lo[6] = 12.;
   W2lo[7] = 16.;
   W2lo[8] = 23.;
   W2lo[9] = 32.;
   W2lo[10] = 45.;
   W2lo[11] = 63.;
   W2lo[12] = 90.;
   W2lo[13] = 125.;
   
   W2hi[0] = 2.,
   W2hi[1] = 3.;
   W2hi[2] = 4.;
   W2hi[3] = 6.;
   W2hi[4] = 8.;
   W2hi[5] = 12.;
   W2hi[6] = 16.;
   W2hi[7] = 23.;
   W2hi[8] = 32.;
   W2hi[9] = 45.;
   W2hi[10] = 63.;
   W2hi[11] = 90.;
   W2hi[12] = 125.;
   W2hi[13] = 225.;

   for ( int i=0; i<nw; ++i )
   {
      aW2[i] = 0.;
      nch[i] = 0.;
      nv[i]  = 0.; 
   }

}

GenPlotsChHadW2::~GenPlotsChHadW2()
{

   Clear();
   delete [] W2lo;
   delete [] W2hi;
   delete [] aW2;
   delete [] nch;
   delete [] nv;

}

void GenPlotsChHadW2::Clear()
{

   GenPlotsBase::Clear();
   
//   fGenResults.clear();
   
   for ( int i=0; i<nw; ++i )
   {
      aW2[i] = 0.;
      nch[i] = 0.;
      nv[i]  = 0.; 
   }
   
   return;

}

void GenPlotsChHadW2::AnalyzeEvent( const EventRecord& event )
{

      // substantially copied over from the original code
      // 
      // input particles and primary (???) lepton
      //

      GHepParticle* probe = event.Probe();
      assert(probe);
      GHepParticle* target = event.Particle(1); 
      assert(target);
      GHepParticle * fsl = event.FinalStatePrimaryLepton();
      assert(fsl);
      //GHepParticle * fsh = event.FinalStateHadronicSystem();
      GHepParticle * fsh = event.Particle(3);
      assert(fsh);
      GHepParticle * hitnucl = event.HitNucleon();
      if (!hitnucl) return; 

      double M = hitnucl->Mass(); //mass of struck nucleon

      TLorentzVector k1 = *(probe->P4());
      TLorentzVector k2 = *(fsl->P4());
      TLorentzVector p1 = *(hitnucl->P4());
      TLorentzVector p2 = *(fsh->P4());
      TLorentzVector q  = k1-k2; 
      double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
      double Q2 = -1 * q.M2();      // momemtum transfer
//      double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
      double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
//      double W  = TMath::Sqrt(W2);
      //double Phs = sqrt(pow(p2.Px(),2)+pow(p2.Py(),2)+pow(p2.Pz(),2));
      v+=M;  //measured v
      
   int ncharged = 0;
   double weight = 1.; // if it's just 1., why would one need it ???
      
   GHepParticle * p = 0;
   TIter event_iter(&event);   
   while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
   {
      if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) 
      {
         if (p->Charge()) ncharged++;
      }     
   }

   int ipos = -1; //W2 bin
   for (int idx = 0; idx<nw; idx++)
   {
      if (W2>=W2lo[idx]&&W2<W2hi[idx])
      {
	  ipos = idx;
      }
   }
   
   if (ipos==-1) return; //out of the energy range

   aW2[ipos] += weight*W2;
   nv[ipos]  += weight;
   nch[ipos] += weight*(ncharged);
   
   return;

}

void GenPlotsChHadW2::EndOfEventPlots()
{
    
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
   fResult = new TH1F( "ChHadE2", "ch.hadrons mult vs W2", nw-1, aW2 );
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
