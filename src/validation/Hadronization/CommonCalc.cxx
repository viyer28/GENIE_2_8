#include "CommonCalc.h"

#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

#include <string>

using namespace genie;
using namespace genie::mc_vs_data;

CommonCalc* CommonCalc::fInstance = 0;

CommonCalc::CommonCalc()
{
   
      fNW = 14;
      W2lo = new double[fNW];
      W2hi = new double[fNW];
// aW2 = new double[fNW];
// nch = new double[fNW];
      nv = new double[fNW];
      nv2 = new double[fNW];
   
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

   for ( int i=0; i<fNW; ++i )
   {
      nv[i] = 0.;
      nv2[i] = 0.;
   }
   
   weight = 1.;
   
   fCurrentW2 = 0.;
   fCurrentW2Bin = -1;
   fEventProcessed = false;

}

CommonCalc::~CommonCalc()
{

   delete [] W2lo;
   delete [] W2hi;
   delete [] nv;
   delete [] nv2;
   delete fInstance;
   fInstance = 0;

}

CommonCalc* CommonCalc::GetInstance()
{

   if ( fInstance == 0 )
   {
      fInstance = new CommonCalc();
   }
   
   return fInstance;

}

void CommonCalc::AnalyzeEvent( const EventRecord& event )
{

   if ( fEventProcessed ) return; // avoid multiple recalculations
                                  // - do it only once per event !

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
      TLorentzVector q = k1-k2;
      double v = (q*p1)/M; // v (E transfer in hit nucleon rest frame)
      double Q2 = -1 * q.M2(); // momemtum transfer
// double y = v*M/(k1*p1); // Inelasticity, y = q*P1/k1*P1
      fCurrentW2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
// double W = TMath::Sqrt(W2);
      //double Phs = sqrt(pow(p2.Px(),2)+pow(p2.Py(),2)+pow(p2.Pz(),2));
      v+=M; //measured v
      
   fCurrentW2Bin = -1; //W2 bin
   for (int idx = 0; idx<fNW; idx++)
   {
      if (fCurrentW2>=W2lo[idx]&&fCurrentW2<W2hi[idx])
      {
fCurrentW2Bin = idx;
      }
   }
   
   fEventProcessed = true;
   
   if ( fCurrentW2Bin == - 1 ) return;
   
   nv[fCurrentW2Bin] += weight;
   nv2[fCurrentW2Bin] += weight*weight;

   return;

}
