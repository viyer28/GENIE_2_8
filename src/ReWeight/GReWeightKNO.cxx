//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 \author   Nathan Mayer <nathan.mayer \at tufts.edu>
          Tufts University

          Anne Norrick <aenorrick \at email.wm.edu
          College of William and Mary

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.
 @ May 17, 2010 - CA
   Code extracted from GReWeightNuXSec and redeveloped in preparation for 
   the Summer 2010 T2K analyses.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightKNO.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightKNO::GReWeightKNO() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightKNO::~GReWeightKNO()
{

}
//_______________________________________________________________________________________
bool GReWeightKNO::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {
   case ( kXSecTwkDial_RvpCC1pi      ) : 
   case ( kXSecTwkDial_RvpCC2pi      ) : 
     
     handle = true;
   break;

   default:
     handle = false;
     break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightKNO::SetSystematic(GSyst_t syst, double twk_dial)
{
   switch(syst) {
   case ( kXSecTwkDial_RvpCC1pi      ) : 
   case ( kXSecTwkDial_RvpCC2pi      ) : 
     
     
     fRTwkDial[syst] = twk_dial;
   break;
   
   default:
     break;
   }
}
//_______________________________________________________________________________________
void GReWeightKNO::Reset(void)
{
  map<GSyst_t,double>::iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    it->second = 0.;
  }

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightKNO::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  map<GSyst_t,double>::iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    GSyst_t syst = it->first;
    double twk_dial = it->second;
    double curr = fRDef[syst] * (1 + twk_dial * fracerr->OneSigmaErr(syst));
    fRCurr[syst] = TMath::Max(0.,curr);
  }
}
//_______________________________________________________________________________________
double GReWeightKNO::CalcWeight(const genie::EventRecord & event) 
{
  
  return 1.;
}
//_______________________________________________________________________________________
double GReWeightKNO::CalcChisq()
{
  double chisq = 0.;
  map<GSyst_t,double>::const_iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    double twk_dial = it->second;
    chisq += TMath::Power(twk_dial, 2.);
  }
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightKNO::Init(void)
{

  // Get the default  parameters 
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * user_config = conf_pool->GlobalParameterList();

  fRDef.insert(map<GSyst_t,double>::value_type( kXSecTwkDial_RvpCC1pi, user_config->GetDouble("DIS-HMultWgt-vp-CC-m2" )));

  map<GSyst_t,double>::const_iterator it = fRDef.begin();
  for( ; it != fRDef.end(); ++it) {
    fRCurr.insert(map<GSyst_t,double>::value_type(it->first, it->second));
  }
}
//_______________________________________________________________________________________
