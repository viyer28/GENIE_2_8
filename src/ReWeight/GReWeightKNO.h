//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightKNO

\brief    Reweighting KNO Hadronization level.

\author   Nathan Mayer <nathan.mayer \at tufts.edu>
          Tufts University

          Anne Norrick <aenorrick \at email.wm.edu
          Colledge of William and Mary

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_KNO
#define _G_REWEIGHT_KNO

#include <map>

#include "ReWeight/GReWeightI.h"

using std::map;

namespace genie {
namespace rew   {

 class GReWeightKNO : public GReWeightI 
 {
 public:
   GReWeightKNO();
  ~GReWeightKNO();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);


 private:

   void Init (void);

   map<GSyst_t, double> fRTwkDial;
   map<GSyst_t, double> fRDef;
   map<GSyst_t, double> fRCurr;
 };

} // rew   namespace
} // genie namespace

#endif

