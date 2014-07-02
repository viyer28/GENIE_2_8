#include "GenPlotsBase.h"

#include <TGraphErrors.h>
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

GenPlotsBase::~GenPlotsBase()
{

   Clear();

}

void GenPlotsBase::Clear()
{

   CommonCalc::GetInstance()->Clear();
   CommonCalc::GetInstance()->ResetForNextEvent();   
   
   if ( fResult ) delete fResult;
   fResult = 0;
         
   return;

}

