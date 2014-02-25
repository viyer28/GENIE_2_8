#include "GenPlotsBase.h"

#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLorentzVector.h>
// #include <TProfile.h>

#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::mc_vs_data;

GenPlotsBase::~GenPlotsBase()
{

   Clear();

}

void GenPlotsBase::Clear()
{

   if ( fResult ) delete fResult;
   fResult = 0;
         
   return;

}

