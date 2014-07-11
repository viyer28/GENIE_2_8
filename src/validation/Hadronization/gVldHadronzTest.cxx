//____________________________________________________________________________
/*!

\program gvld_hadronz_test

\brief   Hadronization validation program
         (comparing GENIE models with neutrino bubble chamber data).

\syntax  gvld_hadronz_test -g genie_inputs.xml [-f fmt]

         Options:

          [] Denotes an optional argument.
          
          -f Output plot format (either `ps', `eps', `gif' or `root').
             If `ps' is selected, all plots are saved in a single document.
             If `eps', `gif' or `root' is selected, then each canvas is 
             stored separately. Default: `ps'.
          -g An input XML file for specifying a GHEP event list for each
             model to be considered in the hadronization benchmark test.
             For info on the XML file format see the GSimFiles class documentation.

\author  Original: Tingjun Yang, Hugh Gallagher, Pauli Kehayias, Costas Andreopoulos
         Re-engineered: Julia Yarba (FNAL), Gabriel Perdue (FNAL)

\created March 1, 2009

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>
#include <string>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"
#include "Utils/CmdLnArgParser.h"

#include "RunConfig.h"
#include "Analyzer.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::mc_vs_data;

void PrintSyntax           (void); // leftover from the original code

//____________________________________________________________________________
int main(int argc, char ** argv)
{

  RunConfig* run      = new RunConfig( argc, argv );
  Analyzer*   analyzer = new Analyzer();
  analyzer->SetOutputFormat( run->GetOutputFormat() );
  analyzer->SetExpDataPtr( run->GetExpData() );
    
  do 
  {
  
     run->Next();
     if ( run->IsDone() ) break;
          
     analyzer->Analyze( run->GetCurrentModelName(), run->GetCurrentGSampleName() );
  
  } while ( !run->IsDone() );

  analyzer->DrawResults( run->GetNModels() ); 

  LOG("gvldtest", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("vldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << " gvld_hadronz_test -g genie_inputs.xml -d exp_data_input.xml [-f format] \n";
}
//____________________________________________________________________________
