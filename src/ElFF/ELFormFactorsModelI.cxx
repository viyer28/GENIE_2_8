//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include <TMath.h>
#include "ElFF/ELFormFactorsModelI.h"
#include "Algorithm/AlgConfigPool.h"
#include "Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI() :
Algorithm()
{

}
//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
ELFormFactorsModelI::~ELFormFactorsModelI()
{

}
