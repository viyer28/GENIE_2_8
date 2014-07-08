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
//____________________________________________________________________________
void ELFormFactorsModelI::ConfigTransEnh()
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * gc = confp->GlobalParameterList();
  fUseTransEnh= gc->GetBoolDef("UseTransverseEnhancement", false);
  fTransEnhA= gc->GetDoubleDef("TransEnhMagneticFF_RT_A", 0.0);
  fTransEnhB= gc->GetDoubleDef("TransEnhMagneticFF_RT_B", 0.0);
}
//____________________________________________________________________________
double ELFormFactorsModelI::GetTransEnhMagFF(double magFF, const Interaction * interaction) const
{
  double A = interaction->InitState().Tgt().A();
  if(A < 12 || !fUseTransEnh)
    return magFF;
  double q2 = interaction->Kine().Q2();
  double rt=1+fTransEnhA*q2*TMath::Exp(-q2/fTransEnhB);
  return TMath::Sqrt(rt)*magFF;
}
