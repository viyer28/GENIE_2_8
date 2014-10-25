//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 31, 2009 - CA
   Was first added in v2.5.1.
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "ElFF/TransverseEnhancementFFModel.h"
#include "Interaction/Interaction.h"
#include "Utils/ConfigIsotopeMapUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::config;

//____________________________________________________________________________
TransverseEnhancementFFModel::TransverseEnhancementFFModel() :
ELFormFactorsModelI("genie::TransverseEnhancementFFModel")
{

}
//____________________________________________________________________________
TransverseEnhancementFFModel::TransverseEnhancementFFModel(string config) :
ELFormFactorsModelI("genie::TransverseEnhancementFFModel", config)
{

}
//____________________________________________________________________________
TransverseEnhancementFFModel::~TransverseEnhancementFFModel()
{

}
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gep(const Interaction * interaction) const
{
  return fElFormFactorsBase->Gep(interaction);
}
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gmp(const Interaction * interaction) const
{
  return GetTransEnhMagFF(fElFormFactorsBase->Gmp(interaction), interaction);
}
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gen(const Interaction * interaction) const
{
  return fElFormFactorsBase->Gen(interaction);
}
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gmn(const Interaction * interaction) const
{
  return GetTransEnhMagFF(fElFormFactorsBase->Gmn(interaction), interaction);
}
//____________________________________________________________________________
double TransverseEnhancementFFModel::GetTransEnhMagFF(
    double magFF, const Interaction * interaction) const
{
  const Target& target = interaction->InitState().Tgt();
  double transEnhA, transEnhB;
  GetTransEnhParams(target, &transEnhA, &transEnhB);
  if (transEnhA == 0) {
    return magFF;
  }
  double Q2 = interaction->Kine().Q2();
  double rt = 1 + transEnhA * Q2 * TMath::Exp(-Q2 / transEnhB);
  return TMath::Sqrt(rt)*magFF;
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::GetTransEnhParams(
    const Target& target, double* teA, double* teB) const {
  if (!GetValueFromNuclearMaps(target, fNucMagFF_RT_A,
                               fRangeMagFF_RT_A, teA) ||
      !GetValueFromNuclearMaps(target, fNucMagFF_RT_B,
                               fRangeMagFF_RT_B, teB) ||
      *teB == 0) {
    *teA = 0;
    *teB = 1;
  }
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  LoadAllIsotopesForKey("MagFF_RT_A", "TransverseEnhancementFFModel", fConfig,
                        &fNucMagFF_RT_A);
  LoadAllNucARangesForKey("MagFF_RT_A", "TransverseEnhancementFFModel", fConfig,
                          &fRangeMagFF_RT_A);
  LoadAllIsotopesForKey("MagFF_RT_B", "TransverseEnhancementFFModel", fConfig,
                        &fNucMagFF_RT_B);
  LoadAllNucARangesForKey("MagFF_RT_B", "TransverseEnhancementFFModel", fConfig,
                          &fRangeMagFF_RT_B);
}
