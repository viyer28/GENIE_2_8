//____________________________________________________________________________
/*!

\class    genie::EffectiveSF

\brief    An effective spectral function to match psi' superscaling.
          Implements the NuclearModelI interface.

\ref      

\author   Brian Coopersmith, University of Rochester

\created  April 02, 2014

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _EFFECTIVE_SF_H_
#define _EFFECTIVE_SF_H_

#include <map>

#include <TH1D.h>
#include "Nuclear/NuclearModelI.h"

using std::map;

namespace genie {

class EffectiveSF : public NuclearModelI {

public:
  EffectiveSF();
  EffectiveSF(string config);
  virtual ~EffectiveSF();

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double p, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmEffSpectralFunc; 
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);
  
private:
  TH1D * ProbDistro (const Target & t) const;
  TH1D * MakeEffectiveSF(const Target & target) const;
  TH1D * MakeEffectiveSF(double bs, double bp, double alpha, double beta, double c1, double c2, double c3, const Target & target) const;
  double ReturnBindingEnergy(const Target & target) const;
  double Returnf1p1h(const Target & target) const;
  void   LoadConfig (void);
  bool   GetDoubleKeyPDG(const char* valName, double & val,const int pdgc) const;

  mutable map<string, TH1D *> fProbDistroMap;
  double fPMax;
  double fPCutOff;
  double fTransEnh1p1hMod;
  map<int, double> fNucRmvE;
  map<int, double> f1p1hMap;
  map<int, std::vector<double> > fProbDistParams;
};

}         // genie namespace
#endif    // _EFFECTIVE_SF_H_

