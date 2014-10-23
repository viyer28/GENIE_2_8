//____________________________________________________________________________
/*!

\class    genie::EffectiveSF

\brief    An effective spectral function to match psi' superscaling.
          Implements the NuclearModelI interface.

\ref      http://arxiv.org/abs/1405.0583

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

  // Initializes the momentum probability distribution based on the parameters
  // described in http://arxiv.org/abs/1405.0583 or supplied in EffectiveSF.xml
  TH1D * MakeEffectiveSF(const Target & target) const;

  // Initializes the momentum probability distribution with the given
  // parameters according to the functional form from
  // http://arxiv.org/abs/1405.0583
  TH1D * MakeEffectiveSF(double bs, double bp, double alpha, double beta,
                         double c1, double c2, double c3,
                         const Target & target) const;

  // Returns the binding energy given in http://arxiv.org/abs/1405.0583 or
  // one supplied in EffectiveSF.xml.
  double ReturnBindingEnergy(const Target & target) const;
  double GetTransEnh1p1hMod(const Target& target) const;
  
  // Returns f1p1h, the probability of interaction via the 1p1h process,
  // given in the reference or supplied in EffectiveSF.xml.
  double Returnf1p1h(const Target & target) const;
  void   LoadConfig (void);

  mutable map<string, TH1D *> fProbDistroMap;
  double fPMax;
  double fPCutOff;

  // Map from PDG code to spectral function parameters
  map<int, double> fNucRmvE;
  map<int, double> f1p1hMap;
  map<int, std::vector<double> > fProbDistParams;
  map<int, double> fTransEnh1p1hMods;
  
  // Map from range of A (pair<lowA, highA> inclusive> to spectral
  // function parameters.
  map<pair<int, int>, double> fRangeNucRmvE;
  map<pair<int, int>, double> fRange1p1hMap;
  map<pair<int, int>, std::vector<double> > fRangeProbDistParams;
  map<pair<int, int>, double> fRangeTransEnh1p1hMods;
};

}         // genie namespace
#endif    // _EFFECTIVE_SF_H_

