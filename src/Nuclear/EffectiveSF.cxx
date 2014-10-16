//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Brian Coopersmith, University of Rochester

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/EffectiveSF.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"
#include "Utils/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
EffectiveSF::EffectiveSF() :
NuclearModelI("genie::EffectiveSF")
{

}
//____________________________________________________________________________
EffectiveSF::EffectiveSF(string config) :
NuclearModelI("genie::EffectiveSF", config)
{

}
//____________________________________________________________________________
EffectiveSF::~EffectiveSF()
{
  map<string, TH1D*>::iterator iter = fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.begin(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//____________________________________________________________________________
bool EffectiveSF::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());
  fCurrf1p1h = 0;
  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  //-- set fermi momentum vector
  //
  TH1D * prob = this->ProbDistro(target);
  if(!prob) {
    LOG("EffectiveSF", pNOTICE)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }
  double p = prob->GetRandom();
  LOG("EffectiveSF", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;  

  fCurrMomentum.SetXYZ(px, py, pz);

  //-- set removal energy 
  //
  fCurrRemovalEnergy = this->ReturnBindingEnergy(target);
  fCurrf1p1h = this->Returnf1p1h(target);
  // Since TE increases the QE peak via a 2p2h process, we decrease f1p1h
  // in order to increase the 2p2h interaction to account for this enhancement.
  if(target.A() >= 12) {
    fCurrf1p1h /= fTransEnh1p1hMod;
  }
  return true;
}
//____________________________________________________________________________
double EffectiveSF::Prob(double p, double w, const Target & target) const
{
  if(w < 0) {
     TH1D * prob = this->ProbDistro(target);
     int bin = prob->FindBin(p);
     double y  = prob->GetBinContent(bin);
     double dx = prob->GetBinWidth(bin);
     double p  = y * dx;
     return p;
  }
  return 1;
}
//____________________________________________________________________________
TH1D * EffectiveSF::ProbDistro(const Target & target) const
{
  //-- return stored /if already computed/
  map<string, TH1D*>::iterator it = fProbDistroMap.find(target.AsString());
  if(it != fProbDistroMap.end()) return it->second;

  LOG("EffectiveSF", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("EffectiveSF", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  return this->MakeEffectiveSF(target);
	
}
//____________________________________________________________________________
TH1D * EffectiveSF::MakeEffectiveSF(const Target & target) const
{
  int pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  map<int,vector<double> >::const_iterator it = fProbDistParams.find(pdgc);
  if(it != fProbDistParams.end()) {
    vector<double> v = it->second;
    return this->MakeEffectiveSF(v[0], v[1], v[2], v[3],
                                 v[4], v[5], v[6], target);
  }
  // Default parameters for nuclei described in http://arxiv.org/abs/1405.0583
  if(target.A() == 2)
    return this->MakeEffectiveSF(0.413475, 1.75629, 8.29029, 0.00362127,
                                 0.186987, 6.24155, 0.000208154, target);
  if(target.A() == 3)
    return this->MakeEffectiveSF(3.05583, 0.901861, 10.932, 6.02976,
                                 199.551, 1.92379, 5.26025e-05, target);
  else if(target.A() > 3 && target.A() < 8)
    return this->MakeEffectiveSF(2.14131, 0.775316, 9.72776, 7.56679,
                                 183.428, 5.53273, 0.000589888, target);
  else if(target.A() > 7 && target.A() < 17)
    return this->MakeEffectiveSF(2.11978, 0.736561, 12.941, 10.6195,
                                 197.014, 9.94218, 4.36074e-05, target);
  else if(target.A() > 16 && target.A() < 26)
    return this->MakeEffectiveSF(1.81985, 0.609983, 6.81283, 6.07869,
                                 25.8688, 0.590602, 0.00220846, target);
  else if(target.A() > 25 && target.A() < 39)
    return this->MakeEffectiveSF(1.7394, 0.620821, 7.20039, 6.73047,
                                 21.0049, 0.592854, 0.00121452, target);
  else if(target.A() > 38 && target.A() < 56)
    return this->MakeEffectiveSF(1.66649, 0.614772, 8.54288, 8.61544,
                                 199.992, 6.24716, 0.00269039, target);
  else if(target.A() > 55 && target.A() < 61)
    return this->MakeEffectiveSF(1.79235, 0.597282, 7.09555, 6.25598,
                                 18.3744, 0.505288, 0.00140556, target);
  else if(target.A() > 60)
    return this->MakeEffectiveSF(1.52067, 0.585172, 11.2361, 13.3336,
                                 174.424, 5.28673, 9.27737e-05, target);
  return NULL;
}
//____________________________________________________________________________
TH1D * EffectiveSF::MakeEffectiveSF(double bs, double bp, double alpha,
                                    double beta, double c1, double c2,
                                    double c3, const Target & target) const
{
  //-- create the probability distribution
  int npbins = (int) (1000 * fPMax);
  
  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);
  prob->SetDirectory(0);

  double dp = fPMax / (npbins-1);
  for(int i = 0; i < npbins; i++) {
    double p  = i * dp;
    double y = p / 0.197;
    double as = c1 * exp(-pow(bs*y,2));
    double ap = c2 * pow(bp * y, 2) * exp(-pow(bp * y, 2));
    double at = c3 * pow(y, beta) * exp(-alpha * (y - 2));
    double rr = (3.14159265 / 4) * (as + ap + at) * pow(y, 2) / 0.197;
    double dP_dp = rr / 1.01691371;
    if(p>fPCutOff)
      dP_dp = 0;
    assert(dP_dp >= 0);
    // calculate probability density : dProbability/dp
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("EffectiveSF", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
    prob->Fill(p, dP_dp);
 }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  //-- store
  fProbDistroMap.insert(
      map<string, TH1D*>::value_type(target.AsString(),prob));
  return prob; 
}
//____________________________________________________________________________
double EffectiveSF::ReturnBindingEnergy(const Target & target) const
{
  int pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  map<int,double>::const_iterator it = fNucRmvE.find(pdgc);
  if(it != fNucRmvE.end()) {
    return it->second;
  }
  
  if(target.A() == 2)
    return 0.000126547;
  if(target.A() == 3)
    return 0.00534258;
  else if(target.A() > 3 && target.A() < 8)
    return 0.0140444;
  else if(target.A() > 7 && target.A() < 17)
    return 0.0124744;
  else if(target.A() > 16 && target.A() < 26)
    return 0.0166338;
  else if(target.A() > 25 && target.A() < 39)
    return 0.0124529;
  else if(target.A() > 38 && target.A() < 56)
    return 0.0205751;
  else if(target.A() > 55 && target.A() < 61)
    return 0.0151136;
  else if(target.A() > 60)
    return 0.0187988;
  return 0;
}
//____________________________________________________________________________
double EffectiveSF::Returnf1p1h(const Target & target) const
{
  int pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  map<int,double>::const_iterator it = f1p1hMap.find(pdgc);
  if(it != f1p1hMap.end()) {
    return it->second;
  }
  
  if(target.A() == 2)
    return 0;
  if(target.A() == 3)
    return 0.312043;
  else if(target.A() > 3 && target.A() < 8)
    return 0.791142;
  else if(target.A() > 7 && target.A() < 17)
    return 0.807758;
  else if(target.A() > 16 && target.A() < 26)
    return 0.764847;
  else if(target.A() > 25 && target.A() < 39)
   return 0.773977;
  else if(target.A() > 38 && target.A() < 56)
    return 0.808667;
  else if(target.A() > 55 && target.A() < 61)
    return 0.822268;
  else if(target.A() > 60)
    return 0.89556;
  return 0;
}
//____________________________________________________________________________
void EffectiveSF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void EffectiveSF::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void EffectiveSF::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * gc = confp->GlobalParameterList();
  fTransEnh1p1hMod = 1.0;
  if(gc->GetBoolDef("UseTransverseEnhancement", false)) {
    fTransEnh1p1hMod = gc->GetDoubleDef("TransEnhf1p1hMod", 1.18);
  }
  fPMax    = fConfig->GetDoubleDef ("MomentumMax", 1.0);

  fPCutOff = fConfig->GetDoubleDef ("MomentumCutOff", 0.65);

  assert(fPMax > 0 && fPCutOff > 0 && fPCutOff <= fPMax);
  
  for(int Z = 1; Z < 140; Z++) {
    for(int A = Z; A < 3*Z; A++) {
      int pdgc = pdg::IonPdgCode(A, Z);
      //get custom specified effective binding energy
      double eb, f;
      if(GetDoubleKeyPDG("BindingEnergy", eb, pdgc)) {
        eb = TMath::Max(eb, 0.);
        LOG("EffectiveSF", pINFO)
          << "Nucleus: " << pdgc << " -> using Eb =  " << eb << " GeV";
        fNucRmvE.insert(map<int,double>::value_type(pdgc, eb));
      }
      //get custom specified f1p1h
      if(GetDoubleKeyPDG("f1p1h", f, pdgc)) {
        if(f < 0 || f > 1){
          f = 1;
        }
        LOG("EffectiveSF", pINFO)
          << "Nucleus: " << pdgc << " -> using f1p1h =  " << f;
        f1p1hMap.insert(map<int,double>::value_type(pdgc, f));
      }
      double bs, bp, alpha, beta, c1, c2, c3;
      if(GetDoubleKeyPDG("bs", bs, pdgc) && GetDoubleKeyPDG("bp", bp, pdgc) && 
        GetDoubleKeyPDG("alpha", alpha, pdgc) && GetDoubleKeyPDG("beta", beta, pdgc) && 
        GetDoubleKeyPDG("c1", c1, pdgc) && GetDoubleKeyPDG("c2", c2, pdgc) && 
        GetDoubleKeyPDG("c3", c3, pdgc)) {
        vector<double> pars = vector<double>();
        pars.push_back(bs);
        pars.push_back(bp);
        pars.push_back(alpha);
        pars.push_back(beta);
        pars.push_back(c1);
        pars.push_back(c2);
        pars.push_back(c3);
        LOG("EffectiveSF", pINFO)
          << "Nucleus: " << pdgc << " -> using bs =  " << bs << "; bp = "<< bp 
          << "; alpha = " << alpha << "; beta = "<<beta<<"; c1 = "<<c1<<"; c2 = "<<c2
          << "; c3 = " << c3;
        fProbDistParams.insert(map<int,vector<double> >::value_type(pdgc, pars));
      }
    }
  }
}
//____________________________________________________________________________
bool EffectiveSF::GetDoubleKeyPDG(const char* valName, double & val,const int pdgc) const
{
	ostringstream s;
	s<<valName<<"@Pdg="<<pdgc;
	RgKey key = s.str();
	if(!this->GetConfig().Exists(key)) {
	  return false;
	}
  val = fConfig->GetDoubleDef(key,0);
  return true;
}
