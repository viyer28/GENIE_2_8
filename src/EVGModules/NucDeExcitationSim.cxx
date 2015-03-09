//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ March 05, 2008 - CA
   This event generation module was added in version 2.3.1. The initial
   implementation handles 16O only.
 @ Sep 15, 2009 - CA
   IsNucleus() is no longer available in GHepParticle. Use pdg::IsIon().
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGModules/NucDeExcitationSim.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim() :
	EventRecordVisitorI("genie::NucDeExcitationSim")
{

}
//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim(string config) :
	EventRecordVisitorI("genie::NucDeExcitationSim", config)
{

}
//___________________________________________________________________________
NucDeExcitationSim::~NucDeExcitationSim()
{

}
//___________________________________________________________________________
void NucDeExcitationSim::ProcessEventRecord(GHepRecord * evrec) const
{
	LOG("NucDeEx", pNOTICE) 
		<< "Simulating nuclear de-excitation gamma rays";

	GHepParticle * nucltgt = evrec->TargetNucleus();
	if (!nucltgt) {
		LOG("NucDeEx", pINFO) 
			<< "No nuclear target found - Won't simulate nuclear de-excitation";
		return;
	}

	if(nucltgt->Z()==8) this->OxygenTargetSim(evrec);
	if(nucltgt->Z()==18) this->ArgonTargetSim(evrec);

	LOG("NucDeEx", pINFO) 
		<< "Done with this event";
}
//__________________________________________________________________________
void NucDeExcitationSim::ArgonTargetSim(GHepRecord * evrec) const
{
	LOG("NucDeEx", pNOTICE)
		<< "Simulating nuclear de-excitation gamma rays for Argon target";

	GHepParticle * hitnuc = evrec->HitNucleon();
	if(!hitnuc) return;

	Interaction * interaction = evrec->Summary();
	double neutrino_energy = interaction->InitState().ProbeE(kRfLab)*1000;


	// const Int_t numev=400;
	const Int_t numev=1;
	TRandom3* r = new TRandom3();

	Double_t vertex[3]; 
	const Int_t num_start_levels=21; 

	Double_t start_level_energy[num_start_levels] = {2.281,2.752,2.937,3.143,3.334,3.569,3.652,3.786,3.861,4.067,4.111,4.267,4.364,4.522,4.655,4.825,5.017,5.080,5.223,5.696,6.006} ;
	Double_t b[num_start_levels] = {0.9,1.5,0.11,0.06,0.04,0.01,0.16,0.26,0.01,0.05,0.11,0.29,3.84,0.31,0.38,0.47,0.36,0.23,0.03,0.11,0.13}; 

	Double_t w=0;//Energy of electron
	Double_t p=0;//Momentum of electron
	Double_t f=0; //Fermi function
	Double_t level_cs[num_start_levels]; 
	Int_t highest_level=0;

	for(Int_t i=0;i<num_start_levels;i++) {
		w= (neutrino_energy-(start_level_energy[i]+1.5))*1000; //Electron energy in keV
		Double_t sigma=0.0; 
		if(neutrino_energy>(start_level_energy[i]+1.5) && w>=511.) {
			highest_level++;
			for(Int_t n=0;n<=i;n++) {
				w= (neutrino_energy-(start_level_energy[n]+1.5)) *1000; //Electron energy for level n to use in cross section
				p= sqrt( pow(w,2) - pow(511.0,2) ); //Electron momentum in keV/c
				f= sqrt(3.0634+(0.6814/(w-1))); //Fermi function approximation
				sigma+= 1.6e-8 *(p*w*f*b[n]); //In cm^2 * 10^-42
			} 
		}
		level_cs[i]=sigma; 
	}

	Double_t tcs=0; 
	for(Int_t i=0;i<num_start_levels;i++) { //Calculating total cross section
		tcs+=level_cs[i]; 
	}

	Double_t start_level_prob[num_start_levels];

	Int_t i=0; //Event counter
	for(Int_t n=0;n<num_start_levels;n++) { //Calculating each level's probability
		if(tcs==0) { 
			i=numev+1;
			break;
		}
		start_level_prob[n] = level_cs[n]/tcs;
	}

	const Int_t numlevels=73; 

	// Branching ratios 

	std::vector<double> br[numlevels];
	// Level to decay to 

	std::vector<int> decayto[numlevels];

	br[0].push_back(1.00);
	decayto[0].push_back(1);

	br[1].push_back(1.00);
	decayto[1].push_back(4);

	br[2].push_back(0.133891);
	br[2].push_back(0.41841);
	br[2].push_back(0.351464);
	br[2].push_back(0.0962343);
	decayto[2].push_back(38);
	decayto[2].push_back(43);
	decayto[2].push_back(59);
	decayto[2].push_back(72);

	br[3].push_back(0.391705);
	br[3].push_back(0.460829);
	br[3].push_back(0.147465);
	decayto[3].push_back(61);
	decayto[3].push_back(65);
	decayto[3].push_back(71);

	br[4].push_back(0.358974);
	br[4].push_back(0.641026);
	decayto[4].push_back(13);
	decayto[4].push_back(58);

	br[5].push_back(0.247808);
	br[5].push_back(0.0468929);
	br[5].push_back(0.213496);
	br[5].push_back(0.11056);
	br[5].push_back(0.381243);
	decayto[5].push_back(40);
	decayto[5].push_back(52);
	decayto[5].push_back(67);
	decayto[5].push_back(71);
	decayto[5].push_back(72);

	br[6].push_back(0.361025);
	br[6].push_back(0.056677);
	br[6].push_back(0.194099);
	br[6].push_back(0.388199);
	decayto[6].push_back(52);
	decayto[6].push_back(55);
	decayto[6].push_back(63);
	decayto[6].push_back(68);

	br[7].push_back(0.0300963);
	br[7].push_back(0.0613965);
	br[7].push_back(0.0152488);
	br[7].push_back(0.0321027);
	br[7].push_back(0.0513644);
	br[7].push_back(0.0682183);
	br[7].push_back(0.073435);
	br[7].push_back(0.0100321);
	br[7].push_back(0.0120385);
	br[7].push_back(0.0922953);
	br[7].push_back(0.152488);
	br[7].push_back(0.401284);
	decayto[7].push_back(17);
	decayto[7].push_back(25);
	decayto[7].push_back(27);
	decayto[7].push_back(32);
	decayto[7].push_back(49);
	decayto[7].push_back(54);
	decayto[7].push_back(56);
	decayto[7].push_back(62);
	decayto[7].push_back(63);
	decayto[7].push_back(67);
	decayto[7].push_back(68);
	decayto[7].push_back(70);

	br[8].push_back(0.0089877);
	br[8].push_back(0.06386);
	br[8].push_back(0.13245);
	br[8].push_back(0.473037);
	br[8].push_back(0.321665);
	decayto[8].push_back(43);
	decayto[8].push_back(56);
	decayto[8].push_back(67);
	decayto[8].push_back(70);
	decayto[8].push_back(71);

	br[9].push_back(0.16129);
	br[9].push_back(0.193548);
	br[9].push_back(0.645161);
	decayto[9].push_back(37);
	decayto[9].push_back(65);
	decayto[9].push_back(72);

	br[10].push_back(0.245826);
	br[10].push_back(0.0816327);
	br[10].push_back(0.20872);
	br[10].push_back(0.463822);
	decayto[10].push_back(40);
	decayto[10].push_back(56);
	decayto[10].push_back(60);
	decayto[10].push_back(71);

	br[11].push_back(0.170984);
	br[11].push_back(0.310881);
	br[11].push_back(0.518135);
	decayto[11].push_back(29);
	decayto[11].push_back(54);
	decayto[11].push_back(66);

	br[12].push_back(0.242424);
	br[12].push_back(0.757576);
	decayto[12].push_back(54);
	decayto[12].push_back(62);

	br[13].push_back(0.159664);
	br[13].push_back(0.840336);
	decayto[13].push_back(48);
	decayto[13].push_back(58);

	br[14].push_back(0.288136);
	br[14].push_back(0.145763);
	br[14].push_back(0.118644);
	br[14].push_back(0.108475);
	br[14].push_back(0.338983);
	decayto[14].push_back(56);
	decayto[14].push_back(66);
	decayto[14].push_back(70);
	decayto[14].push_back(71);
	decayto[14].push_back(72);

	br[15].push_back(0.00869263);
	br[15].push_back(0.087274);
	br[15].push_back(0.0973574);
	br[15].push_back(0.288595);
	br[15].push_back(0.347705);
	br[15].push_back(0.170376);
	decayto[15].push_back(40);
	decayto[15].push_back(64);
	decayto[15].push_back(65);
	decayto[15].push_back(68);
	decayto[15].push_back(70);
	decayto[15].push_back(71);

	br[16].push_back(1);
	decayto[16].push_back(65);

	br[17].push_back(0.341763);
	br[17].push_back(0.0567327);
	br[17].push_back(0.164046);
	br[17].push_back(0.239234);
	br[17].push_back(0.198223);
	decayto[17].push_back(35);
	decayto[17].push_back(39);
	decayto[17].push_back(51);
	decayto[17].push_back(67);
	decayto[17].push_back(69);

	br[18].push_back(0.106406);
	br[18].push_back(0.254103);
	br[18].push_back(0.0465855);
	br[18].push_back(0.529381);
	br[18].push_back(0.0635257);
	decayto[18].push_back(60);
	decayto[18].push_back(61);
	decayto[18].push_back(63);
	decayto[18].push_back(70);
	decayto[18].push_back(72);

	br[19].push_back(0.0905797);
	br[19].push_back(0.228261);
	br[19].push_back(0.181159);
	br[19].push_back(0.137681);
	br[19].push_back(0.362319);
	decayto[19].push_back(43);
	decayto[19].push_back(45);
	decayto[19].push_back(52);
	decayto[19].push_back(70);
	decayto[19].push_back(71);

	br[20].push_back(0.0316414);
	br[20].push_back(0.0415293);
	br[20].push_back(0.0362558);
	br[20].push_back(0.0481213);
	br[20].push_back(0.0903098);
	br[20].push_back(0.0929466);
	br[20].push_back(0.659196);
	decayto[20].push_back(30);
	decayto[20].push_back(32);
	decayto[20].push_back(46);
	decayto[20].push_back(61);
	decayto[20].push_back(64);
	decayto[20].push_back(66);
	decayto[20].push_back(70);

	br[21].push_back(0.310757);
	br[21].push_back(0.398406);
	br[21].push_back(0.290837);
	decayto[21].push_back(64);
	decayto[21].push_back(66);
	decayto[21].push_back(70);

	br[22].push_back(0.152542);
	br[22].push_back(0.847458);
	decayto[22].push_back(67);
	decayto[22].push_back(71);

	br[23].push_back(0.168367);
	br[23].push_back(0.321429);
	br[23].push_back(0.510204);
	decayto[23].push_back(52);
	decayto[23].push_back(70);
	decayto[23].push_back(71);

	br[24].push_back(0.0276437);
	br[24].push_back(0.078982);
	br[24].push_back(0.0245722);
	br[24].push_back(0.162352);
	br[24].push_back(0.184291);
	br[24].push_back(0.438789);
	br[24].push_back(0.0833699);
	decayto[24].push_back(36);
	decayto[24].push_back(53);
	decayto[24].push_back(62);
	decayto[24].push_back(64);
	decayto[24].push_back(70);
	decayto[24].push_back(71);
	decayto[24].push_back(72);

	br[25].push_back(0.0163934);
	br[25].push_back(0.336276);
	br[25].push_back(0.226986);
	br[25].push_back(0.420345);
	decayto[25].push_back(43);
	decayto[25].push_back(67);
	decayto[25].push_back(68);
	decayto[25].push_back(70);

	br[26].push_back(0.184332);
	br[26].push_back(0.0460829);
	br[26].push_back(0.460829);
	br[26].push_back(0.078341);
	br[26].push_back(0.230415);
	decayto[26].push_back(53);
	decayto[26].push_back(54);
	decayto[26].push_back(60);
	decayto[26].push_back(61);
	decayto[26].push_back(71);

	br[27].push_back(0.0147145);
	br[27].push_back(0.0170689);
	br[27].push_back(0.0500294);
	br[27].push_back(0.329606);
	br[27].push_back(0.588582);
	decayto[27].push_back(36);
	decayto[27].push_back(46);
	decayto[27].push_back(56);
	decayto[27].push_back(67);
	decayto[27].push_back(68);

	br[28].push_back(1);
	decayto[28].push_back(70);

	br[29].push_back(0.280702);
	br[29].push_back(0.0935673);
	br[29].push_back(0.0409357);
	br[29].push_back(0.584795);
	decayto[29].push_back(63);
	decayto[29].push_back(66);
	decayto[29].push_back(68);
	decayto[29].push_back(70);

	br[30].push_back(0.393701);
	br[30].push_back(0.283465);
	br[30].push_back(0.188976);
	br[30].push_back(0.133858);
	decayto[30].push_back(61);
	decayto[30].push_back(67);
	decayto[30].push_back(71);
	decayto[30].push_back(72);

	br[31].push_back(0.0477737);
	br[31].push_back(0.194805);
	br[31].push_back(0.0245826);
	br[31].push_back(0.269017);
	br[31].push_back(0.463822);
	decayto[31].push_back(45);
	decayto[31].push_back(60);
	decayto[31].push_back(65);
	decayto[31].push_back(71);
	decayto[31].push_back(72);

	br[32].push_back(0.105769);
	br[32].push_back(0.129808);
	br[32].push_back(0.0528846);
	br[32].push_back(0.480769);
	br[32].push_back(0.230769);
	decayto[32].push_back(46);
	decayto[32].push_back(56);
	decayto[32].push_back(66);
	decayto[32].push_back(70);
	decayto[32].push_back(71);

	br[33].push_back(0.21875);
	br[33].push_back(0.78125);
	decayto[33].push_back(67);
	decayto[33].push_back(71);

	br[34].push_back(0.102011);
	br[34].push_back(0.179598);
	br[34].push_back(0.718391);
	decayto[34].push_back(43);
	decayto[34].push_back(61);
	decayto[34].push_back(66);

	br[35].push_back(0.00826446);
	br[35].push_back(0.393546);
	br[35].push_back(0.334514);
	br[35].push_back(0.263676);
	decayto[35].push_back(64);
	decayto[35].push_back(67);
	decayto[35].push_back(68);
	decayto[35].push_back(70);

	br[36].push_back(0.056338);
	br[36].push_back(0.704225);
	br[36].push_back(0.239437);
	decayto[36].push_back(51);
	decayto[36].push_back(70);
	decayto[36].push_back(71);

	br[37].push_back(0.21875);
	br[37].push_back(0.78125);
	decayto[37].push_back(67);
	decayto[37].push_back(70);

	br[38].push_back(0.181818);
	br[38].push_back(0.757576);
	br[38].push_back(0.0606061);
	decayto[38].push_back(66);
	decayto[38].push_back(71);
	decayto[38].push_back(72);

	br[39].push_back(0.157258);
	br[39].push_back(0.403226);
	br[39].push_back(0.237903);
	br[39].push_back(0.201613);
	decayto[39].push_back(62);
	decayto[39].push_back(70);
	decayto[39].push_back(71);
	decayto[39].push_back(72);

	br[40].push_back(0.0740741);
	br[40].push_back(0.925926);
	decayto[40].push_back(52);
	decayto[40].push_back(72);

	br[41].push_back(0.0535714);
	br[41].push_back(0.35119);
	br[41].push_back(0.595238);
	decayto[41].push_back(67);
	decayto[41].push_back(68);
	decayto[41].push_back(70);

	br[42].push_back(0.00816803);
	br[42].push_back(0.0583431);
	br[42].push_back(0.350058);
	br[42].push_back(0.583431);
	decayto[42].push_back(49);
	decayto[42].push_back(62);
	decayto[42].push_back(71);
	decayto[42].push_back(72);

	br[43].push_back(0.0961538);
	br[43].push_back(0.423077);
	br[43].push_back(0.480769);
	decayto[43].push_back(66);
	decayto[43].push_back(67);
	decayto[43].push_back(68);

	br[44].push_back(0.450549);
	br[44].push_back(0.549451);
	decayto[44].push_back(69);
	decayto[44].push_back(72);

	br[45].push_back(0.469925);
	br[45].push_back(0.0836466);
	br[45].push_back(0.446429);
	decayto[45].push_back(61);
	decayto[45].push_back(65);
	decayto[45].push_back(72);

	br[46].push_back(0.0408163);
	br[46].push_back(0.510204);
	br[46].push_back(0.44898);
	decayto[46].push_back(67);
	decayto[46].push_back(70);
	decayto[46].push_back(71);

	br[47].push_back(1);
	decayto[47].push_back(72);

	br[48].push_back(0.641026);
	br[48].push_back(0.358974);
	decayto[48].push_back(58);
	decayto[48].push_back(69);

	br[49].push_back(0.0458015);
	br[49].push_back(0.954198);
	decayto[49].push_back(66);
	decayto[49].push_back(70);

	br[50].push_back(0.401639);
	br[50].push_back(0.188525);
	br[50].push_back(0.409836);
	decayto[50].push_back(61);
	decayto[50].push_back(69);
	decayto[50].push_back(72);

	br[51].push_back(0.0188679);
	br[51].push_back(0.173585);
	br[51].push_back(0.754717);
	br[51].push_back(0.0528302);
	decayto[51].push_back(61);
	decayto[51].push_back(67);
	decayto[51].push_back(71);
	decayto[51].push_back(72);

	br[52].push_back(0.0100849);
	br[52].push_back(0.00796178);
	br[52].push_back(0.530786);
	br[52].push_back(0.451168);
	decayto[52].push_back(59);
	decayto[52].push_back(68);
	decayto[52].push_back(70);
	decayto[52].push_back(71);

	br[53].push_back(0.0379902);
	br[53].push_back(0.0490196);
	br[53].push_back(0.612745);
	br[53].push_back(0.300245);
	decayto[53].push_back(67);
	decayto[53].push_back(70);
	decayto[53].push_back(71);
	decayto[53].push_back(72);

	br[54].push_back(0.106557);
	br[54].push_back(0.819672);
	br[54].push_back(0.0737705);
	decayto[54].push_back(59);
	decayto[54].push_back(68);
	decayto[54].push_back(70);

	br[55].push_back(0.699301);
	br[55].push_back(0.300699);
	decayto[55].push_back(64);
	decayto[55].push_back(70);

	br[56].push_back(1);
	decayto[56].push_back(71);

	br[57].push_back(1);
	decayto[57].push_back(72);

	br[58].push_back(0.888099);
	br[58].push_back(0.111901);
	decayto[58].push_back(69);
	decayto[58].push_back(72);

	br[59].push_back(0.00647298);
	br[59].push_back(0.752672);
	br[59].push_back(0.165588);
	br[59].push_back(0.0752672);
	decayto[59].push_back(65);
	decayto[59].push_back(70);
	decayto[59].push_back(71);
	decayto[59].push_back(72);

	br[60].push_back(0.0708556);
	br[60].push_back(0.668449);
	br[60].push_back(0.260695);
	decayto[60].push_back(65);
	decayto[60].push_back(71);
	decayto[60].push_back(72);

	br[61].push_back(0.166667);
	br[61].push_back(0.833333);
	decayto[61].push_back(69);
	decayto[61].push_back(72);

	br[62].push_back(0.0898551);
	br[62].push_back(0.57971);
	br[62].push_back(0.330435);
	decayto[62].push_back(67);
	decayto[62].push_back(68);
	decayto[62].push_back(70);

	br[63].push_back(0.813008);
	br[63].push_back(0.186992);
	decayto[63].push_back(71);
	decayto[63].push_back(72);

	br[64].push_back(0.29078);
	br[64].push_back(0.70922);
	decayto[64].push_back(70);
	decayto[64].push_back(71);

	br[65].push_back(0.05);
	br[65].push_back(0.08);
	br[65].push_back(0.5);
	br[65].push_back(0.37);
	decayto[65].push_back(69);
	decayto[65].push_back(70);
	decayto[65].push_back(71);
	decayto[65].push_back(72);

	br[66].push_back(0.398406);
	br[66].push_back(0.310757);
	br[66].push_back(0.290837);
	decayto[66].push_back(70);
	decayto[66].push_back(71);
	decayto[66].push_back(72);

	br[67].push_back(0.819672);
	br[67].push_back(0.180328);
	decayto[67].push_back(70);
	decayto[67].push_back(71);

	br[68].push_back(0.186992);
	br[68].push_back(0.813008);
	decayto[68].push_back(70);
	decayto[68].push_back(71);

	br[69].push_back(1);
	decayto[69].push_back(72);

	br[70].push_back(1);
	decayto[70].push_back(71);

	br[71].push_back(1);
	decayto[71].push_back(72);

	Int_t ground_level=numlevels-1;

	while (i<numev) {	//Event loop 
		Double_t rsl = r->Rndm();
		Double_t tprob=0;
		Int_t chosen_start_level = -1;
		for(Int_t n=0;n<highest_level;n++) { //Picking a starting level
			if(rsl<(start_level_prob[n]+tprob)) {
				chosen_start_level=n;
				break;
			}
			tprob+= start_level_prob[n];
		}	

		vertex[0]=(r->Rndm())*256.35;
		vertex[1]=((r->Rndm())*233)-116.5;
		vertex[2]=(r->Rndm())*1036;

		Int_t particle_c=0;

		Int_t lastlevel= -1;
		Int_t level= -1;

		// Energy levels
		Double_t level_energy[numlevels] = {7.4724,6.2270,5.06347,4.99294,4.8756,4.87255,4.78865,4.744093,4.53706,4.47299,4.41936,4.39588,4.3840,4.3656,4.28052,4.25362,4.21307,4.18003,4.14901,4.11084,4.10446,4.02035,3.92390,3.88792,3.86866,3.840228,3.82143,3.79757,3.76779,3.73848,3.663739,3.62995,3.59924,3.55697,3.48621,3.439144,3.41434,3.39363,3.36803,3.22867,3.15381,3.14644,3.12836,3.109721,3.1002,3.02795,2.98587,2.9508,2.87901,2.80788,2.7874,2.786644,2.75672,2.74691,2.730372,2.625990,2.57593,2.558,2.54277,2.419171,2.397165,2.290493,2.289871,2.26040,2.103668,2.069809,2.047354,1.959068,1.643639,0.891398,0.8001427,0.0298299,0.0};

		//Time delays
		Double_t level_delay[numlevels];
		for(Int_t n=0;n<numlevels;n++) {
			level_delay[n] = 0.0;
		}

		Int_t highesthigher=0; //The highest n for which start_level_energy[chosen_start_level] is higher than level_energy[n]
		Int_t lowestlower=0; //The lowest n for which start_level_energy[chosen_start_level] is lower than level_energy[n]

		for(Int_t n=0;n<numlevels;n++) { //Finding lowestlower and highesthigher
			if(start_level_energy[chosen_start_level]<level_energy[n]) {
				lowestlower=n;
			}
			if(start_level_energy[chosen_start_level]>level_energy[n]) {
				highesthigher=n;
				break;
			}
		}

		
		Double_t llimit = start_level_energy[chosen_start_level]-level_energy[lowestlower];
		Double_t ulimit = start_level_energy[chosen_start_level]-level_energy[highesthigher];


		Double_t abslvl = llimit * ( (llimit<0) * (-1) + (llimit>0));
		Double_t abslvu = ulimit * ( (ulimit<0) * (-1) + (ulimit>0));	

		if(abslvl<abslvu) {
			lastlevel = lowestlower;
			level = lowestlower;
		} //If the chosen start level energy is closest to the lowest level energy that it's lower than than the highest level energy that it's higher than, it starts at the level of the lowest level energy that it's lower than.

		if(abslvu<abslvl) {
			lastlevel = highesthigher;
			level = highesthigher;
		} //If the chosen start level energy is closest to the highest level energy that it's higher than than the lowest level energy that it's lower than, it starts at the level of the highest level energy that it's higher than.

		Double_t ttime = 0;
		Int_t nomoredecay = 0;

		while(level != ground_level){	//Level loop
			Double_t rl = r->Rndm();

			Int_t decaynum = 0;

			tprob=0; //Used this variable above for cross section

			for (unsigned int ilevel=0; ilevel < br[level].size(); ilevel++) { //Decay loop
				if (rl<(br[level][ilevel]+tprob)) {

					// We have a decay

					particle_c++;
					level = decayto[level][decaynum];

					Double_t gamma_energy = level_energy[lastlevel]-level_energy[level];
					Double_t gamma_energy_G = gamma_energy / 1000;

					Double_t time = (-TMath::Log(r->Rndm())/(1/(level_delay[lastlevel]))) + ttime;

					double dt=-1;

					this->AddPhoton(evrec, gamma_energy_G, dt);
					LOG("NucDeEx", pNOTICE)
                                << "Added photon";
					lastlevel = level;
					ttime = time;

					break;

				}

				if((tprob+br[level][ilevel])>1) {
					nomoredecay = 1; //If it doesn't do any more gamma decay
					break;
				}

				decaynum++;
				tprob += br[level][ilevel];

			}	//End of decay loop

			if(nomoredecay == 1) {
				break;
			}

		}	//End of level loop
        i++;
	}
}
//___________________________________________________________________________
void NucDeExcitationSim::OxygenTargetSim(GHepRecord * evrec) const
{
	LOG("NucDeEx", pNOTICE) 
		<< "Simulating nuclear de-excitation gamma rays for Oxygen target";

	//LOG("NucDeEx", pNOTICE) << *evrec;

	GHepParticle * hitnuc = evrec->HitNucleon();
	if(!hitnuc) return;

	bool p_hole = (hitnuc->Pdg() == kPdgProton);
	double dt   = -1;

	RandomGen * rnd = RandomGen::Instance();

	//
	// ****** P-Hole 
	//
	if (p_hole) {
		// 
		// * Define all the data required for simulating deexcitations of p-hole states
		//

		// > probabilities for creating a p-hole in the P1/2, P3/2, S1/2 shells
		double Pp12 = 0.25;              // P1/2 
		double Pp32 = 0.47;              // P3/2 
		double Ps12 = 1. - Pp12 - Pp32;  // S1/2 

		// > excited state energy levels & probabilities for P3/2-shell p-holes
		const int np32 = 3;
		double p32Elv[np32] = { 0.00632, 0.00993, 0.01070 };
		double p32Plv[np32] = { 0.872,   0.064,   0.064   }; 
		// - probabilities for deexcitation modes of P3/2-shell p-hole state '1' 
		double p32Plv1_1gamma  = 0.78;  // prob to decay via 1 gamma
		double p32Plv1_cascade = 0.22;  // prob to decay via gamma cascade

		// > excited state energy levels & probabilities for S1/2-shell p-holes
		const int ns12 = 11;
		double s12Elv[ns12] = { 
			0.00309, 0.00368, 0.00385, 0.00444, 0.00492,
			0.00511, 0.00609, 0.00673, 0.00701, 0.00703, 0.00734 }; 
		double s12Plv[ns12] = { 
			0.0625,  0.1875,  0.075,   0.1375,  0.1375,
			0.0125,  0.0125,  0.075,   0.0563,  0.0563,  0.1874  };
		// - gamma energies and probabilities for S1/2-shell p-hole excited 
		//   states '2','7' and '10' with >1 deexcitation modes
		const int ns12lv2 = 3;
		double s12Elv2[ns12lv2]    = { 0.00309, 0.00369, 0.00385 };
		double s12Plv2[ns12lv2]    = { 0.013,   0.360,   0.625   };
		const int ns12lv7 = 2;
		double s12Elv7[ns12lv7]    = { 0.00609, 0.00673 };
		double s12Plv7[ns12lv7]    = { 0.04,    0.96    };
		const int ns12lv10 = 3;
		double s12Elv10[ns12lv10]  = { 0.00609, 0.00673, 0.00734 };
		double s12Plv10[ns12lv10]  = { 0.050,   0.033,   0.017   };

		// Select one of the P1/2, P3/2 or S1/2
		double rshell = rnd->RndDec().Rndm();
		//
		// >> P1/2 shell
		//
		if(rshell < Pp12) {
			LOG("NucDeEx", pNOTICE) 
				<< "Hit nucleon left a P1/2 shell p-hole. Remnant is at g.s.";
			return;
		} 
		//
		// >> P3/2 shell
		//
		else
			if(rshell < Pp12 + Pp32) {
				LOG("NucDeEx", pNOTICE) 
					<< "Hit nucleon left a P3/2 shell p-hole";
				// Select one of the excited states 
				double rdecmode  = rnd->RndDec().Rndm();        
				double prob_sum  = 0;
				int    sel_state = -1;
				for(int istate=0; istate<np32; istate++) {
					prob_sum += p32Plv[istate];
					if(rdecmode < prob_sum) {
						sel_state = istate;
						break;
					}
				}
				LOG("NucDeEx", pNOTICE) 
					<< "Selected P3/2 excited state = " << sel_state;

				// Decay that excited state
				// >> 6.32 MeV state
				if(sel_state==0) { 
					this->AddPhoton(evrec, p32Elv[0], dt);
				} 
				// >> 9.93 MeV state
				else 
					if(sel_state==1) {    
						double r = rnd->RndDec().Rndm();        
						// >>> emit a single gamma 
						if(r < p32Plv1_1gamma) {
							this->AddPhoton(evrec, p32Elv[1], dt);
						}
						// >>> emit a cascade of gammas 
						else 
							if(r < p32Plv1_1gamma + p32Plv1_cascade) {
								this->AddPhoton(evrec, p32Elv[1],           dt);
								this->AddPhoton(evrec, p32Elv[1]-p32Elv[0], dt);
							}
					}
				// >> 10.7 MeV state 
					else 
						if(sel_state==2) {
							// Above the particle production threshold - need to emit
							// a 0.5 MeV kinetic energy proton.
							// Will neglect that given that it is a very low energy
							// kinetic energy nucleon and the intranuke break-up nucleon
							// cross sections are already tuned.
							return;
						}
			} //p3/2
		//
		// >> S1/2 shell
		//
			else if (rshell < Pp12 + Pp32 + Ps12) {
				LOG("NucDeEx", pNOTICE) 
					<< "Hit nucleon left an S1/2 shell p-hole";
				// Select one of the excited states caused by a S1/2 shell hole
				double rdecmode  = rnd->RndDec().Rndm();  
				double prob_sum  = 0;
				int    sel_state = -1;
				for(int istate=0; istate<ns12; istate++) {
					prob_sum += s12Plv[istate];
					if(rdecmode < prob_sum) {
						sel_state = istate;
						break;
					}
				}
				LOG("NucDeEx", pNOTICE) 
					<< "Selected S1/2 excited state = " << sel_state;

				// Decay that excited state
				bool multiple_decay_modes = 
					(sel_state==2 || sel_state==7 || sel_state==10);
				if(!multiple_decay_modes) {
					this->AddPhoton(evrec, s12Elv[sel_state], dt); 
				} else {
					int ndec = -1;
					double * pdec = 0, * edec = 0;
					switch(sel_state) {
						case(2) : 
							ndec = ns12lv2;  pdec = s12Plv2;  edec = s12Elv2;
							break;
						case(7) : 
							ndec = ns12lv7;  pdec = s12Plv7;  edec = s12Elv7;
							break;
						case(10) : 
							ndec = ns12lv10; pdec = s12Plv10; edec = s12Elv10;
							break;
						default:
							return;
					}
					double r = rnd->RndDec().Rndm();  
					double decmode_prob_sum = 0;
					int sel_decmode = -1;
					for(int idecmode=0; idecmode < ndec; idecmode++) {
						decmode_prob_sum += pdec[idecmode];
						if(r < decmode_prob_sum) {
							sel_decmode = idecmode;
							break;
						}
					}
					if(sel_decmode == -1) return;
					this->AddPhoton(evrec, edec[sel_decmode], dt);  
				}//mult.dec.ch 

			} // s1/2
			else {
			}
	} // p-hole

	//
	// ****** n-hole
	//
	else {
		// 
		// * Define all the data required for simulating deexcitations of n-hole states
		//

		// > probabilities for creating a n-hole in the P1/2, P3/2, S1/2 shells
		double Pp12 = 0.25;  // P1/2 
		double Pp32 = 0.44;  // P3/2 
		double Ps12 = 0.09;  // S1/2 
		//>
		double p32Elv = 0.00618;
		//>
		double s12Elv = 0.00703;
		double s12Plv = 0.222;

		// Select one of the P1/2, P3/2 or S1/2
		double rshell = rnd->RndDec().Rndm();
		//
		// >> P1/2 shell
		//
		if(rshell < Pp12) {
			LOG("NucDeEx", pNOTICE) 
				<< "Hit nucleon left a P1/2 shell n-hole. Remnant is at g.s.";
			return;
		} 
		//
		// >> P3/2 shell
		//
		else
			if(rshell < Pp12 + Pp32) {
				LOG("NucDeEx", pNOTICE) 
					<< "Hit nucleon left a P3/2 shell n-hole";
				this->AddPhoton(evrec, p32Elv, dt); 
			} 
		//
		// >> S1/2 shell
		//
			else
				if(rshell < Pp12 + Pp32 + Ps12) {
					LOG("NucDeEx", pNOTICE) 
						<< "Hit nucleon left a S1/2 shell n-hole";
					// only one of the deexcitation modes involve a (7.03 MeV) photon
					double r = rnd->RndDec().Rndm();
					if(r < s12Plv) this->AddPhoton(evrec, s12Elv,dt);
				}
				else {
				}
	} //n-hole
}
//___________________________________________________________________________
void NucDeExcitationSim::AddPhoton(
		GHepRecord * evrec, double E0, double dt) const
{
	// Add a photon at the event record & recoil the remnant nucleus so as to
	// conserve energy/momenta
	//
	double E = (dt>0) ? this->PhotonEnergySmearing(E0, dt) : E0;

	LOG("NucDeEx", pNOTICE) 
		<< "Adding a " << E/units::MeV << " MeV photon from nucl. deexcitation";

	GHepParticle * target  = evrec->Particle(1);
	GHepParticle * remnant = 0;
	int iremnant = -1;
	for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
		remnant  = evrec->Particle(i);
		iremnant = i;
		if(pdg::IsIon(remnant->Pdg())) break;
	}

	TLorentzVector x4(0,0,0,0);
	TLorentzVector p4 = this->Photon4P(E);
	// GHepParticle gamma(kPdgGamma, kIStStableFinalState,iremnant,-1,-1,-1, p4, x4);
	GHepParticle gamma(kPdgGamma, kIStStableFinalState, 1, -1, -1, -1, p4, x4);
	evrec->AddParticle(gamma);  


	remnant->SetPx     ( remnant->Px() - p4.Px() );
	remnant->SetPy     ( remnant->Py() - p4.Py() );
	remnant->SetPz     ( remnant->Pz() - p4.Pz() );
	remnant->SetEnergy ( remnant->E()  - p4.E()  );
}
//___________________________________________________________________________
double NucDeExcitationSim::PhotonEnergySmearing(double E0, double dt) const
{
	// Returns the smeared energy of the emitted gamma
	// E0 : energy of the excited state (GeV)
	// dt: excited state lifetime (sec)
	// 
	double dE = kPlankConstant / (dt*units::s);

	RandomGen * rnd = RandomGen::Instance();
	double E = rnd->RndDec().Gaus(E0 /*mean*/, dE /*sigma*/);        

	LOG("NucDeEx", pNOTICE) 
		<< "<E> = " << E0 << ", dE = " << dE << " -> E = " << E;

	return E;
}
//___________________________________________________________________________
TLorentzVector NucDeExcitationSim::Photon4P(double E) const
{
	// Generate a photon 4p 

	RandomGen * rnd = RandomGen::Instance();

	double costheta = -1. + 2. * rnd->RndDec().Rndm();
	double sintheta = TMath::Sqrt(TMath::Max(0., 1.-TMath::Power(costheta,2)));
	double phi      = 2*kPi * rnd->RndDec().Rndm();
	double cosphi   = TMath::Cos(phi);
	double sinphi   = TMath::Sin(phi);

	double px = E * sintheta * cosphi;
	double py = E * sintheta * sinphi;
	double pz = E * costheta;

	TLorentzVector p4(px,py,pz,E);
	return p4;
}
//___________________________________________________________________________
