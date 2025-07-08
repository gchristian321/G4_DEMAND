#include <fstream>
#include <cassert>
#include <map>
#include <G4SystemOfUnits.hh>
#include "g4gen/PhaseSpace.hh"
#include "g4gen/NuclearMasses.hh"
#include "g4gen/Rng.hh"
#include "g4gen/Error.hh"
#include "g4gen/Particle.hh"
#include "g4gen/NeutronDecay.hh"

namespace {

const G4double kNeutronMass =
	g4gen::NuclearMasses::GetNuclearMass(0, 1)*MeV;

g4gen::RngUniform kRngUniform(0,1); // RNG uniform 0->1
g4gen::RngUniform kRngUniform_neg1_pos1(-1,1); // RNG uniform -1->1
g4gen::RngUniform kRngUniform_0_2pi(0, 2*CLHEP::pi); // RNG uniform 0->2pi

template<typename T> T POW2(const T& t) { return t*t; }
} //namespace


/////////////////////////////////////////////////////////////////
// g4gen::NeutronDecayIntermediate
//

g4gen::NeutronDecayIntermediate::NeutronDecayIntermediate(G4int number_of_neutrons_emitted):
	mNumberOfNeutronsEmitted(number_of_neutrons_emitted),
	mFinal(number_of_neutrons_emitted + 2),
	fVerb(2),
	mFinalFragMass(0),
	mInitial(0)
{ }

g4gen::NeutronDecayIntermediate::~NeutronDecayIntermediate()
{ }

void g4gen::NeutronDecayIntermediate::SetInputParticle(const g4gen::Particle* p)
{
	mInitial = p;
	mFinalFragMass = 
		g4gen::NuclearMasses::GetNuclearMass(mInitial->Z(), mInitial->A() - mNumberOfNeutronsEmitted)*MeV;
}

void g4gen::NeutronDecayIntermediate::SetFinal(G4int indx, const G4LorentzVector& v)
{
	if(size_t(indx) < mFinal.size())	{
		mFinal[indx].set(v.x(), v.y(), v.z(), v.t());
	}
	else {
		if(GetVerboseLevel() > 0) {
			G4GERR << "Invalid index to g4gen::NeutronDecayIntermediate::SetFinal:: " << indx << G4endl;
		}
		exit(1);
	}
}

const G4LorentzVector& g4gen::NeutronDecayIntermediate::GetFinal(G4int indx) const
{
	if(size_t(indx) < mFinal.size())	{
		return mFinal[indx]; 
	}
	else {
		if(GetVerboseLevel() > 0) {
			G4GERR << "Invalid index to g4gen::NeutronDecayIntermediate::GetFinal:: " << indx << G4endl;
		}
		exit(1);
	}
}

void g4gen::NeutronDecayIntermediate::SetDecayParam(const G4String& par, G4double val)
{
	mParams[par] = val;
}

G4double g4gen::NeutronDecayIntermediate::GetDecayParam(const G4String& par)
{
	std::map<G4String, G4double>::iterator it = mParams.find(par);
	if(it != mParams.end()) { return it->second; }

	if(GetVerboseLevel() > 0) {
		G4GERR << "g4gen::NeutronDecayIntermediate::GetDecayParam:: Invalid Parameter "
					 << par << G4endl;
	}
	exit(1);
}

std::unique_ptr<g4gen::Rng> g4gen::NeutronDecayBreitWigner::CreateRngEx()
{
	const G4double ex = GetDecayParam("ex"); // excitation energy
	const G4double width = GetDecayParam("width"); // BW width
	g4gen::Rng* rng = 0;
	if(width < 1e-8) { // treat width < 1e-8 as constant
		rng = new g4gen::RngConstant(ex);
	} else {
		rng = new g4gen::RngBreitWigner(ex, width);
	}
	mRngEx = rng;
	return std::unique_ptr<g4gen::Rng>(rng);
}




/////////////////////////////////////////////////////////////////
// g4gen::OneNeutronDecay
//

g4gen::OneNeutronDecay::OneNeutronDecay():
	g4gen::NeutronDecayBreitWigner(1)
{ }

g4gen::OneNeutronDecay::~OneNeutronDecay()
{ }


G4bool g4gen::OneNeutronDecay::Generate()
{
	assert (mInitial);
	assert (mRngEx);

	G4LorentzVector lvF, lvN;
	g4gen::NeutronEvaporation evap(mInitial->MplusEx(), mFinalFragMass, kNeutronMass);
	evap(&lvF, &lvN);

	G4ThreeVector t3boost = mInitial->Momentum().boostVector();
	lvN.boost(t3boost); // neutron 1
	lvF.boost(t3boost); // final fragment

	SetFinal(0, mInitial->Momentum());
	SetFinal(1, lvF);
	SetFinal(2, lvN);

	return true;
}


/////////////////////////////////////////////////////////////////
// g4gen::TwoNeutronDecayPhaseSpace
//

g4gen::TwoNeutronDecayPhaseSpace::TwoNeutronDecayPhaseSpace(G4bool fsi):
	g4gen::NeutronDecayBreitWigner(2),
	fFSI(fsi)
{
	SetDecayParam("r0", 2.4);
}

g4gen::TwoNeutronDecayPhaseSpace::~TwoNeutronDecayPhaseSpace()
{ }


namespace {
//Function from Marques for FSI three body
//email from MT to JKS on 3/18/2014 (later forwarded by
//JKS to GAC June 2017)
double Cnn0(double X, double r0)
{
	double hc,f0,d0;
	hc = 197.3;   // hc [MeVfm]
	f0 =  18.5;   // nn scatt. length [fm]
	d0 =   2.8;   // effective range  [fm]
	double    a,b,Ref,Imf,f_2,B0,Bi,pi,F1,F2,CNN0;

	if( r0 == 0. )
	{
		CNN0 = 1.;
		return CNN0;
	}

	pi = acos(-1.);
	a = 1./f0+d0*pow(X/hc,2)/2.; //
	b = -X/hc;                 // f = 1/(a+ib)
	Ref = a/(a*a+b*b);       //
	Imf = -b/(a*a+b*b);      //
	f_2 = Ref*Ref+Imf*Imf;
	B0 = -.5*exp(-4.*pow(r0*X/hc,2));    // t0 = 0
	if( X == 0. )
	{
		F1 = 1.;
		F2 = 0.;
	}
	else
	{
		F1 = (erf(2.*r0*X/hc)*sqrt(pi)/2.)*exp(-4.*pow(r0*X/hc,2))/(2.*r0*X/hc);
		F2 = (1.-exp(-4.*pow(r0*X/hc,2)))/(2.*r0*X/hc);
	}
	Bi = .25*f_2/pow(r0,2)*(1.-d0/(2.*sqrt(pi)*r0))+Ref*F1/(sqrt(pi)*r0)-Imf*F2/(2.*r0);

	CNN0 = 1.+B0+Bi;

	return CNN0;
} }

G4bool g4gen::TwoNeutronDecayPhaseSpace::Generate()
{
	assert (mInitial);
	assert (mRngEx);

	G4double mOut[3] = { mFinalFragMass, kNeutronMass, kNeutronMass };

	g4gen::PhaseSpace gen;
	G4bool possible = gen.SetDecay(mInitial->Momentum(), 3, mOut);
	if(!possible) {
		if(GetVerboseLevel() > 1) {
			G4GERR << "g4gen::TwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
			G4GERR << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial->MplusEx() << ", "
						 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		}
		return false;
	}

	if(fFSI) { // Use Final State Interaction
		/**
		 * \note The 'r0' parameter is the "source size" from
		 * PLB 476, 219 (2000). It is set in the 'reaction.dat' 
		 * input file. The PLB paper gives respective source sizes for 
		 * 6He, 11Li, 14Be of r0 = 2.4, 2.7, and 2.2 fm.
		 */
		const G4double r0 = GetDecayParam("r0");
		G4double CnnM = Cnn0(0, r0);
		while(1) {
			G4double rel_weight = gen.Generate() / gen.GetWtMax();
			G4double ran = kRngUniform.Generate();

			G4LorentzVector ptemp1 = *gen.GetDecay(0);
			G4LorentzVector ptemp2 = *gen.GetDecay(1);
			G4LorentzVector ptemp3 = *gen.GetDecay(2);
			G4LorentzVector n2sys = ptemp2 + ptemp3;
			G4ThreeVector n2vel = n2sys.boostVector();
			ptemp2.boost(-n2vel);
			ptemp3.boost(-n2vel);

			// n.b. already in MeV here, but JKS code was ptemp2.P()*1000
			G4double Cnn = Cnn0(ptemp2.v().mag(), r0);
			G4double ran2 = g4gen::RngUniform(0, CnnM).Generate();

			if( (ran < rel_weight) && (ran2 < Cnn) ) { break; }
		}
	}	else { // NO FSI
		G4double relwt, ran;
		do {
			relwt = gen.Generate() / gen.GetWtMax();
			ran = kRngUniform.Generate();
		} while(relwt < ran);
	}

	SetFinal(0, mInitial->Momentum()); // recoil
	SetFinal(1, *(gen.GetDecay(0)));  // fragment
	SetFinal(2, *(gen.GetDecay(1)));  // neutron 1
	SetFinal(3, *(gen.GetDecay(2)));  // neutron 1

	return true;
}


#if 0
/////////////////////////////////////////////////////////////////
// g4gen::TwoNeutronDecaySequential
//

g4gen::TwoNeutronDecaySequential::TwoNeutronDecaySequential():
	g4gen::NeutronDecayIntermediate(2),
	mIntermediateFragMass(0)
{ }

g4gen::TwoNeutronDecaySequential::~TwoNeutronDecaySequential()
{ }

void g4gen::TwoNeutronDecaySequential::SetInputParticle(const g4gen::Particle* p)
{
	g4gen::NeutronDecayIntermediate::SetInputParticle(p);
	mIntermediateFragMass =
		g4gen::NuclearMasses::GetNuclearMass(mInitial->Z(), mInitial->A() - 1)*MeV;
}

std::unique_ptr<g4gen::Rng> g4gen::TwoNeutronDecaySequential::CreateRngEx()
{
	const G4double ex = GetDecayParam("ex"); // excitation energy
	const G4double width = GetDecayParam("width"); // BW width
	const G4int Z_i = GetDecayParam("Z_i");
	const G4int A_i = GetDecayParam("A_i");
	const G4double Mass_i = g4gen::NuclearMasses::GetNuclearMass(Z_i, A_i - 1);
	const G4double Mass_f = g4gen::NuclearMasses::GetNuclearMass(Z_i, A_i - 2);
	const G4double Ev  = Mass_i - Mass_f - kNeutronMass;
	const G4double sI  = GetDecayParam("S_i");
	const G4double sF  = GetDecayParam("S_f");
	const G4double ell = GetDecayParam("L");
	g4gen::RngVolyaSequentialEx* rng = 
		new g4gen::RngVolyaSequentialEx(ex,Ev,sI,sF,ell,width,A_i);
	mRngEx = rng;
	return std::unique_ptr<g4gen::Rng>(rng);
}


G4bool g4gen::TwoNeutronDecaySequential::Generate()
{
/** Code taken from st_reaction.cc out of st_mona simulation in
 *  use by the MoNA collaboration.
 *
 *  \note Using 'Volya' sequential decay.
 *  Generates total *excitation* energy, and relative n-n energy
 *  The generation happens in g4gen::RngTwoBodyGenerator, so we
 *  need to read the last generated values from the
 *  RNG used to pick the recoil excitation energy, when
 *  the reaction generator was run.
 */
 	assert (mInitial);
	assert (mRngEx);
	
	// Choose decay energies
	//
	// 2n separation energy
	const G4double s2n = 
		mInitial->M() - mFinalFragMass - GetNumberOfNeutrons()*kNeutronMass;
	const G4double edecay2n = s2n + mRngEx->GetRng2d()->GetLast().first;
	const G4double erel = mRngEx->GetRng2d()->GetLast().second;
	const G4double edecay1 = (edecay2n + erel)/2; // decay energy n1
	const G4double edecay2 = (edecay2n - erel)/2; // decay energy n2
	// Check both decay energies are valid
	if(edecay1 <= 0 || edecay2 <= 0) {
		if(GetVerboseLevel() > 1) {
			G4GWAR << "g4gen::TwoNeutronDecaySequential :: Not enough energy for decay "
						 << "(edecay1, edecay2 :: "  << edecay1 << ", " << edecay2 << G4endl;
		}
		return false;
	}
	
	////////////////////////////////////////
	// Evaportation of neutrons

	G4LorentzVector lvN1, lvF1, lvN2, lvF2;
	G4double intermediateMass = mInitial->MplusEx() - kNeutronMass - edecay1;

	// Do first neutron evaporation in COM frame
	{
		g4gen::NeutronEvaporation evap1(mInitial->MplusEx(), intermediateMass, kNeutronMass);
		evap1(&lvF1, &lvN1);
	}

	// Second neutron evaporation (COM frame)
	{
		g4gen::NeutronEvaporation evap2(intermediateMass, mFinalFragMass, kNeutronMass);
		evap2(&lvF2, &lvN2);
	}

	
  ////////////////////////////////////////
	// Boost into lab frame

	// First boost to frame where frag is at rest after first decay
	// n.b. neutron 1 (lvN1) is already in this frame
	lvF2.boost(lvF1.boostVector()); // FINAL fragment
	lvN2.boost(lvF1.boostVector()); // neutron from SECOND decay

	// Now boost all three into lab frame
	lvN1.boost(mInitial->Momentum().boostVector()); // neutron 1
	lvN2.boost(mInitial->Momentum().boostVector()); // neutron 2
	lvF2.boost(mInitial->Momentum().boostVector()); // final fragment

	
	SetFinal(0, mInitial->Momentum()); // initial nucleus, before decay
	SetFinal(1, lvF2);     // fragment
	SetFinal(2, lvN1);     // neutron 1
	SetFinal(3, lvN2);     // neutron 2

	return true;
}


/////////////////////////////////////////////////////////////////
// g4gen::TwoNeutronDecayDiNeutron
//

g4gen::TwoNeutronDecayDiNeutron::TwoNeutronDecayDiNeutron():
	g4gen::NeutronDecayIntermediate(2)
{ }

g4gen::TwoNeutronDecayDiNeutron::~TwoNeutronDecayDiNeutron()
{ }

std::unique_ptr<g4gen::Rng> g4gen::TwoNeutronDecayDiNeutron::CreateRngEx()
{
	const G4double ex = GetDecayParam("ex"); // excitation energy
	const G4double width = GetDecayParam("width"); // BW width
	const G4int A_i = GetDecayParam("A_i");
	g4gen::RngVolyaDiNeutronEx* rng = new g4gen::RngVolyaDiNeutronEx(ex, width, -18.7, A_i);
	mRngEx = rng;
	return std::unique_ptr<g4gen::Rng>(rng);
}

G4bool g4gen::TwoNeutronDecayDiNeutron::Generate()
{
	/** Code taken from st_reaction.cc out of st_mona simulation in
	 *  use by the MoNA collaboration. Adapted for this simulation.
	 */
	assert (mInitial);
	assert (mRngEx);

	////////////////////////////////////////////////
  // figure out decay energy, etc
	//
	/** \note Using 'Volya' dineutron decay.
	 *  Generates intrinsic and kinetic energies of dineutron.
	 *  The generation happens in g4gen::RngTwoBodyGenerator, so we
	 *  need to read the last generated values from the
	 *  RNG used to pick the recoil excitation energy, when
	 *  the reaction generator was run.
	 */
	// 2n separation energy
	const G4double s2n =
		mInitial->M() - mFinalFragMass - GetNumberOfNeutrons()*kNeutronMass;
	// total DECAY energy:  A -> [(A-2) + (2n)]
	const G4double ei2n = mRngEx->GetRng2d()->GetLast().first;
	// dineutron INTRINSIC energy
	const G4double ed2n = s2n + mRngEx->GetRng2d()->GetLast().second;
	// make sure edecay is valid ( > 0)
	if(ed2n <= 0) {
		if(GetVerboseLevel() > 1) {
			G4GERR << "g4gen::TwoNeutronDecayDiNeutron:: Not enough energy for decay!" << G4endl
						 << "MASSES (MBeam, ex, MF, 2*MN):: " << mInitial->M() << ", " 
						 << mInitial->Ex() << ", "
						 << mFinalFragMass << ", " << 2*kNeutronMass << G4endl
						 << "Edecay (2n), s2n, Eintrinsic (2n):" 
						 << ed2n << ", " << s2n << ", " << ei2n << G4endl;
		}
		return false;
	}

	////////////////////////////////////////////////
  // initial unbound state decay: A -> (A-2) + (2n)
	//
	G4LorentzVector lvFrag, lv2N; // fragment, dineutron
	{
		G4double eCM = 2*kNeutronMass + mFinalFragMass + ed2n; // total initial energy
		g4gen::NeutronEvaporation evap(eCM, mFinalFragMass, 2*kNeutronMass);
		evap(&lvFrag, &lv2N);
	}

	////////////////////////////////////////////////
  // di-neutron decay: 2n -> n+n
	//
	G4LorentzVector lvN1, lvN2; // neutron 1, neutron 2
	{
		G4double eCM = 2*kNeutronMass + ei2n; // total CM energy
		g4gen::NeutronEvaporation evap(eCM, kNeutronMass, kNeutronMass);
		evap(&lvN1, &lvN2);
	}
	
	////////////////////////////////////////////////
  // Boosting into the frame where the fragment is
	// at rest after the first decay
	//
  lvN1.boost(lv2N.boostVector());
  lvN2.boost(lv2N.boostVector());

	////////////////////////////////////////////////
  // final boost into lab frame + set return values
	//
  const G4ThreeVector& labBoost = mInitial->Momentum().boostVector();
  lvFrag.boost(labBoost);
  lvN1.boost(labBoost);
  lvN2.boost(labBoost);

	SetFinal(0, mInitial->Momentum()); // initial nucleus, before decay
	SetFinal(1, lvFrag);   // fragment
	SetFinal(2, lvN1);     // neutron 1
	SetFinal(3, lvN2);     // neutron 2

	return true;
}
#endif



///////////////////////////////////////////////////////////////
// g4gen::NeutronEvaporation
// 
void g4gen::NeutronEvaporation::operator() (G4LorentzVector* Frag, G4LorentzVector* Neut)
{
	// Method for simulating neutron evaporation in center of mass frame
	// Copied from st_mona code. Results still need to be boosted into the
	// lab frame.

	// total CM energy
	G4double eCM = mM0;

	// total neutron, frag energies
	G4double eN = (POW2(eCM) + POW2(mMn) - POW2(mMf)) / (2*eCM);
	G4double eF = (POW2(eCM) - POW2(mMn) + POW2(mMf)) / (2*eCM);

	// neutron, frag momenta
	G4double pN = POW2(eN) - POW2(mMn);
	pN = pN < 0 ? 0 : sqrt(pN);
	
	G4double pF = POW2(eF) - POW2(mMf);
	pF = pF < 0 ? 0 : -1*sqrt(pF); // n.b.: opposite direction

	// Momentum along z-direction
	Frag->set(0,0,pF,eF);
	Neut->set(0,0,pN,eN);

	// Generate & set angles
	G4double cosTheta = kRngUniform_neg1_pos1.Generate();
	G4double theta = acos(cosTheta);
	G4double phi = kRngUniform_0_2pi.Generate();

	Frag->setTheta(theta);
	Frag->setPhi(phi);
	Neut->setTheta(CLHEP::pi - theta);
	Neut->setPhi(phi + CLHEP::pi);
}


///////////////////////////////////////////////////////////////
// g4gen::NeutronDecayFactory
//

G4double g4gen::NeutronDecayFactory::GetDecayOption(G4String option) const
{
	std::map<G4String, G4double>::const_iterator it = mOptions.find(option);
	return it == mOptions.end() ? 0 : it->second;
}

namespace {
struct ToLower { void operator() (char& c) {
	const std::string up = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	const std::string lo = "abcdefghijklmnopqrstuvwxyz";
	size_t i = up.find(c);
	if(i  < up.size() ) c = lo[i];
} };
	
void to_lower(std::string& str) {
	std::for_each(str.begin(), str.end(), ::ToLower());
} }

g4gen::NeutronDecay* g4gen::NeutronDecayFactory::Create()
{
	g4gen::NeutronDecay* decay = 0;
	G4String decayType = this->GetDecayType();
	::to_lower(decayType);
	
	if(decayType == "1n") {
		decay = new g4gen::OneNeutronDecay();      
	}
	else if(decayType == "2nphasespace") {
		decay = new g4gen::TwoNeutronDecayPhaseSpace(FALSE); 
	}
	else if(decayType == "2nphasespacefsi") {
		decay = new g4gen::TwoNeutronDecayPhaseSpace(TRUE); 
	}
	// else if(decayType == "2ndineutron") {
	// 	decay = new g4gen::TwoNeutronDecayDiNeutron();
	// }
	// else if(decayType == "2nsequential") {
	// 	decay = new g4gen::TwoNeutronDecaySequential();
	// }
	else {
		G4GERR << "g4gen::NeutronDecayFactory::Create:: Invalid GetDecayType():: " << GetDecayType()
					 << G4endl;
		exit(1);
	}

	g4gen::NeutronDecayIntermediate *di = dynamic_cast<g4gen::NeutronDecayIntermediate*>(decay);
	if(di) {
		for(const auto& it : mOptions) {
			di->SetDecayParam(it.first, it.second);
		}
	}
	return decay;
}
