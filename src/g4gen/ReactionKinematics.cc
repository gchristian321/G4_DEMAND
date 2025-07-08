#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "g4gen/NuclearMasses.hh"
#include "g4gen/ReactionKinematics.hh"


g4gen::ReactionKinematics::ReactionKinematics():
	fHaveInputs(FALSE),
	fBeam(0,0,0,0), fTarget(0,0,0,0),
	fOutputMasses(0), fOutputs(0)
{ }

g4gen::ReactionKinematics::ReactionKinematics(const G4LorentzVector& beam, 
																							const G4LorentzVector& target,
																							const std::vector<G4double>& finalProductMasses)
{
	SetInputs(beam,target,finalProductMasses);
}

g4gen::ReactionKinematics::~ReactionKinematics()
{ }

void g4gen::ReactionKinematics::SetInputs(const G4LorentzVector& beam, 
																					const G4LorentzVector& target,
																					const std::vector<G4double>& finalProductMasses)
{
	fHaveInputs = TRUE;
	fBeam = beam;
	fTarget = target;
	fOutputMasses = finalProductMasses;
	fOutputs.resize(finalProductMasses.size());
	for(auto& x : fOutputs) { x.set(0,0,0,0); }
}

G4double g4gen::ReactionKinematics::GetProductMass(size_t i)
{
	try { 
		return fOutputMasses.at(i);
	} catch(std::range_error& e) {
		G4cerr << "WARNING:: g4gen::ReactionKinematics::GetProductMass:: Index " << i
					 << " is out of range (max: " << fOutputMasses.size() << ")" <<G4endl;
		return -1;
	}
}

const G4LorentzVector& g4gen::ReactionKinematics::GetProduct(size_t i) const
{
	try { 
		return fOutputs.at(i);
	} catch(std::range_error& e) {
		G4cerr << "WARNING:: g4gen::ReactionKinematics::GetProduct:: Index " << i
					 << " is out of range (max: " << fOutputs.size() << ")" <<G4endl;
		static G4LorentzVector dummy(0,0,0,0);
		return dummy;
	}
}

void g4gen::ReactionKinematics::BoostToCM(G4LorentzVector* v) const
{
	G4ThreeVector bv = (GetBeam() + GetTarget()).boostVector();
	v->boost(-bv);
}

G4LorentzVector g4gen::ReactionKinematics::GetProductCM(size_t i) const
{
	G4LorentzVector out = GetProduct(i);
	BoostToCM(&out);
	return out;
}

G4LorentzVector g4gen::ReactionKinematics::GetBeamCM() const
{
	G4LorentzVector out = GetBeam();
	BoostToCM(&out);	
	return out;
}

G4LorentzVector g4gen::ReactionKinematics::GetTargetCM() const
{
	G4LorentzVector out = GetTarget();
	BoostToCM(&out);
	return out;
}


G4bool g4gen::ReactionKinematics::SetProduct(size_t i, const G4LorentzVector& v)
{
	try { 
		fOutputs.at(i) = v;
		return true;
	} catch(std::range_error& e) {
		G4cerr << "WARNING:: g4gen::ReactionKinematics::SetProduct:: Index " << i
					 << " is out of range (max: " << fOutputs.size() << ")" <<G4endl;
		return false;
	}
}


G4bool g4gen::TwoBodyReactionKinematics::Calculate(G4double thetaCM, G4double phiCM)
{
	// Check inputs okay
	if(HaveInputs() == FALSE) {
		G4cerr << "ERROR:: g4gen::TwoBodyReactionKinematics::Calculate:: "
					 << "Inputs not set!" << G4endl;
		return FALSE;
	}
	if(GetNumProducts() != 2) {
		G4cerr << "ERROR:: g4gen::TwoBodyReactionKinematics::Calculate:: Wrong number of "
					 << "output product masses: got " << GetNumProducts()
					 << ", should have gotten 2" << G4endl;
		return FALSE;
	}
	G4double mSum = GetProductMass(0) + GetProductMass(1);
	if(GetBeam().e() + GetTarget().e() < mSum) {
		G4cerr << "ERROR:: g4gen::ReactionKinematics::Calculate:: "
					 << "Not enough energy for reaction!" << G4endl;
		G4cerr << "\tBeam, Target, Ejectile, Recoil:: " << GetBeam().e() << ", " << GetTarget().e()
					 << ", " << GetProductMass(0) << ", " << GetProductMass(1) << G4endl;
		return FALSE;
	}

	G4LorentzVector v1 = GetBeam();   // beam
	G4LorentzVector v2 = GetTarget(); // target
	G4LorentzVector vinv = v1+v2; // invariant
	G4ThreeVector bv = vinv.boostVector(); // boost to lab frame

	// Boost to center of mass
	G4LorentzVector v1cm = v1.boost(-bv);
	G4LorentzVector v2cm = v1.boost(-bv);

	// Reaction in CM frame
	G4double S = vinv.m2(); // invariant mass squared
	G4double M3 = GetProductMass(0);
	G4double M4 = GetProductMass(1);
	G4double pcm = sqrt((pow(S - M3*M3 - M4*M4, 2) - 4*M3*M3*M4*M4) / (4*S));

	G4LorentzVector v3(pcm*sin(CLHEP::pi - thetaCM)*cos(phiCM),
										 pcm*sin(CLHEP::pi - thetaCM)*sin(phiCM),
										 pcm*cos(CLHEP::pi - thetaCM),
										 sqrt(pcm*pcm + M3*M3));

	G4LorentzVector v4(pcm*sin(thetaCM)*cos(CLHEP::pi - phiCM),
										 pcm*sin(thetaCM)*sin(CLHEP::pi - phiCM),
										 pcm*cos(thetaCM),
										 sqrt(pcm*pcm + M4*M4));

	// Boost to LAB frame
	v3.boost(bv);
	v4.boost(bv);

	SetProduct(0, v3);
	SetProduct(1, v4);

	return TRUE;
}																	

