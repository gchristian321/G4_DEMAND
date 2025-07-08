#include <cassert>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "g4gen/PhaseSpace.hh"
#include "g4gen/Error.hh"
#include "g4gen/ReactionKinematics.hh"
#include "g4gen/ReactionGenerator.hh"


namespace {
template<typename T> T POW2(const T& t) { return t*t; }


// Helper function to generate beam + target distributions
// shared implementation common to different classes
void generate_beam_target(g4gen::Rng* fRngEbeam,
													g4gen::Particle& fP1,
													g4gen::Particle& fP2,
													g4gen::BeamEmittance *fEmX,
													g4gen::BeamEmittance *fEmY)
{
	// Generate incoming beam particle energy
	//
	G4double ebeam = fRngEbeam->GenerateAbove(0); // beam kinetic energy
	G4double pbeam = sqrt( POW2(fP1.M()+ebeam) - fP1.M2() ); // beam momentum;
	
	// Generate Beam Angle & Position //
	//
	G4ThreeVector pos(0,0,0);
	G4ThreeVector pbeam3(0,0,pbeam);
	G4int iloop = 0; 
	for (auto* em : { fEmX, fEmY }) {
		if(em) {
			g4gen::RngGaus2d rng(em->GetSigmaX(), em->GetSigmaTX(), em->GetRho());
			auto xtx = rng.Generate();
			G4double x  = xtx.first * mm  +  em->GetX0() * mm; // position
			pos[iloop] = x;

			G4double tx = xtx.second * mrad; // angle
			pbeam3[iloop] = pbeam*sin(tx/rad);
		}
		++iloop;
	}	
	pbeam3[2] = sqrt( POW2(pbeam) - POW2(pbeam3[0]) - POW2(pbeam3[1]) );
	fP1.SetPosition(pos);
	fP1.SetP3(pbeam3);
	assert( fabs(fP1.E() - fP1.M()) - ebeam < 1e-8 );
	assert( fabs(fP1.P() - pbeam) < 1e-8 );
	assert( fabs(fP1.Momentum().m() - fP1.M()) < 1e-8 );
	
	// Set target
	//
	fP2.SetP3XYZ(0,0,0);
	fP2.SetPosition(fP1.Position());
}

} // namespace


g4gen::TwoBodyReactionGenerator::TwoBodyReactionGenerator():
	fRngEbeam(0), fRngEx3(0), fRngEx4(0), fRngTheta(0), fRngPhi(0),
	fEmX(0), fEmY(0)
{ fTheta=fPhi=0; }

g4gen::TwoBodyReactionGenerator::~TwoBodyReactionGenerator()
{ }

void g4gen::TwoBodyReactionGenerator::SetBeamTargetEjectile(const G4String& beam, 
																														const G4String& target, 
																														const G4String& ejectile)
{
	fP1.SetNucleus(beam);
	fP2.SetNucleus(target);
	fP3.SetNucleus(ejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4);
}

void g4gen::TwoBodyReactionGenerator::SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																														G4int Ztrgt, G4int Atrgt,
																														G4int Zejectile, G4int Aejectile)
{
	fP1.SetNucleus(Zbeam, Abeam);
	fP2.SetNucleus(Ztrgt, Atrgt);
	fP3.SetNucleus(Zejectile, Aejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4);
}

void g4gen::TwoBodyReactionGenerator::SetEmittance(g4gen::BeamEmittance* emX, g4gen::BeamEmittance* emY)
{
	fEmX = emX;
	fEmY = emY;
}

void g4gen::TwoBodyReactionGenerator::SetRNGs ( std::initializer_list<g4gen::Rng*> rngList )
{
	/** List items are:	rngEbeam, rngEx3, rngEx4, rngTheta, rngPhi
	 */
	int i=0;
	g4gen::Rng** fRng[5] =
		{ &fRngEbeam, &fRngEx3, &fRngEx4, &fRngTheta, &fRngPhi };
	for(const auto& rng : rngList) { *(fRng[i++]) = rng; }
}

const g4gen::Particle& g4gen::TwoBodyReactionGenerator::GetReactant(G4int i) const
{
	switch(i) {
	case 1: return fP1;
	case 2: return fP2;
	case 3: return fP3;
	case 4: return fP4;
	default:
		{
			G4cerr << "ERROR:: g4gen::TwoBodyReactionGenerator::GetReactant:: Invalid indx: " << i
						 << ". Valid arguments are 1,2,3,4. Returning DUMMY vector: (0,0,0,0)" << G4endl;
			static g4gen::Particle dummy;
			return dummy;
		}
	}
}

G4bool g4gen::TwoBodyReactionGenerator::Generate()
{
	// Generate beam and target
	//
	generate_beam_target(fRngEbeam,fP1,fP2,fEmX,fEmY);

	// Generate Excitation Energies
	//
	G4double ex3 = fRngEx3 ? fRngEx3->GenerateAbove(0) : 0;
	fP3.SetEx(ex3);
	G4double ex4 = fRngEx4 ? fRngEx4->GenerateAbove(0) : 0;
	fP4.SetEx(ex4);
	
	// Generate CM angles
	//
	fTheta = fRngTheta ? fRngTheta->Generate() : 0;
	fPhi = fRngPhi ? fRngPhi->Generate() : 0;

	// Now calculate reaction
	G4double masses[2] = { fP3.MplusEx() , fP4.MplusEx() };
	g4gen::TwoBodyReactionKinematics 
		reacKin(fP1.Momentum(), fP2.Momentum(), std::vector<G4double>(masses, masses+2));
	G4bool isGood = reacKin.Calculate(fTheta, fPhi);
	if(isGood) {
		fP3.SetPosition(fP1.Position());
		fP4.SetPosition(fP1.Position());
		fP3.SetP3(reacKin.GetProduct(0).vect());
		fP4.SetP3(reacKin.GetProduct(1).vect());
	} else {
		fP3.SetPosition(0,0,0);
		fP4.SetPosition(0,0,0);
		fP3.SetP3XYZ(0,0,0);
		fP4.SetP3XYZ(0,0,0);
	}
	
	return isGood;
}



//===============================================
//== g4gen::NeutronPhaseSpaceReactionGenerator
//==

g4gen::NeutronPhaseSpaceReactionGenerator::NeutronPhaseSpaceReactionGenerator(G4int n_neut):
	fNeutrons(n_neut), 
	fRngEbeam(0),
	fEmX(0), 
	fEmY(0)
{ fTheta=fPhi=0; }
	
g4gen::NeutronPhaseSpaceReactionGenerator::~NeutronPhaseSpaceReactionGenerator()
{ }

void g4gen::NeutronPhaseSpaceReactionGenerator::SetBeamTargetEjectile(const G4String& beam, 
																																			const G4String& target, 
																																			const G4String& ejectile)
{
	fP1.SetNucleus(beam);
	fP2.SetNucleus(target);
	fP3.SetNucleus(ejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4 - fNeutrons.size());
	for(g4gen::Particle& P : fNeutrons) { P.SetNucleus(0, 1); }
}

void g4gen::NeutronPhaseSpaceReactionGenerator::SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																																			G4int Ztrgt, G4int Atrgt,
																																			G4int Zejectile, G4int Aejectile)
{
	fP1.SetNucleus(Zbeam, Abeam);
	fP2.SetNucleus(Ztrgt, Atrgt);
	fP3.SetNucleus(Zejectile, Aejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4 - fNeutrons.size());
	for(g4gen::Particle& P : fNeutrons) { P.SetNucleus(0, 1); }
}

void g4gen::NeutronPhaseSpaceReactionGenerator::SetEmittance(g4gen::BeamEmittance* emX, g4gen::BeamEmittance* emY)
{
	fEmX = emX;
	fEmY = emY;
}

void g4gen::NeutronPhaseSpaceReactionGenerator::SetRNGs( std::initializer_list<g4gen::Rng*> rngList )
{
	fRngEbeam = *rngList.begin();
}

const g4gen::Particle& g4gen::NeutronPhaseSpaceReactionGenerator::GetReactant(G4int i) const
{
	switch(i) {
	case 1: return fP1;
	case 2: return fP2;
	case 3: return fP3;
	case 4: return fP4;
	default:
		try { return fNeutrons.at(i-5); }
		catch (std::exception&)
		{
			G4cerr << "ERROR:: g4gen::TwoBodyReactionGenerator::GetReactant:: Invalid indx: " << i
						 << ". Valid arguments are 1,2,3,4. Returning DUMMY vector: (0,0,0,0)" << G4endl;
			static g4gen::Particle dummy;
			return dummy;
		}
	}
}


G4bool g4gen::NeutronPhaseSpaceReactionGenerator::Generate()
{
	generate_beam_target(fRngEbeam, fP1, fP2, fEmX, fEmY);

	// Generate phase space
	//
	g4gen::PhaseSpace gen;
	std::vector<G4double> finals = { fP3.M() , fP4.M() };
	for(const g4gen::Particle& p : fNeutrons) { finals.push_back(p.M()); }
	G4bool isGood = gen.SetDecay(fP1.Momentum() + fP2.Momentum(), finals.size(), &finals[0]);
	if(!isGood) 
	{ // not enough energy for decay!
		fP3.SetPosition(0,0,0);
		fP4.SetPosition(0,0,0);
		fP3.SetP3XYZ(0,0,0);
		fP4.SetP3XYZ(0,0,0);
		for(g4gen::Particle& P : fNeutrons) {
			P.SetPosition(0,0,0);
			P.SetP3XYZ(0,0,0);
		}
		G4GERR << "g4gen::NeutronPhaseSpaceReactionGenerator::Generate() :: "
					 << "Not enough energy for decay!\nEnergies are (in | {out} | out sum):\n\t"
			;
		G4double sum = 0; for(G4double mf : finals) { sum += mf; }
		G4cerr << fP1.Momentum().e() << " | { ";
		for(G4double mf : finals) { G4cerr << mf << " "; }
		G4cerr << "} | " << sum << G4endl;
	}
	else
	{	
		G4double relwt, ran;
		do {
			relwt = gen.Generate() / gen.GetWtMax();
			ran = g4gen::RngUniform().Generate();
		} while(relwt < ran);
	
		// Calculate CM angles
		//
		G4LorentzVector v3CM = *gen.GetDecay(0);
		G4ThreeVector boost3 = v3CM.boostVector();
		v3CM.boost(-boost3);
		fTheta = v3CM.theta();
		fPhi = v3CM.phi();

		// Set outputs
		//
		fP3.SetPosition(fP1.Position());
		fP4.SetPosition(fP1.Position());
		fP3.SetP3(gen.GetDecay(0)->vect());
		fP4.SetP3(gen.GetDecay(1)->vect());
		int iP = 2;
		for(g4gen::Particle& P : fNeutrons) {
			P.SetPosition(fP1.Position());
			P.SetP3(gen.GetDecay(iP++)->vect());
		}
	}
	
	return isGood;
}
