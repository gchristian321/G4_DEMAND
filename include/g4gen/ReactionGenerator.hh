#ifndef G4GEN_REACTION_HH
#define G4GEN_REACTION_HH
#include <initializer_list>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Rng.hh"
#include "Particle.hh"
#include "BeamEmittance.hh"

namespace g4gen {

class ReactionGenerator {
public:
	ReactionGenerator() { }
	virtual ~ReactionGenerator() { }
	virtual void SetBeamTargetEjectile(const G4String& beam, 
																		 const G4String& target, 
																		 const G4String& ejectile) = 0;
	virtual void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																		 G4int Ztrgt, G4int Atrgt,
																		 G4int Zejectile, G4int Aejectile) = 0;
	virtual void SetEmittance(BeamEmittance* emX, BeamEmittance* emY) = 0;
	virtual void SetRNGs( std::initializer_list<Rng*> rngList ) = 0;
	virtual G4double GetThetaCM() const = 0;
	virtual G4double GetPhiCM() const = 0;
	virtual const Particle& GetReactant(G4int i) const = 0;
	virtual const BeamEmittance* GetEmittanceX() const = 0;
	virtual const BeamEmittance* GetEmittanceY() const = 0;
	virtual const Rng* GetRngEbeam() const  = 0;
	virtual const Rng* GetRngEx3()   const  = 0;
	virtual const Rng* GetRngEx4()   const  = 0;
	virtual const Rng* GetRngTheta() const  = 0;
	virtual const Rng* GetRngPhi()   const  = 0;
	virtual G4bool Generate() = 0;
};


/// Class to GENERATE Two Body nuclear reactions
/// Includes beam emittance and excitation energy of recoil
class TwoBodyReactionGenerator : public ReactionGenerator {
public:
	TwoBodyReactionGenerator();
	virtual ~TwoBodyReactionGenerator();
	void SetBeamTargetEjectile(const G4String& beam, 
														 const G4String& target, 
														 const G4String& ejectile);
	void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
														 G4int Ztrgt, G4int Atrgt,
														 G4int Zejectile, G4int Aejectile);
	void SetEmittance(BeamEmittance* emX, BeamEmittance* emY);
	/// rngEbeam, rngEx3, rngEx4, rngTheta, rngPhi
	void SetRNGs( std::initializer_list<Rng*> rngList );
	G4double GetThetaCM() const { return fTheta; }
	G4double GetPhiCM() const { return fPhi; }
	const Particle& GetReactant(G4int i) const;
	const BeamEmittance* GetEmittanceX() const { return fEmX; }
	const BeamEmittance* GetEmittanceY() const { return fEmY; }
	const Rng* GetRngEbeam() const { return fRngEbeam; }
	const Rng* GetRngEx3()   const { return fRngEx3;   } /// EJECTILE
	const Rng* GetRngEx4()   const { return fRngEx4;   } /// RECOIL
	const Rng* GetRngTheta() const { return fRngTheta; }
	const Rng* GetRngPhi()   const { return fRngPhi;   }

	G4bool Generate();
private:
	Particle fP1, fP2, fP3, fP4; /// Beam, target, ejectile, recoil
	/// Beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	Rng *fRngEbeam, *fRngEx3, *fRngEx4, *fRngTheta, *fRngPhi;
	/// Beam emittance x, y
	BeamEmittance *fEmX, *fEmY;
	/// Theta, Phi from last generated Event
	G4double fTheta, fPhi;
};


/// N-body "neutron" phase space reaction generator
/** Generates final states using a n-body phase space
 *  calculation (e.g. TGenPhaseSpace from ROOT).
 *  The size of N is set based on the arguments to
 *  SetBeamTargetEjectile()
 */
class NeutronPhaseSpaceReactionGenerator : public ReactionGenerator {
public:
	NeutronPhaseSpaceReactionGenerator(G4int n_neut);
	virtual ~NeutronPhaseSpaceReactionGenerator();
	virtual void SetBeamTargetEjectile(const G4String& beam, 
																		 const G4String& target, 
																		 const G4String& ejectile);
	virtual void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																		 G4int Ztrgt, G4int Atrgt,
																		 G4int Zejectile, G4int Aejectile);
	virtual void SetEmittance(BeamEmittance* emX, BeamEmittance* emY);
	/// rngEbeam
	void SetRNGs( std::initializer_list<Rng*> rngList );
	virtual G4double GetThetaCM() const { return fTheta; }
	virtual G4double GetPhiCM() const { return fPhi;   }
	virtual const Particle& GetReactant(G4int i) const;
	const BeamEmittance* GetEmittanceX() const { return fEmX; }
	const BeamEmittance* GetEmittanceY() const { return fEmY; }
	virtual const Rng* GetRngEbeam() const  { return fRngEbeam; }
	virtual const Rng* GetRngEx3()   const  { return 0; }
	virtual const Rng* GetRngEx4()   const  { return 0; }
	virtual const Rng* GetRngTheta() const  { return 0; }
	virtual const Rng* GetRngPhi()   const  { return 0; }

	virtual G4bool Generate();
private:
  /// Beam, target, ejectile, recoil
	Particle fP1, fP2, fP3, fP4;
	/// Neutrons
	std::vector<Particle> fNeutrons;
	/// Beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	Rng *fRngEbeam;
	/// Beam emittance x, y
	BeamEmittance *fEmX, *fEmY;
	/// Theta, Phi from last generated Event
	G4double fTheta, fPhi;
};

}

#endif
