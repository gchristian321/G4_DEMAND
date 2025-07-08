#ifndef G4GEN_REACTION_KINEMATICS_HH
#define G4GEN_REACTION_KINEMATICS_HH
#include <memory>
#include <vector>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Rng.hh"
#include "BeamEmittance.hh"


namespace g4gen {

/// Class to calculate nuclear reaction kinematics
class ReactionKinematics {
public:
	ReactionKinematics();
	
	/// Set Initial Beam + Excited State Four-Vectors (LAB frame), and output masses [MeV]
	/// Outpus masses should include excitation energy
	ReactionKinematics(const G4LorentzVector& beam, 
												const G4LorentzVector& target,
												const std::vector<G4double>& finalProductMasses);
	/// Empty
	virtual ~ReactionKinematics();
	
	/// Set Initial Beam + Excited State Four-Vectors (LAB frame), and output masses [MeV]
	/// Outpus masses should include excitation energy
	void SetInputs(const G4LorentzVector& beam,
								 const G4LorentzVector& target, 
								 const std::vector<G4double>& finalProductMasses);

	/// Get Output Reaction product (LAB FRAME) With index i, corresponding to the index in the input mass vector
	const G4LorentzVector& GetProduct(size_t i) const;
	/// Get Beam (LAB FRAME)
	const G4LorentzVector& GetBeam() const { return fBeam; }
	/// Get Target (LAB FRAME)
	const G4LorentzVector& GetTarget() const { return fTarget; }
	/// Get Output Reaction product (CM FRAME) With index i, corresponding to the index in the input mass vector
	G4LorentzVector GetProductCM(size_t i) const;
	/// Get Beam (CM FRAME)
	G4LorentzVector GetBeamCM() const;
	/// Get Target (CM FRAME)
	G4LorentzVector GetTargetCM() const;
	// Have Inputs?
	G4bool HaveInputs() const { return fHaveInputs; }
	/// Get Output Product Masses
	G4double GetProductMass(size_t i);
	/// Get Number of output products
	G4int GetNumProducts() const { return fOutputs.size(); }

	/// Calculate Output Lorentz Vectors, for a given theta, phi (both in center of mass)
	/// Pure virtual, must be implemented in child class
	virtual G4bool Calculate(G4double thetaCM, G4double phiCM) = 0;

protected:
	G4bool SetProduct(size_t i, const G4LorentzVector& v);
	void BoostToCM(G4LorentzVector* v) const;
	
private:
	G4bool fHaveInputs;
	G4LorentzVector fBeam, fTarget;
	std::vector<G4double> fOutputMasses;
	std::vector<G4LorentzVector> fOutputs;
};


/// Class for Two Body Reaction Kinematics
class TwoBodyReactionKinematics : public ReactionKinematics {
public:
	TwoBodyReactionKinematics() :
		ReactionKinematics() { }
	TwoBodyReactionKinematics(const G4LorentzVector& beam,
															 const G4LorentzVector& target, 
															 const std::vector<G4double>& finalProductMasses) 
		:	ReactionKinematics(beam, target, finalProductMasses) {  }
	virtual ~TwoBodyReactionKinematics() { }
	virtual G4bool Calculate(G4double thetaCM, G4double phiCM);
};

}

#endif
