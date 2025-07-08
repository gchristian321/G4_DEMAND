#ifndef G4GenPhaseSpace_HEADER_INCLUDE_GUARD
#define G4GenPhaseSpace_HEADER_INCLUDE_GUARD

#include "G4LorentzVector.hh"


namespace g4gen {

/// Porting of ROOT's TGenPhaseSpace code to the GEANT4 architecture
/** Motivations are twofold:
 *   - Allow phase space calcs w/o ROOT
 *   - Use GEANT4 RNGs for phase space calcs, for consistency and for
 *     common seed w/ the rest of the simulation
 *  \attention The input units are all taken to be MeV, the GEANT4 default
 *   This is *DIFFERENT FROM THE ROOT VERSION* which uses GeV.
 */
class PhaseSpace {
private:
	G4int        fNt;             // number of decay particles
	G4double     fMass[18];       // masses of particles
	G4double     fBeta[3];        // betas of decaying particle
	G4double     fTeCmTm;         // total energy in the C.M. minus the total mass
	G4double     fWtMax;          // maximum weigth
	G4LorentzVector  fDecPro[18];  //kinematics of the generated particles

	G4double PDK(G4double a, G4double b, G4double c);

public:
	PhaseSpace(): fNt(0), fMass(), fBeta(), fTeCmTm(0.), fWtMax(0.) {}
	PhaseSpace(const PhaseSpace &gen);
	virtual ~PhaseSpace() {}
	PhaseSpace& operator=(const PhaseSpace &gen);

	bool          SetDecay(const G4LorentzVector &P, G4int nt, const G4double *mass, const char *opt="");
	G4double        Generate();
	G4LorentzVector *GetDecay(G4int n);

	G4int    GetNt()      const { return fNt;}
	G4double GetWtMax()   const { return fWtMax;}
};

}

#endif

