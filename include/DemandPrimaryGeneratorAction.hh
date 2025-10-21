//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
// 
/// \file DemandPrimaryGeneratorAction.hh
/// \brief Definition of the DemandPrimaryGeneratorAction class

#ifndef DemandPrimaryGeneratorAction_h
#define DemandPrimaryGeneratorAction_h 1

#include <array>
#include <memory>
#include <utility>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

class DemandPrimaryGeneratorMessenger;
class G4ParticleGun;
class G4Event;
namespace g4gen {
class ReactionKinematics;
}

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class G4ParticleDefinition;

class DemandPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  DemandPrimaryGeneratorAction();    
  virtual ~DemandPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);

	void SetReactionPosition(const G4ThreeVector& pos);
  G4ThreeVector GetReactionPosition() const;
  
  
  // set methods
  void SetRandomFlag(G4bool value);
	void SetBeamA(int A) { fBeamA = A; }
	void SetBeamZ(int Z) { fBeamZ = Z; }
	void SetBeamEnergy(double E) { fBeamEnergy = E; }
	void SetReaction(const G4String& reaction) { fReaction = reaction; }
	void SetReactionAngdist(const G4String& angdist);
	void SetExRecoil(double ex) { fExRecoil = ex; }
	void SetBeamSigmaX(double sig) { fSigma[0] = sig; }
	void SetBeamSigmaY(double sig) { fSigma[1] = sig; }
	void SetBeamSigmaThetaX(double sig) { fSigma[2] = sig; }
	void SetBeamSigmaThetaY(double sig) { fSigma[3] = sig; }
	void SetSourceEnergyLimits(double low, double high)
		{ fSourceEnergies.reset(new std::pair<double, double>(low,high)); }
	
	G4double GetBeamEnergy() const { return fBeamEnergy; }
	G4double GetBeamSigmaX() const { return fSigma[0]; }
	G4double GetBeamSigmaY() const { return fSigma[1]; }
	G4double GetBeamSigmaThetaX() const { return fSigma[2]; }
	G4double GetBeamSigmaThetaY() const { return fSigma[3]; }
	
	G4ParticleDefinition* GetBeamDefinition()     const
		{return fBeamDefinition;		 }
	G4ParticleDefinition* GetTargetDefinition()   const
		{return fTargetDefinition;	 }
	G4ParticleDefinition* GetEjectileDefinition() const
		{return fEjectileDefinition;}
	G4ParticleDefinition* GetRecoilDefinition()   const
		{return fRecoilDefinition;  }
	g4gen::ReactionKinematics* GetReactionKinematics() const
		{ return fReactionGenerator; }

	const G4String& GetReaction() const { return fReaction; }
	
private:
	bool SetupBeam();
	bool SetupReaction(bool);
	void ShootBeam(G4Event*);
	void ShootReaction(G4Event*);
	void ShootGeant3(G4Event*);
	bool CheckThetaLimits(
		const G4LorentzVector& lv, const G4ThreeVector& pos);
	
private:
  G4ParticleGun*  fParticleGun; // G4 particle gun
	G4ThreeVector fReactionPosition;
	int fBeamA;
	int fBeamZ;
	double fBeamEnergy;
	G4String fReaction;
	G4String fReactionAngdist;
	double fExRecoil;
	std::array<G4double,4> fSigma;
	std::unique_ptr<std::pair<double,double> > fSourceEnergies;
	
	G4ParticleDefinition* fBeamDefinition;
	G4ParticleDefinition* fTargetDefinition;
	G4ParticleDefinition* fEjectileDefinition;
	G4ParticleDefinition* fRecoilDefinition;
	
	g4gen::ReactionKinematics* fReactionGenerator;

	DemandPrimaryGeneratorMessenger *fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
