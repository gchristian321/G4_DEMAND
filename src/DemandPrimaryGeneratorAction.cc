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
/// \file DemandPrimaryGeneratorAction.cc
/// \brief Implementation of the DemandPrimaryGeneratorAction class

#include <fstream>

#include "DemandPrimaryGeneratorAction.hh"
#include "DemandPrimaryGeneratorMessenger.hh"
#include "DemandDetectorConstruction.hh"
#include "DemandAnalysis.hh"
#include "DemandRunAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "Randomize.hh"

#include "ReactionKinematics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandPrimaryGeneratorAction::DemandPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr),
	 fBeamA(1),
	 fBeamZ(0),
	 fBeamEnergy(1*MeV),
	 fReaction("none"),
	 fReactionAngdist("flat"),
	 fExRecoil(0),
	 fSourceEnergies(nullptr),
	 fBeamDefinition(0),
	 fTargetDefinition(0),
	 fEjectileDefinition(0),
	 fRecoilDefinition(0),
	 fReactionGenerator(0),
	 fMessenger(new DemandPrimaryGeneratorMessenger(this))
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
	SetupBeam();
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(1.*MeV);
  // Set gun position
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("DemandPrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  fParticleGun
    ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandPrimaryGeneratorAction::~DemandPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
	bool beamChanged = SetupBeam();
	fReaction.toLower();
	
	if(fReaction.empty() || fReaction == "none") {
		ShootBeam(anEvent);
	}	else if(fReaction == "geant3") {
		ShootGeant3(anEvent);
	} else {

		if(dynamic_cast<const DemandDetectorConstruction&>(
				 *(G4RunManager::GetRunManager()->GetUserDetectorConstruction())).
			 GetHaveTarget())
		{
			G4LogicalVolume* trgtLV =
				G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
			if ( trgtLV ) {			beamChanged = true;    }
		}

		if(fabs(GetBeamSigmaX()) >= 1e-5*mm &&
			 fabs(GetBeamSigmaY()) >= 1e-5*mm &&
			 fabs(GetBeamSigmaThetaX()) >= 1e-5*mrad &&
			 fabs(GetBeamSigmaThetaY()) >= 1e-5*mrad
			)
		{   beamChanged = true;   }
						
		
		SetupReaction(beamChanged);
		ShootReaction(anEvent);
	}
}

void DemandPrimaryGeneratorAction::ShootGeant3(G4Event* anEvent)
{
	G3Event g3evt;
	auto runAction = static_cast<const DemandRunAction*>(
		G4RunManager::GetRunManager()->GetUserRunAction() );
	runAction->GetGeant3Event( anEvent->GetEventID() , &g3evt );

	fParticleGun->SetParticlePosition( g3evt.fPosition );
	fParticleGun->SetParticleEnergy( g3evt.fKineticEnergy );
	fParticleGun->SetParticleMomentumDirection( g3evt.fMomentumDirection );

	fParticleGun->GeneratePrimaryVertex(anEvent);

	const G4double mass   = fParticleGun->GetParticleDefinition()->GetPDGMass();
	const G4double pmag   = sqrt(pow(g3evt.fKineticEnergy+mass,2) - mass*mass);
	DemandAnalysis::Instance()->SetGeneratedNeutron(
		G4LorentzVector(pmag*g3evt.fMomentumDirection, g3evt.fKineticEnergy+mass)
		);
	DemandAnalysis::Instance()->SetGeneratedRecoil(
		g3evt.fRecoilLorentzVector
		);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandPrimaryGeneratorAction::ShootBeam(G4Event* anEvent)
{
	G4LorentzVector momentum;

	if(fSourceEnergies.get()) {
		//
		// isotropic source w/ randomized energy
		fBeamEnergy = G4RandFlat::shoot(
			fSourceEnergies->first, fSourceEnergies->second);
		const G4ThreeVector pos = fParticleGun->GetParticlePosition();
		const G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
		const G4double cos_theta = 2*G4RandFlat::shoot() - 1;
		const G4double theta  = acos(cos_theta);
		const G4double phi = G4RandFlat::shoot()*2*CLHEP::pi;
		const G4double pmag = sqrt(pow(fBeamEnergy+mass,2) - mass*mass);
		momentum.set(pmag*sin(theta)*cos(phi),
								 pmag*sin(theta)*sin(phi),
								 pmag*cos_theta,
								 fBeamEnergy+mass);
	}
	else {
#if 1
		// diffuse source over detector surface
		//
		const G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
		const G4double pmag = sqrt(pow(fBeamEnergy+mass,2) - mass*mass);
		momentum.set(0,0,pmag,fBeamEnergy+mass);
		const DemandDetectorConstruction::Module_t *module = 0;
		if(!module) {
			module = dynamic_cast<const DemandDetectorConstruction&>(
				*(G4RunManager::GetRunManager()->GetUserDetectorConstruction())
				).GetModuleByOriginalSpecification(0);
			if(!module) {
				throw std::logic_error(
					"DemandPrimaryGeneratorAction::CheckThetaLimits -- "
					"Could not get module 0 from detector construction.");
			}
		}
		double xWidth =
			module->m_Voxnum[0] * module->m_Voxsize[0];
		double yWidth =
			module->m_Voxnum[1] * module->m_Voxsize[1];
		fParticleGun->SetParticlePosition(
			G4ThreeVector(
				G4RandFlat::shoot(-xWidth/2, +xWidth/2),
				G4RandFlat::shoot(-yWidth/2, +yWidth/2),
				0) );																			
	
#else
		// isotropic source over detector surface
		//
		const G4ThreeVector pos = fParticleGun->GetParticlePosition();
		do {
			const G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
			const G4double cos_theta = 2*G4RandFlat::shoot() - 1;
			const G4double theta  = acos(cos_theta);
			const G4double phi = G4RandFlat::shoot()*2*CLHEP::pi;
			const G4double pmag = sqrt(pow(fBeamEnergy+mass,2) - mass*mass);
			momentum.set(pmag*sin(theta)*cos(phi),
									 pmag*sin(theta)*sin(phi),
									 pmag*cos_theta,
									 fBeamEnergy+mass);
		} while(CheckThetaLimits(momentum, pos) == false);
#endif
	}
	
	fParticleGun->SetParticleEnergy(fBeamEnergy);
	fParticleGun->SetParticleMomentumDirection(
		momentum.vect().unit());

	fParticleGun->GeneratePrimaryVertex(anEvent);

	DemandAnalysis::Instance()->SetGeneratedNeutron(momentum);
	DemandAnalysis::Instance()->SetGeneratedRecoil(G4LorentzVector(0,0,0,0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandPrimaryGeneratorAction::ShootReaction(G4Event* anEvent)
{
	if(!fReactionGenerator) {
		throw std::logic_error(
			"fReactionGenerator not set in ShootBeamAndReaction");
	}

	G4LorentzVector momentum, recoil_momentum;
	auto gen_isotropic =
		[&]() {
			const G4double cos_theta = 2*G4RandFlat::shoot() - 1;
			const G4double theta  = acos(cos_theta);
			const G4double phi = G4RandFlat::shoot()*2*CLHEP::pi;

			bool success = fReactionGenerator->Calculate(theta, phi);
			if(!success){
				throw std::logic_error(
					"not enough energy for reaction!");
			}

			momentum = fReactionGenerator->GetProduct(0);
			recoil_momentum = fReactionGenerator->GetProduct(1);
		};									 
	if(fReactionAngdist == "isotropic_surface") {
		// isotropic reaction generation over detector surface
		//
		const G4ThreeVector pos = fParticleGun->GetParticlePosition();
		do {gen_isotropic();} while(CheckThetaLimits(momentum, pos) == false);
	}
	else if(fReactionAngdist == "flat" || fReactionAngdist == "isotropic") {
		gen_isotropic();
	}
	else {
		std::ifstream ifs(fReactionAngdist.c_str());
		if(!ifs.good()) {
			char buf[4096];
			sprintf(buf,"Angular distribution file \"%s\" does not exist.",
							fReactionAngdist.c_str());
			throw std::invalid_argument(buf);
		}
		ifs.close();

		throw std::runtime_error(
			"Non-isotropic angular distribution still needs implementation");
	}
	
	fParticleGun->SetParticleMomentumDirection(
		momentum.vect().unit());

	fParticleGun->SetParticleEnergy(
		momentum.e() - momentum.m());

	fParticleGun->GeneratePrimaryVertex(anEvent);

	DemandAnalysis::Instance()->SetGeneratedNeutron(momentum);
	DemandAnalysis::Instance()->SetGeneratedRecoil(recoil_momentum);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool DemandPrimaryGeneratorAction::CheckThetaLimits(
	const G4LorentzVector& lv, const G4ThreeVector& pos)
{
	// Check if the generated neutron momentum, at the generated
	// position will intersect the first layer of the detector
	// if so, return true, else return false.
	//
	// The intent is to avoid generating events where the neutron
	// misses the detector - to keep run time down while collecting
	// good statistics.
	//
	bool retval = false;

	const DemandDetectorConstruction::Module_t *module = 0;
	if(!module) {
		module = dynamic_cast<const DemandDetectorConstruction&>(
			*(G4RunManager::GetRunManager()->GetUserDetectorConstruction())
			).GetModuleByOriginalSpecification(0);
		if(!module) {
			throw std::logic_error(
				"DemandPrimaryGeneratorAction::CheckThetaLimits -- "
				"Could not get module 0 from detector construction.");
		}
	}
	
	G4ThreeVector vHit(0,0,0);

	double xWidth =
		module->m_Voxnum[0] * module->m_Voxsize[0];
	double yWidth =
		module->m_Voxnum[1] * module->m_Voxsize[1];
	
	// propagate lv to detector using method given here:
	// https://math.stackexchange.com/questions/100439/
	// determine-where-a-vector-will-intersect-a-plane
	// (first answer - David Mitra)
	//
	if(fabs(module->m_Phi) > 1e-6) {
		throw std::runtime_error(
			"Phi > 0 not yet supported for CheckThetaLimits()");
	}
	G4ThreeVector N ( // normal to detector plane
		sin(module->m_Theta),
		0,
		cos(module->m_Theta));
			
	G4ThreeVector a ( // position of detector centre
		N.x() * module->m_Dist,
		N.y() * module->m_Dist,
		N.z() * module->m_Dist);
	
	G4ThreeVector o = pos;       // position of neutron vector origin
	G4ThreeVector d = lv.vect(); // direction of neutron vector

	double tFactor = N[0]*d[0] + N[1]*d[1] + N[2]*d[2];
	double RHS = -(N[0]*(o[0]-a[0]) + N[1]*(o[1]-a[1]) + N[2]*(o[2]-a[2]));
	if(fabs(tFactor) < 1e-6)
	{
		// if 0 == 0 --> true (intersects at all points)
		// otherwise, 0 != 0 --> false (no intersection - parallel ray)
		retval = fabs(RHS) < 1e-6;
	}
	else
	{
		// we have an actual equation to solve for t
		double tVal = RHS/tFactor;
		if(tVal > 0) // ray traveling TOWARDS detector
		{
			vHit = G4ThreeVector (
				o[0] + d[0]*tVal, o[1] + d[1]*tVal, o[2] + d[2]*tVal);
			if(vHit.x() > -xWidth/2 && vHit.x() < +xWidth/2 &&
				 vHit.y() > -yWidth/2 && vHit.y() < +yWidth/2)
			{
				retval = true; // intersects in physical detecto volume
			}
		}
	}
	
	return retval;
}

bool DemandPrimaryGeneratorAction::SetupBeam()
{
	static int lastA = -1;
	static int lastZ = -1;
	if(fBeamA == lastA && fBeamZ == lastZ) {
		return false;
	}

	if(fBeamA == 1 && fBeamZ == 0) {
    fBeamDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
	}
	else if (fBeamA == 1 && fBeamZ == 1) {
    fBeamDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
	}
	else {
		fBeamDefinition = G4ParticleTable::GetParticleTable()->
			GetIonTable()->GetIon(fBeamZ,fBeamA,0.);
	}
	if(!fBeamDefinition) {
		throw std::invalid_argument("Bad beam specification");
	}
  fParticleGun->SetParticleDefinition(fBeamDefinition);
	lastA = fBeamA;
	lastZ = fBeamZ;
	return true;
}

bool DemandPrimaryGeneratorAction::SetupReaction(bool beamChanged)
{
	static G4String lastReaction = "dummy123";
	static double lastEx = -1;
	if(beamChanged == false && lastReaction == fReaction && lastEx == fExRecoil) {
		return false;
	}

	fReaction.toLower();
	if(fReaction.empty() || fReaction == "none") {
		lastReaction = fReaction;
		lastEx = fExRecoil;
		return true;
	}
	
	if(!fReactionGenerator) {
		fReactionGenerator = new g4gen::TwoBodyReactionKinematics();
	}

	if(false) {
	}

	else if(fReaction == "p,n" || fReaction == "(p,n)") {
		fTargetDefinition = G4Proton::Definition();
		fEjectileDefinition = G4Neutron::Definition();
	}
	else if(fReaction == "p,p" || fReaction == "(p,p)") {
		fTargetDefinition = G4Proton::Definition();
		fEjectileDefinition = G4Proton::Definition();
	}
	else if(fReaction == "p,d" || fReaction == "(p,d)") {
		fTargetDefinition = G4Proton::Definition();
		fEjectileDefinition = G4Deuteron::Definition();		
	}
	else if(fReaction == "p,t" || fReaction == "(p,t)") {
		fTargetDefinition = G4Proton::Definition();
		fEjectileDefinition = G4Triton::Definition();				
	}

	else if(fReaction == "d,n" || fReaction == "(d,n)") {
		fTargetDefinition = G4Deuteron::Definition();
		fEjectileDefinition = G4Neutron::Definition();
	}
	else if(fReaction == "d,p" || fReaction == "(d,p)") {
		fTargetDefinition = G4Deuteron::Definition();
		fEjectileDefinition = G4Proton::Definition();
	}
	else if(fReaction == "d,d" || fReaction == "(d,d)") {
		fTargetDefinition = G4Deuteron::Definition();
		fEjectileDefinition = G4Deuteron::Definition();		
	}
	else if(fReaction == "d,t" || fReaction == "(d,t)") {
		fTargetDefinition = G4Deuteron::Definition();
		fEjectileDefinition = G4Triton::Definition();				
	}
	
	else if(fReaction == "3He,n" || fReaction == "(3He,n)") {
		fTargetDefinition = G4He3::Definition();
		fEjectileDefinition = G4Neutron::Definition();
	}	
	else if(fReaction == "3He,p" || fReaction == "(3He,p)") {
		fTargetDefinition = G4He3::Definition();
		fEjectileDefinition = G4Proton::Definition();
	}	
	else if(fReaction == "3He,d" || fReaction == "(3He,d)") {
		fTargetDefinition = G4He3::Definition();
		fEjectileDefinition = G4Deuteron::Definition();		
	}	
	else if(fReaction == "3He,t" || fReaction == "(3He,t)") {
		fTargetDefinition = G4He3::Definition();
		fEjectileDefinition = G4Triton::Definition();				
	}
	
	else if(fReaction == "a,n" || fReaction == "(a,n)") {
		fTargetDefinition = G4Alpha::Definition();
		fEjectileDefinition = G4Neutron::Definition();				
	}
	else if(fReaction == "a,p" || fReaction == "(a,p)") {
		fTargetDefinition = G4Alpha::Definition();
		fEjectileDefinition = G4Proton::Definition();				
	}
	else if(fReaction == "a,d" || fReaction == "(a,d)") {
		fTargetDefinition = G4Alpha::Definition();
		fEjectileDefinition = G4Deuteron::Definition();				
	}
	else if(fReaction == "a,t" || fReaction == "(a,t)") {
		fTargetDefinition = G4Alpha::Definition();
		fEjectileDefinition = G4Triton::Definition();				
	}
	else {
		// attempt to figure out
		bool success = false;
		std::string target = "", ejectile = "";
		size_t icomma = fReaction.size();
		for(size_t ii=0; ii< fReaction.size(); ++ii) {
			if(fReaction[ii] == '(') ;
			else if(fReaction[ii] == ',') {icomma=ii; break;}
			else {target += fReaction[ii];}
		}
		if(icomma < fReaction.size()) {
			for(size_t ii=icomma+1; ii< fReaction.size(); ++ii) {
				if(fReaction[ii] != ')') ejectile += fReaction[ii];
			}
			
		}
		if(target.size() && ejectile.size()) {
			auto set_particle =
				[&](const std::string& name, G4ParticleDefinition*& def) {
					bool retval = true;
					if     (name == "n")  { def = G4Neutron::Definition();  }
					else if(name == "p")  { def = G4Proton::Definition();   }
					else if(name == "d")  { def = G4Deuteron::Definition(); }
					else if(name == "t")  { def = G4Triton::Definition();   }
					else if(name == "3He"){ def = G4He3::Definition();      }
					else if(name == "a")  { def = G4Alpha::Definition();    }
					else                  { retval = false;                 }
					return retval;
				};
			success = set_particle(target,fTargetDefinition) && 
				set_particle(ejectile,fEjectileDefinition);
			if(success) { G4cerr << target << " -->> " << ejectile << G4endl; }
			throw(1);
		}
		if(!success) {
			throw std::invalid_argument(
				std::string("unsupported reaction type: \"") +
				fReaction.c_str() + "\"");
		}
	}
	
	fParticleGun->SetParticleDefinition(fEjectileDefinition);

	// energy loss of beam
	G4double beam_energy_at_reaction_point = fBeamEnergy;

	if(dynamic_cast<const DemandDetectorConstruction&>(
			 *(G4RunManager::GetRunManager()->GetUserDetectorConstruction())).
		 GetHaveTarget())
	{
		G4LogicalVolume* trgtLV =
			G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
		if ( trgtLV ) {
			double thickness = dynamic_cast<G4Box&>(*(trgtLV->GetSolid())).
				GetZHalfLength()*2.;
		
			static G4EmCalculator* emCalc = 0;
			if(!emCalc) { emCalc = new G4EmCalculator(); }
			G4double dedx = emCalc->ComputeTotalDEDX(
				fBeamEnergy, fBeamDefinition, trgtLV->GetMaterial());
			G4double dx = G4RandFlat::shoot(0., thickness);
			beam_energy_at_reaction_point -= (dedx * dx);
		}
	}

	
	double mBeam = fBeamDefinition->GetPDGMass();	
	double pBeam = sqrt(pow(beam_energy_at_reaction_point + mBeam,2) - mBeam*mBeam);
	G4LorentzVector vBeam(0,0,pBeam,mBeam+beam_energy_at_reaction_point);
	// theta spread
	if(fabs(GetBeamSigmaThetaX()) >= 1e-5*mrad &&
		 fabs(GetBeamSigmaThetaY()) >= 1e-5*mrad)
	{
		const G4double thx = G4RandGauss::shoot(0, GetBeamSigmaThetaX());
		const G4double thy = G4RandGauss::shoot(0, GetBeamSigmaThetaY());
		const G4ThreeVector pdir(tan(thx), tan(thy), 1);
		vBeam = G4LorentzVector(pBeam*sin(pdir.theta())*cos(pdir.phi()),
														pBeam*sin(pdir.theta())*sin(pdir.phi()),
														pBeam*cos(pdir.theta()),
														mBeam+beam_energy_at_reaction_point);
	}

	// position spread
	if(fabs(GetBeamSigmaX()) >= 1e-5*mm &&
		 fabs(GetBeamSigmaY()) >= 1e-5*mm)
	{
		const G4double x = G4RandGauss::shoot(0, GetBeamSigmaX());
		const G4double y = G4RandGauss::shoot(0, GetBeamSigmaY());
		const G4double z = fParticleGun->GetParticlePosition().z();
		fParticleGun->SetParticlePosition(
			G4ThreeVector(x,y,z));
	}

	double mTarget = fTargetDefinition->GetPDGMass();
	G4LorentzVector vTarget(0,0,0,mTarget);

	double mEjectile = fEjectileDefinition->GetPDGMass();

	double recoilA =
		fBeamA + fTargetDefinition->GetAtomicMass() - fEjectileDefinition->GetAtomicMass();
	double recoilZ =
		fBeamZ + fTargetDefinition->GetAtomicNumber() - fEjectileDefinition->GetAtomicNumber();
	fRecoilDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(
		recoilZ,recoilA,0); // ground state of recoil
	double mRecoil = fRecoilDefinition->GetPDGMass() + fExRecoil;
	
	fReactionGenerator->SetInputs(
		vBeam, vTarget, {mEjectile, mRecoil});

	lastReaction = fReaction;
	lastEx = fExRecoil;
	return true;
}

void DemandPrimaryGeneratorAction::SetReactionAngdist(
	const G4String& angdist)
{
	fReactionAngdist = angdist;
}

void DemandPrimaryGeneratorAction::SetReactionPosition(const G4ThreeVector& pos) {
    fReactionPosition = pos;
}

G4ThreeVector DemandPrimaryGeneratorAction::GetReactionPosition() const {
    return fReactionPosition;
}
