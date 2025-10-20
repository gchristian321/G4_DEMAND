#include "DemandSD.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

DemandSD::DemandSD(G4String name) :
  G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1) {
  G4String HCname = "demandCollection";
  collectionName.insert(HCname);
}

DemandSD::~DemandSD() {

}

void DemandSD::Initialize(G4HCofThisEvent* hce) {
	fHitsCollection =
		new DemandHitsCollection(
			SensitiveDetectorName,collectionName[0]);
	
	if (fHCID<0) { 
		fHCID =
			G4SDManager::GetSDMpointer()->
			GetCollectionID(fHitsCollection); 
	}
	hce->AddHitsCollection(fHCID,fHitsCollection);
}

G4bool DemandSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
	G4double edep = step->GetTotalEnergyDeposit();
	if (fabs(edep) < 1e-6) { return true; }
	
	// edep = CalculateQuenching(
	// 	edep, step->GetTrack()->GetParticleDefinition());
	// if (fabs(edep) < 1e-6) { return true; }
	
	G4StepPoint* preStepPoint = step->GetPreStepPoint();
	G4TouchableHistory* touchable
		= (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4int copyNo = touchable->GetVolume()->GetCopyNo();
	G4double hitTime = preStepPoint->GetGlobalTime();
	G4ThreeVector position_actual = preStepPoint->GetPosition();
	G4ThreeVector position = touchable->GetVolume()->GetTranslation();
	
	bool alreadyHaveHitInVolume = false;
	for(DemandHit* existingHit : *(fHitsCollection->GetVector())) {
		if(existingHit->GetID() == copyNo) {
			if(touchable->GetVolume() != existingHit->GetPhysicalVolume()) {
				throw std::logic_error(
					"PhysicalVolumes are different but IDs are the same");
			}
			existingHit->AddAdditionalHitInVolume(
				position,position_actual,edep,hitTime,
				step->GetTrack()->GetParticleDefinition());
			alreadyHaveHitInVolume = true;
		}
	}

	if(!alreadyHaveHitInVolume) {
		DemandHit* hit = new DemandHit(
			copyNo,position,position_actual,edep,hitTime,
			step->GetTrack()->GetParticleDefinition(),
			touchable->GetVolume());
		fHitsCollection->insert(hit);
	}

#if 0
	std::cout << "TEXNEUT HIT!!!\n";
	std::cout << "VOLUME:: " << touchable->GetVolume()->GetName() << std::endl;
	std::cout << "copyNo " << copyNo << "\n";
	std::cout << "position: (" << position.x() << ", " << position.y() << ", " << position.z() << ")\n";
	std::cout << "SD:: " << this->GetName() << std::endl;
#endif	

	return true;
}

G4double DemandSD::CalculateQuenching(
	G4double edep, const G4ParticleDefinition* particle)
{
	int A = particle->GetAtomicMass();
	int Z = particle->GetAtomicNumber();
	if(particle == G4Electron::Definition() ||
		 particle == G4Positron::Definition()) {
		A = 0; Z = 1001;
	}

	return CalculateQuenching(edep,A,Z);
}

G4double DemandSD::CalculateQuenching(
	G4double edep, G4int A, G4int Z)
{
	G4double light = edep;
	if(A == 0 && abs(Z) == 1001) { // electron or positron
		light = edep;
	}
	else if(A >= 1 && Z == 1) { // p,d,t (NOTE - only correct for p!)
// Parameters from original menate_R
// these may have a problem at low energies (in fact, Ee/Ep becomes negative below 80 keV!)
//			light = 0.83*edep-2.82*(1-exp(-0.25*pow(edep,0.93)));

// Parameters from doi.org/10.1016/j.nima.2014.03.028 [also NE-213]
		// double edepKeV = edep*1e3;
		// double a0=0.80,a1=-2519,a2=3.68e-4,a3=0.96;
		// light = a0*edepKeV + a1*(1-exp(-a2*pow(edepKeV,a3))); // keVee
		// light*=1e-3; // MeVee

// Parameters from NIMA 792, p. 74 (2015)
// This one is specific to p-Terphenyl (valid ONLY up to ~6 MeV proton energy)
// Use original menate_R formula (for NE-213) above this			
		if(edep <= 6.) {
			light = 0.0122 + 0.0886*edep + 0.0772*pow(edep,2) - 8.27e-4*pow(edep,3) - 7.04e-4*pow(edep,4);
		} else {
			light = 0.83*edep-2.82*(1-exp(-0.25*pow(edep,0.93)));
		}
	}
	else if(A >= 3 && Z == 2) {  // 3He, alpha
		light = 0.41*edep-5.9*(1-exp(-0.065*pow(edep,1.01)));
	}
	else if(Z == 3)  { // Li
		light = 0.1795*( edep );   // Obtained from EXP fit of measured leading coeffs
	}
	else if(Z == 4) { // Be
		light = 0.0821*( edep );   // Obtained from EXP fit of measured leading coeffs
	}
	else if(Z == 5) { // B
		light = 0.0375*( edep );  // Obtained from EXP fit of measured leading coeffs
	}
	else if(Z == 6) { // Carbon
		light = 0.017*( edep );
	}
	return light > 0 ? light : 0;
}
