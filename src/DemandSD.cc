#include "DemandSD.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4RunManager.hh"
#include "DemandAnalysis.hh"
#include "DemandRunAction.hh"

#include "TGraph.h"

DemandSD::DemandSD(G4String name) :
  G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1) {
  G4String HCname = "demandCollection";
  collectionName.insert(HCname);

	// Read file of LaPlace Quenching data
	// T.A. Laplace et al 2020 JINST 15 P11020
	// http://doi.org/10.1088/1748-0221/15/11/P11020
	G4String quench_fname = "../Quenching_LaPlace.dat";
	fProtonQuenchingData = TGraph(quench_fname.c_str());
	if(fProtonQuenchingData.IsZombie() || fProtonQuenchingData.GetN() == 0){
		G4Exception("DemandSD", "QuenchFileMissing", FatalException,
                ("Could not read quenching file: " + quench_fname).c_str()
			);
	}	else {
		for(int i=0; i< fProtonQuenchingData.GetN(); ++i){
			fProtonQuenchingData.SetPointY(
				i, fProtonQuenchingData.GetPointY(i) * 0.477
				);
		}
	}
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
	// ---> GAC New 2026/01/22
	// ---> Need to account for "non ionizing" energy deposition
	//      where the particle is a neutron. This is highly quenched
	//      (effectively zero scintillatio light), so remove it from
	//      the total.
	//
	G4double edep_total = step->GetTotalEnergyDeposit();
	G4double edep = edep_total - step->GetNonIonizingEnergyDeposit();
	if (edep < 1 * CLHEP::eV) { return true; }
	
	G4double edep_quenched =
		this->CalculateQuenching(edep, step);
	
	G4StepPoint* preStepPoint = step->GetPreStepPoint();
	G4TouchableHistory* touchable
		= (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4double hitTime = preStepPoint->GetGlobalTime();
	G4double preEnergy = preStepPoint->GetKineticEnergy();
	G4ThreeVector position_actual = preStepPoint->GetPosition();
	G4int parentID = step->GetTrack()->GetParentID();

	G4int volumeDepth = 0;
	if(touchable->GetVolume()->GetName() == "DEMAND_scintPV") {
		volumeDepth = 2;
	}
	G4int copyNo = touchable->GetVolume(volumeDepth)->GetCopyNo();
	G4ThreeVector position = touchable->GetVolume(volumeDepth)->GetTranslation();
	
	bool alreadyHaveHitInVolume = false;
	for(DemandHit* existingHit : *(fHitsCollection->GetVector())) {
		if(existingHit->GetID() == copyNo) {
			if(touchable->GetVolume() != existingHit->GetPhysicalVolume()) {
				throw std::logic_error(
					"PhysicalVolumes are different but IDs are the same");
			}
			existingHit->AddAdditionalHitInVolume(
				position,position_actual,edep,edep_quenched,hitTime,
				step->GetTrack()->GetParticleDefinition(),
				parentID, preEnergy
				);
			alreadyHaveHitInVolume = true;
		}
	}

	if(!alreadyHaveHitInVolume) {
		DemandHit* hit = new DemandHit(
			copyNo,position,position_actual,edep,edep_quenched,hitTime,
			step->GetTrack()->GetParticleDefinition(),
			touchable->GetVolume(), parentID, preEnergy);
		fHitsCollection->insert(hit);
	}

	bool PrintAllHits = false;
	if(PrintAllHits){
		static G4int eventAboveZeroID = -1;
		static G4int eventIDLast = -1;
		const G4int eventID =
			G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
		if(eventID != eventIDLast){
			eventIDLast = eventID;
			++eventAboveZeroID;
		}
		int isProton = 0;
		if(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding() == 2212){
			isProton = 1;
		}
		G4int isRecdet = static_cast<const DemandRunAction*>(
			G4RunManager::GetRunManager()->GetUserRunAction())
			->GetRecdet();
		
		// event detector particle energy
		G4cout << "HIT: " << eventID << " " <<
			eventAboveZeroID << " " << copyNo << " " <<
			step->GetTrack()->GetParticleDefinition()->GetParticleName() << " " <<
			isProton << " " << edep << " " << edep_quenched << " " << isRecdet << G4endl;
	}

	auto analysisManager = DemandAnalysis::Instance();
	if(analysisManager->GetSaveStepTree()){
		auto particle = step->GetTrack()->GetParticleDefinition();
		analysisManager->AddStep(
			copyNo, edep, edep_quenched, step->GetStepLength(),
			particle->GetParticleName(), particle->GetPDGEncoding()
			);
	}



	
	
#if 0
	std::cout << "DEMAND HIT!!!\n";
	std::cout << "VOLUME:: " << touchable->GetVolume()->GetName() << std::endl;
	std::cout << "copyNo " << copyNo << "\n";
	std::cout << "position: (" << position.x() << ", " << position.y() << ", " << position.z() << ")\n";
	std::cout << "SD:: " << this->GetName() << std::endl;
#endif	

	return true;
}

G4double DemandSD::CalculateQuenching(
	G4double edep, const G4Step* step) const
{
	G4double light = 0;
	auto particle = step->GetTrack()->GetParticleDefinition();
	//
	// energy where protons switch from semi-emperical fit to
	// birks
	const G4double proton_crossover = 25.5084 * CLHEP::MeV;
	if (particle == G4Electron::Definition() ||
			particle == G4Positron::Definition())  {
		// --> no quenching for e-, e+
		light = edep;
	}
	else if (particle == G4Proton::Definition() && edep < proton_crossover) {
		const double RangeMin = fProtonQuenchingData.GetPointX(0);
		const double RangeMax = fProtonQuenchingData.GetPointX(
			fProtonQuenchingData.GetN() - 1 );

		bool UseQuenchingInterp = false; // --> true to use interpolated quenching
		
		if(UseQuenchingInterp && edep >= RangeMin && edep < RangeMax){
			// --> use interpolation of LaPlace data
			light = fProtonQuenchingData.Eval(edep);
		}
		else {
			// --> use semi-emperical fit
			//  a*E - b*(1.0 - np.exp(-c*E))
			const double a = 0.7787, b = 1.62564, c = 0.417876;
			light = a*edep - b*(1 - exp(-c*edep));
		}
	}
	else {
		// use birks eqn with birks constant set from fit
		// to LaPlace data [see GetScintillatorMaterial() in
		// DemandDetectorConstruction.cc ]
	 
		auto emSat = G4LossTableManager::Instance()->EmSaturation();
		if (!emSat) {
			throw std::runtime_error("DemandSD::CalculateQuenching --> no emSat...");
		}

    // Returns "visible" edep using Birks constant of the current material
    light = emSat->VisibleEnergyDepositionAtAStep(step);
	}
	
	return light > 0 ? light : 0;
}

G4double DemandSD::CalculateEnergyResolution(G4double e_MeVee)
{
	// GAC --> energy resolution from fits of DEMAND detectors
	// to 22Na, 7Be, 137Cs compton edges (see DEMAND NIM paper)
	//
	// FWHM/E(ee) v. E(MeVee) given by sqrt(A^2/x + B^2/x^2 + C^2)
	// A=   0.00(7), B=0.0618(23), C=0.1354(22)

	const G4double A = 0, B = 0.0618, C = 0.1354;
	const G4double FWHM = e_MeVee*sqrt(A*A/e_MeVee + pow(B/e_MeVee, 2) + C*C);
	return G4RandGauss::shoot(e_MeVee, FWHM/2.355);	
}


#if 0
// old quenching codes no longer used
//
G4double DemandSD::CalculateQuenching( // static //
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
	else if(Z == 14) {// Silicon
		light = 0;
	}
	else {
		light = 0;
		G4cerr << "WARNING: unrecognized particle (A, Z) = (" << A << ", " << Z << "), setting light = 0!" << G4endl;
	}
	return light > 0 ? light : 0;
}
#endif
