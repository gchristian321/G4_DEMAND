#include "DemandHit.hh"
#include "DemandSD.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleTable.hh"

G4ThreadLocal G4Allocator<DemandHit>* DemandHitAllocator;

DemandHit::DemandHit(
	G4int ID, G4ThreeVector pos, G4ThreeVector actualPos,
	G4double energy, G4double energy_quenched, G4double time,
	const G4ParticleDefinition* particle,
	G4VPhysicalVolume* physicalvolume) :
  G4VHit(), fID(ID), fPos(pos), fActualPos(actualPos),
	fEnergy(0), fEnergyQuenched(0), fTime(time),
	fParticleA(-1), fParticleZ(-1),
	fPhysicalVolume(physicalvolume) {

	AppendEnergy(energy, energy_quenched, particle);
	
}

DemandHit::~DemandHit() {

}


void DemandHit::AppendEnergy(
	G4double edep, G4double edep_quenched, const G4ParticleDefinition* particle)
{
	auto it = fEnergyByParticle.find(particle->GetPDGEncoding());
	if(it == fEnergyByParticle.end()) {
		fEnergyByParticle.emplace(
			particle->GetPDGEncoding(), std::make_pair(edep, edep_quenched));
	}
	else {
		it->second.first += edep;
		it->second.second += edep_quenched;
	}

	fEnergy += edep;
	fEnergyQuenched += edep_quenched;

	// get particle A, Z of max energy deposition
	G4double maxDeposit = 0;
	for(const auto& p : fEnergyByParticle) {
		auto theParticle = G4ParticleTable::GetParticleTable()->
			FindParticle(p.first);
		if(p.second.second > maxDeposit){
			maxDeposit = p.second.second;
			fParticleA = theParticle->GetAtomicMass();
			fParticleZ = theParticle->GetAtomicNumber();			
		}
	}
	
#if 0
	// re-calculate all quenched energies
	fEnergyQuenched = 0;
	G4double maxDeposit = 0;
	for(const auto& p : fEnergyByParticle) {
		auto theParticle = G4ParticleTable::GetParticleTable()->
			FindParticle(p.first);
		auto e_deposit = p.second;
		
		G4double eQuench = DemandSD::CalculateQuenching(
			e_deposit, theParticle);
		if(eQuench > 0) { fEnergyQuenched += eQuench; }
		if(eQuench > maxDeposit) {
			maxDeposit = eQuench;	
			fParticleA = theParticle->GetAtomicMass();
			fParticleZ = theParticle->GetAtomicNumber();
		}
#endif
#if 0
		G4cout << "Quenched Energies\n";
		G4cout << "num hits in volume: " << fEnergyByParticle.size() << "\n";
		G4cout << "Current particle: " << theParticle->GetParticleName() << "\n";
		G4cout << "MAX hit particle (A,Z): " << fParticleA << ", " << fParticleZ << "\n";
		G4cout << "eQuench, edep, fEnergyQuenched: " << eQuench << ", " << p.second << ", " << fEnergyQuenched << "\n";
		G4cout << "------------" << G4endl;
	}
#endif
}


void DemandHit::AddAdditionalHitInVolume(
	G4ThreeVector pos, G4ThreeVector actualPos,
	G4double edep, G4double edep_quenched, G4double time,
	const G4ParticleDefinition* particle)
{
	// take time of earliest hit
	if(time < fTime) { 
		fTime = time;
		fPos = pos;
		fActualPos = actualPos;
	}

	// append deposited energy
	AppendEnergy(edep, edep_quenched, particle);
}
