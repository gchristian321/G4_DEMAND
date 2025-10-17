#include "DemandHit.hh"
#include "DemandSD.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

G4ThreadLocal G4Allocator<DemandHit>* DemandHitAllocator;

DemandHit::DemandHit(
	G4int ID, G4ThreeVector pos, G4ThreeVector actualPos, G4double energy, G4double time,
	const G4ParticleDefinition* particle,
	G4VPhysicalVolume* physicalvolume) :
  G4VHit(), fID(ID), fPos(pos), fActualPos(actualPos), fEnergy(0), fTime(time),
	fParticleA(-1), fParticleZ(-1),
	fPhysicalVolume(physicalvolume) {

	AppendEnergy(energy, particle);
	
}

DemandHit::~DemandHit() {

}


void DemandHit::AppendEnergy(G4double edep, const G4ParticleDefinition* particle)
{
	auto particleAZ =
		std::make_pair(particle->GetAtomicMass(), particle->GetAtomicNumber());

	if(particle == G4Electron::Definition() ||
		 particle == G4Positron::Definition()) {

		particleAZ.first = 0;
		particleAZ.second = 1001;
	}

	auto it = fEnergyByParticle.find(particleAZ);
	if(it == fEnergyByParticle.end()) {
		fEnergyByParticle.emplace(particleAZ, edep);
	}
	else {
		it->second += edep;
	}

	fEnergy += edep;

	// re-calculate all quenched energies
	fEnergyQuenched = 0;
	G4double maxDeposit = 0;
	for(const auto& p : fEnergyByParticle) {
		G4double eQuench = DemandSD::CalculateQuenching(
			p.second, p.first.first, p.first.second);
		if(eQuench > 0) { fEnergyQuenched += eQuench; }
		if(eQuench > maxDeposit) {
			maxDeposit = eQuench;
			fParticleA = particleAZ.first;
			fParticleZ = particleAZ.second;
		}
	}
}


void DemandHit::AddAdditionalHitInVolume(
	G4ThreeVector pos, G4ThreeVector actualPos, G4double edep, G4double time,
	const G4ParticleDefinition* particle)
{
	// take time of earliest hit
	if(time < fTime) { 
		fTime = time;
		fPos = pos;
		fActualPos = actualPos;
	}

	// append deposited energy
	AppendEnergy(edep, particle);
}
