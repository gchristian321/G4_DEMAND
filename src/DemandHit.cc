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
	G4VPhysicalVolume* physicalvolume, G4int parentID, G4double parentKE) :
  G4VHit(), fID(ID), fPos(pos), fActualPos(actualPos),
	fEnergy(0), fEnergyQuenched(0), fTime(time),
	fPhysicalVolume(physicalvolume),
	fParticleA(-1), fParticleZ(-1), fParticleID(-1),
	fParentID(parentID), fParentEnergy(parentKE)
{

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
			particle->GetPDGEncoding(),
			HitInfo_t(edep, edep_quenched, 0)
			);
	}
	else {
		it->second.edep += edep;
		it->second.edep_quenched += edep_quenched;
		it->second.num_hits += 1;
	}

	fEnergy += edep;
	fEnergyQuenched += edep_quenched;
#if 0
	{
		// G4double light = 0;
		// if (particle->GetPDGEncoding() == 2212) {
		// 	const double a = 0.7787, b = 1.62564, c = 0.417876;
		// 	light = a*fEnergy - b*(1 - exp(-c*fEnergy));
		// 	fEnergyQuenched = light;
		// }
		
		G4cout << "Det " << fID << ", " <<particle->GetParticleName() << " "
					 << fEnergyByParticle[particle->GetPDGEncoding()].num_hits << ": "
					 << edep << ", " << edep_quenched << " --> "
					 << ", " << fEnergyByParticle[particle->GetPDGEncoding()].edep << ", "
					 << fEnergyByParticle[particle->GetPDGEncoding()].edep_quenched
					 << " --> " << fEnergy << ", " << fEnergyQuenched
//					 << " --> " << light
					 << G4endl;
	}
#endif
	// get particle A, Z of max energy deposition
	G4double maxDeposit = 0;
	for(const auto& p : fEnergyByParticle) {
		auto theParticle = G4ParticleTable::GetParticleTable()->
			FindParticle(p.first);
		if(p.second.edep_quenched > maxDeposit){
			maxDeposit = p.second.edep_quenched;
			fParticleA = theParticle->GetAtomicMass();
			fParticleZ = theParticle->GetAtomicNumber();
			fParticleID = theParticle->GetPDGEncoding();
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
	const G4ParticleDefinition* particle,
	G4int parentID, G4double parentKE
	)
{
	// take time of earliest hit
	if(time < fTime) { 
		fTime = time;
		fPos = pos;
		fActualPos = actualPos;
		fParentID = parentID;
		fParentEnergy = parentKE;
	}
	
	// append deposited energy
	AppendEnergy(edep, edep_quenched, particle);
}
