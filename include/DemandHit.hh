#ifndef DemandHit_h
#define DemandHit_h

#include <map>
#include <array>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4ParticleDefinition;


struct HitInfo_t {
	double edep;
	double edep_quenched;
	int num_hits;
	HitInfo_t(
		double edep_ = 0,
		double edep_quenched_ = 0,
		int num_hits_ = 0
		): edep(edep_),
			 edep_quenched(edep_quenched_),
			 num_hits(num_hits_)
		{}
};


class DemandHit : public G4VHit {
public:
  DemandHit(G4int,G4ThreeVector,G4ThreeVector,G4double,G4double,G4double,
					 const G4ParticleDefinition*,
						G4VPhysicalVolume*, G4int, G4double);
  virtual ~DemandHit();

  G4int GetID() const {return fID;};
  G4double GetTime() const {return fTime;};
  G4double GetEnergy() const {return fEnergy;};
  G4double GetEnergyQuenched() const {return fEnergyQuenched;};
  const G4ThreeVector& GetPosition() const {return fPos;};
	std::array<G4double,3> GetPositionAsArray() const
		{return {fPos.x(),fPos.y(),fPos.z()};}
  const G4ThreeVector& GetActualPosition() const {return fActualPos;};
	std::array<G4double,3> GetActualPositionAsArray() const {
		return {fActualPos.x(),fActualPos.y(),fActualPos.z()};}
	const G4VPhysicalVolume* GetPhysicalVolume() const {return fPhysicalVolume;};
	G4int GetParticleA() const { return fParticleA;}
	G4int GetParticleZ() const { return fParticleZ;}
	G4int GetParticleID() const { return fParticleID;}
	G4int GetParentID() const { return fParentID; }
	G4double GetParentEnergy() const { return fParentEnergy; }

	void AddAdditionalHitInVolume(
		G4ThreeVector pos, G4ThreeVector actualPos,
		G4double edep, G4double edep_quenched, G4double time,
		const G4ParticleDefinition*, G4int parentID, G4double parentEnergy);
	
  inline void* operator new(size_t);
  inline void operator delete(void*);

private:
	void AppendEnergy(G4double edep, G4double edep_quenched, const G4ParticleDefinition*);
	
private:
  G4int fID;
  G4ThreeVector fPos; // detector center
	G4ThreeVector fActualPos; // actual coordinates of the interaction
	//<pd->GetPDGEncoding(), <edep, edep_quenched> >
	std::map<G4int, HitInfo_t > fEnergyByParticle;
	G4double fEnergy;
	G4double fEnergyQuenched;
  G4double fTime;
	G4VPhysicalVolume* fPhysicalVolume;
	G4int fParticleA;
	G4int fParticleZ;
	G4int fParticleID;
	G4int fParentID;
	G4double fParentEnergy;
};

typedef G4THitsCollection<DemandHit> DemandHitsCollection;

extern G4ThreadLocal G4Allocator<DemandHit>* DemandHitAllocator;

inline void* DemandHit::operator new(size_t)
{
    if (!DemandHitAllocator) DemandHitAllocator =
																new G4Allocator<DemandHit>;
    return (void*)DemandHitAllocator->MallocSingle();
}

inline void DemandHit::operator delete(void* aHit)
{
    DemandHitAllocator->FreeSingle((DemandHit*) aHit);
}

#endif
