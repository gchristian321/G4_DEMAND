#ifndef PterpHit_h
#define PterpHit_h

#include <map>
#include <array>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4ParticleDefinition;

class PterpHit : public G4VHit {
public:
  PterpHit(G4int,G4ThreeVector,G4ThreeVector,G4double,G4double,
					 const G4ParticleDefinition*,
					 G4VPhysicalVolume*);
  virtual ~PterpHit();

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

	void AddAdditionalHitInVolume(
		G4ThreeVector pos, G4ThreeVector actualPos,
		G4double edep, G4double time,
		const G4ParticleDefinition*);
	
  inline void* operator new(size_t);
  inline void operator delete(void*);

private:
	void AppendEnergy(G4double edep, const G4ParticleDefinition*);
	
private:
  G4int fID;
  G4ThreeVector fPos; // detector center
	G4ThreeVector fActualPos; // actual coordinates of the interaction
	std::map<std::pair<G4int, G4int>, G4double> fEnergyByParticle; //<<A,Z>, edep>
	G4double fEnergy;
	G4double fEnergyQuenched;
  G4double fTime;
	G4VPhysicalVolume* fPhysicalVolume;
	G4int fParticleA;
	G4int fParticleZ;
};

typedef G4THitsCollection<PterpHit> PterpHitsCollection;

extern G4ThreadLocal G4Allocator<PterpHit>* PterpHitAllocator;

inline void* PterpHit::operator new(size_t)
{
    if (!PterpHitAllocator) PterpHitAllocator =
																new G4Allocator<PterpHit>;
    return (void*)PterpHitAllocator->MallocSingle();
}

inline void PterpHit::operator delete(void* aHit)
{
    PterpHitAllocator->FreeSingle((PterpHit*) aHit);
}

#endif
