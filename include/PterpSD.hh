#ifndef PterpSD_h
#define PterpSD_h

#include "G4VSensitiveDetector.hh"
#include "PterpHit.hh"
#include <map>

class G4ParticleDefinition;

class PterpSD : public G4VSensitiveDetector {
public:
  PterpSD(G4String name);
  virtual ~PterpSD();
  
  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);

	static G4double CalculateQuenching(
		G4double edep, const G4ParticleDefinition* particle);
	
	static G4double CalculateQuenching(
		G4double edep, G4int A, G4int Z);

	void AddThreshold(G4int copyNo, G4double thresh)
		{
			auto emp = fThresholdMap.emplace(copyNo, thresh);
			if(emp.second == false) {
				throw std::logic_error(
					"PterpSD::AddThreshold -- Duplicate Copy Number!");
			}
		}

	G4double GetThreshold(G4int copyNo)
		{
			auto it = fThresholdMap.find(copyNo);
			if(it == fThresholdMap.end()){
				throw std::logic_error(
					"PterpSD::AddThreshold -- Looked for bad Copy Number!");
			}
			return it->second;
		}
			
	
private:
  PterpHitsCollection* fHitsCollection;
  G4int fHCID;
	std::map<G4int, G4double> fThresholdMap;
};

#endif
