#ifndef DemandSD_h
#define DemandSD_h

#include "G4VSensitiveDetector.hh"
#include "DemandHit.hh"
#include <map>

class G4ParticleDefinition;

class DemandSD : public G4VSensitiveDetector {
public:
  DemandSD(G4String name);
  virtual ~DemandSD();
  
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
					"DemandSD::AddThreshold -- Duplicate Copy Number!");
			}
		}

	G4double GetThreshold(G4int copyNo)
		{
			auto it = fThresholdMap.find(copyNo);
			if(it == fThresholdMap.end()){
				throw std::logic_error(
					"DemandSD::AddThreshold -- Looked for bad Copy Number!");
			}
			return it->second;
		}
			
	
private:
  DemandHitsCollection* fHitsCollection;
  G4int fHCID;
	std::map<G4int, G4double> fThresholdMap;
};

#endif
