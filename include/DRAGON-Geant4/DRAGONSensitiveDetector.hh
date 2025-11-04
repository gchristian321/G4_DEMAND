#ifndef DRAGONSENSITIVEDETECTOR_HH
#define DRAGONSENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"

#include "DRAGONHit.hh"

namespace DRAGON
{ 

class DRAGONDetectorConstruction;

class DRAGONSensitiveDetector : public G4VSensitiveDetector {
public:
    DRAGONSensitiveDetector(const G4String& name,DRAGONDetectorConstruction*);
    virtual ~DRAGONSensitiveDetector();

    virtual void Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    virtual void EndOfEvent(G4HCofThisEvent* hce) override;

    const G4String& GetHitsCollectionName(size_t i = 0) const {return collectionName[i];}
	
	void ghidet(const G4Step*);
    void ghipmt(const G4Step*);

private:
    G4THitsCollection<DRAGONHit>* fHitsCollection;
    G4int fHitCollectionID;
	
	G4int in_new_vol, inwvol, itrtyp, n_detmate, numed, number, idtype, istop, nlevel, noCopy;
	G4String chname_nlevel;
	G4double hits[5],edep, tofg, tofin, tofout, x1[3], x2[3], xd[3], vect[3];
	bool savehit;
	
	DRAGONDetectorConstruction* DRAGON_det = nullptr;
	
	
};

}

#endif





