#ifndef DRAGONION_HH
#define DRAGONION_HH

#include "G4Ions.hh"

namespace DRAGON
{

class DRAGONIon : public G4Ions {
public:
    DRAGONIon(const G4String& name,
              G4double mass,
              G4double width,
              G4double charge,
              G4int iSpin,
              G4int iParity,
              G4int iConjugation,
              G4int iIsospin,
              G4int iIsospin3,
              G4int gParity,
              const G4String& type,
              G4int lepton,
              G4int baryon,
              G4int encoding,
              G4bool stable,
              G4double lifetime,
              G4DecayTable* decaytable,
              G4bool shortlived,
              G4int Z,
              G4int A);

    virtual ~DRAGONIon() = default;
};

}

#endif

