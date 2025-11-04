#include "DRAGONIon.hh"

namespace DRAGON
{

DRAGONIon::DRAGONIon(const G4String& name,
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
                     G4int A)
    : G4Ions(name, mass, width, charge, iSpin, iParity, iConjugation,
             iIsospin, iIsospin3, gParity, type, lepton, baryon,
             encoding, stable, lifetime, decaytable, shortlived)
{
    SetAtomicNumber(Z);
    SetAtomicMass(A);
}

}
