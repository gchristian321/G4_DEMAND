#ifndef beamcom_h
#define beamcom_h 1

#include "G4Types.hh"

namespace DRAGON
{
// Common for beam energy and emittance initialization
extern G4double e0beam, e0recoil, ex, ey, el;
extern G4double beamvel;
extern G4double recoilmom, beamo;
extern G4double beamenerg;

// Constant for amu in MeV
//extern G4double amumev = 0.93149432e3;  // Equivalent to 0.93149432 E+03
extern G4double amumev; 
}

#endif

