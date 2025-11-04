#ifndef params_h
#define params_h 1

#include "G4Types.hh"

//----67----

namespace DRAGON {

G4double beam_mass_excess, recoil_mass_excess;
G4double part_width, gam_width, spin_stat_fac;
G4double level[16], life[16], ztarg;
G4double br[15][10], beamlifetime;

G4int rstate, md[15][10];

G4String beamtyp;
G4String rectyp;

}

#endif
