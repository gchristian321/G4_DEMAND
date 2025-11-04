#ifndef history_h
#define history_h 1

#include "G4Types.hh"

namespace DRAGON {
extern G4double E_int, E_rec, E_g[15], E_gp[15], cost_g[15], phi_g[15],
       cost_gp[15], cost_r, cosp_r, x_r, y_r, z_r, thet_r,
       xstop, ystop, zstop, xint, yint, zint, x, y, xp, yp, xtest[10], 
       ytest[10], 
       etest[10], beamtof;
extern G4int Nodec, react, recdet, dsssdpos;
}

#endif


