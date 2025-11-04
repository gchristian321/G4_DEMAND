#ifndef gammahit_h
#define gammahit_h 1

#include "G4Types.hh"

namespace DRAGON {
// flag for each event, equals 1 if recoil reached ENDV
extern G4int recoil_hit_ENDV, 
      // num of BGOs triggered during event
      num_bgos_hit, 
      // ID# of leading BGO
      num_bgo_first, 
      // ID# of second leading BGO
      num_bgo_second, 
      // number of pairs produced in event
      pair_productions, 
      // number of BGOs hit by adding back nearest neighbors
      num_bgos_hit_ab, 
      // ID# of leading BGO after adding back nearest neighbors
      num_bgo_first_ab, 
      // ID# of second leading BGO after adding back nearest neighbors
      num_bgo_second_ab;

// summed energy deposited in all BGOs
extern G4double e_bgos_total, 
         // energy deposited in leading BGO
         e_bgo_first, 
         // energy deposited in second leading BGO
         e_bgo_second, 
         // energy deposited in leading BGO after adding back nearest neighbors
         e_bgo_first_ab, 
         // energy deposited in second leading BGO after adding back nearest neighbors
         e_bgo_second_ab, 
         // gamma time of flight
         gammatof, 
         // conversion energy
         e0_conv;

// adjacency_matrix, element (i,j) = 1 if BGOs i and j are neighbors
extern G4int adjacency_matrix[30][30];
}

#endif
