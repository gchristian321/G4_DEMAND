#include "DRAGONEventAction.hh"      //local

//C.
//C.    ******************************************************************
//C.    *                                                                *
//C.    *      GEANT3 user routine to control tracking of one event      *
//C.    *                                                                *
//C.    *   ==>Called by GRUN                                            *
//C.    *                                                                *
//C.    ******************************************************************
//C.

namespace DRAGON {


void DRAGONEventAction::gutrev()
     {
      n_flag = 0;
      e_detect = 0;
      for (G4int i = 0; i < 6; ++i) 
          {nloss[i] = 0;}
      }
}