#include "DRAGONTrackingAction.hh"    //Geant4
#include "DRAGONEventAction.hh"

namespace DRAGON {


void DRAGONTrackingAction::gutrak()
     {
//C.
//C.    ******************************************************************
//C.    *                                                                *
//C.    *       User routine to control tracking of one track            *
//C.    *                                                                *
//C.    *    ==>Called by : GTREVE                                       *
//C.    *                                                                *
//C.    ******************************************************************
//C.
      fEventAction->jslit = 0;
      fEventAction->jstop = 0;
	  }
	  
}