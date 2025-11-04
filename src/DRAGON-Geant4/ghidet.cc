#include "DRAGONSteppingAction.hh"

#include "geant3functions.hh"

namespace DRAGON {


void DRAGONSteppingAction::ghidet()
{
 std::cout << "CALLING ghidet" << std::endl;

//C.
//C *** Description: Give actions when the current particle is in detector
//C.

      G4int i, ihit;
      
      G4double hits[5], edep, tofin, tofout, x1[3], x2[3], xd[3];
      
      G4String chname;                                          //OJO Eliminar

      if(idtype <= 0) goto _999_;
      if (inwvol == 1 && istop != 0) goto _999_;
//C.
//C *** INWVOL = 1     store entrance quantities
//C.
      if(inwvol == 1)
        {
         edep  = 0.0;
         tofin = tofg;
         ucopy(vect,x1,3);  
         chname = names;		 
         gmtod(x1,xd,1,chname);  
//C this is to record depth instead of z        x1[2] = xd[2]
         goto _999_;
         }
//C.
//C *** INWVOL = 0    Accumulate EDEP
//C.
      if(inwvol == 0)
        {
         if(istop != 0)goto _10_;
         edep = edep + destep;
          goto _999_;
        }
//C.
//C *** INWVOL = 2   Close hit (also when ISTOP # 0)
//C.
   _10_: 

      edep = edep + destep;
      tofout = tofg;
      ucopy(vect,x2,3);       
      gmtod(x2,xd,1,chname);  
//C this is to record depth instead of z       x2[2] = xd[2]
//C.

      if(idtype == 1)
        {
         hits[0] = (x1[0]+x2[0])/2.;
         hits[1] = (x1[1]+x2[1])/2.;
         hits[2] = (x1[2]+x2[2])/2.;
         hits[3] = (tofin + tofout)/2. * 1.E9;
         hits[4] = 1000. * edep;
//C.
//C.--> Note: Use the factor used in file ugeom.f; 
//C.-->       subroutine UDET; variable fact_dedx
//C.
      }
  _999_:
  ;
}

}
