#include "DRAGONEMField.hh"


namespace DRAGON {

void DRAGONEMField::mitray_zone(G4double ZB, G4double ZC, G4double Z11, G4double Z12, G4double Z21, G4double Z22, G4int& IZONE) const 
     {
//C
//C     This subroutine determines which zone of the magnet we are in
//C     (entrance fringe field, uniform region, exit fringe field, or 
//c     overlapping entrance/exit fringe fields).
//C
//C     INPUT (REAL*8)
//C     --------------
//C     ZB, ZC        z-coordinate in the B and C axis systems
//C     Z11, Z12      boundaries of entrance fringe field in B axis system
//C     Z21, Z22      boundaries of exit fringe field in C axis system
//c     
//C     OUTPUT (INTEGER*4)
//C     ------------------
//C     IZONE         = 0 for entrance/exit far field
//C                   = 1 for pure entrance fringe field (no overlap with 
//c                       exit fringe field)
//C                   = 2 for interior ("uniform field") region
//C                   = 3 for pure exit fringe field (no overlap with
//C                       entrance fringe field
//C                   = 4 for case where the entrance and exit fringe fields
//C                       overlap
//C                   = -1  error
//C
//C   
      bool LENTR,LEXIT;
//C
//C     DEFAULT VALUE OF -1 FOR IZONE
      IZONE=-1;
//C
//C----------------------------------------------------------------------------
//C
//C     Now figure out which zone we are in, from the B and C coordinates
//c
//C     Check for entrance/exit far zone
   
      if (ZB > Z11 || ZC > Z22) 
         {
          IZONE = 0;
          if (ldiag) 
             {std::cout << "FAR ENTRANCE/EXIT REGION" << std::endl;}
          return;
          }
//C
//C     Check for uniform field zone
      if (ZB <= Z12 && ZC <= Z21) 
         {
          IZONE = 2;
          if (ldiag) 
             {std::cout << "UNIFORM FIELD REGION" << std::endl;}
          return;
          }
//C
//C     Check for entrance or exit fringe field
      if (ZB <= Z11 && ZB > Z12) 
         {
          LENTR = true;
          IZONE = 1;
          } 
      else 
         {LENTR = false;}

      if (ZC <= Z22 && ZC > Z21) 
         {
          LEXIT = true;
          IZONE = 3;
          } 
      else 
         {LEXIT = false;}
//C      
//C     For the special case of a very short magnet, where the entrance
//c     and exit fringe fields overlap, and there is NO uniform field
//c     region, both LENTR and LEXIT will be true from the above tests.
//c     We signal this by setting IZONE=4
//c
      if (LENTR && LEXIT) 
         {IZONE = 4;}
//C
//C     LDIAG=.TRUE. means we want diagnostic printout at each step
//C
      if (IZONE == 1) 
         {
          if (ldiag) 
             {std::cout << "ENTRANCE FRINGE FIELD REGION" << std::endl;}
          }     
      else if (IZONE == 3) 
              {
               if (ldiag) 
                  {std::cout << "EXIT FRINGE FIELD REGION" << std::endl;}
               } else if (IZONE == 4) 
                         {
                          if (ldiag) 
                             {
                              std::cout << "OVERLAPPING ENTRANCE & EXIT FRINGE FIELDS (SHORT MAGNET)\n";
                              std::cout << "TOTAL FIELD = ENTRANCE + EXIT - UNIFORM FIELD" << std::endl;
                              }
                          } 
                       else if (IZONE == -1) 
                               {
                                std::cout << " **ERROR** IN GU_MITRAY_ZONE" << std::endl;
                                std::cout << "ERROR DETERMINING WHICH ZONE WE ARE IN" << std::endl;
                                std::cout << "ZB=" << std::setw(12) << std::fixed << std::setprecision(4) << ZB
                                          << " ZC=" << std::setw(12) << ZC
                                          << " Z11=" << std::setw(12) << Z11
                                          << " Z12=" << std::setw(12) << Z12
                                          << " Z21=" << std::setw(12) << Z21
                                          << " Z22=" << std::setw(12) << Z22 << std::endl;
                                //jstop  = 1;
                                //ieotri = 1;
                                exit(1);
                                }
     }
}
