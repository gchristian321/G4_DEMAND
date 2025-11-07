//************************************************************************
//*                                                                      *
//*                          Define Geometry                             *
//*                                                                      *
//************************************************************************


#include "G4SystemOfUnits.hh"                //Geant4 
#include "G4RunManager.hh"

#include "DRAGONDetectorConstruction.hh"     //local


namespace DRAGON {

void DRAGONDetectorConstruction::ugeo_defin() 
     {
//C.
//C.
//C.
//C.--> Define the geometrical dimensions of the different parts
//C.                      of the apparatus.
//C.
       if (tubetype == 0 || (tubetype > 1 && tubetype < 7)) 
         {
          Rrms    =  6.0;         //! gas volume radius
          TLrms   = 88.0;         //! DETE Volume
//C.
//C.                          Apertures
//C.
          rent = 0.8;
          lent = 0.8;
          len1 = 7.6;
          riren1 = 0.4;
          rilen1 = 0.5;

          len2 = 8.7;
          riren2 = 0.495;
          rilen2 = 0.59;

          len3 = 10.;
          riren3 = 0.65;
          rilen3 = 0.75;
          
          lex1   =  0.34;
          rilex1 =  0.45;
          rirex1 =  0.59;

          lex2   =  1.65;
          rilex2 =  0.675;
          rirex2 =  0.675;

          lex3   =  3.9;
          rilex3 =  0.895;
          rirex3 =  1.03;

          lex4   = 15.25;
          rilex4 =  1.25;
          rirex4 =  1.80;
//C
//C.                              Z Positions
//C.
          zent[0] = -15.6;
          zent[1] =   0.0;
          zent[2] = -36.05;
          zent[3] =   0.0;
          zent[4] = -63.75;
          zent[5] =   0.0;
          
          zex[0] = 15.6;
          zex[1] =  0.0;
          zex[2] = 29.95;
          zex[3] =   0.0;
          zex[4] = 39.7;
          zex[5] =   0.0;
          zex[6] =   68.75;
          }
      else
         {
		  if (tubetype == 1) 
		     {
              Rrms    =  6.0;         //! gas volume radius
              TLrms   = 88.0;        //! DETE Volume

//C.
//C.                          Apertures
//C.
              rent = 0.8;
              lent = 0.8;
              len1 = 7.6;
              riren1 = 0.4;
              rilen1 = 0.5;

              len2 = 8.7;
              riren2 = 0.495;
              rilen2 = 0.59;

              len3 = 10.;
              riren3 = 0.65;
              rilen3 = 0.75;
              
              lex1 = 0.34;
              rilex1 = 0.64135;
              rirex1 = 0.71755;

              lex2 = 7.15;
              rilex2 = 1.16713;
              rirex2 = 1.51003;

              lex3 = 10.25;
              rilex3 = 1.83134;
              rirex3 = 2.36474;

              lex4 = 11.45;
              rilex4 = 2.72161;
              rirex4 = 3.26644;
//C
//C.                              Z Positions
//C.
              zent[0] =  -15.6;
              zent[1] =    0.0;
              zent[2] = -36.05;
              zent[3] =    0.0;
              zent[4] = -63.75;
              zent[5] =    0.0;
              
              zex[0] =  10.25;
              zex[1] =    0.0;
              zex[2] =  38.35;
              zex[3] =    0.0;
              zex[4] =  67.95;
              zex[5] =    0.0;
              zex[6] = 106.35;
              }
          else
		     {
			  G4cerr << "UNKNOWN GEOMETRY" << G4endl;
			  G4RunManager::GetRunManager()->AbortRun();
		      }
		  }
			 G4cout << "END ugeo_defin() [TLrms]: " << TLrms << G4endl;
		 }

}
