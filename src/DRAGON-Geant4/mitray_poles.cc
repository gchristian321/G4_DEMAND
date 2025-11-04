#include "DRAGONEMField.hh"
#include "G4RunManager.hh"


//C***********************************************************************
//C                          MULTI-POLES SUBROUTINE
//C***********************************************************************


namespace DRAGON {


void DRAGONEMField::mitray_poles(G4double* DATA, G4double* XPOS, G4double* BFLD) const 
     {
      //C
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C                                                                      C
      //C     Subroutine for multipoles, GEANT implementation of MIT-RAYTRACE  C
      //C     adapted from:                                                    C
      //C     Subroutine POLES (NO, NP, T, TP, NUM) by S. Kowalski             C
      //C                                                                      C
      //C      TC(1) to  TC(6) =  (  X,  Y,  Z, VX, VY, VZ )                   C
      //C     DTC(1) to DTC(6) =  ( VX, VY, VZ, VXDOT, VYDOT, VZDOT )          C
      //C                                                                      C
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C
      //C     Modification History
      //C     --------------------
      //C     Mar 26, 1998   S. Yen
      //C     Apr  5, 1999   S. Yen  Calculate overlapping entrance and exit fringe
      //C                            fields for very short multipole.
      //C 

      G4double LF1, LF2, LU1, L;
      G4double GRAD1C, GRAD2C, GRAD3C, GRAD4C, GRAD5C;
      G4double XA, YA, ZA;
      G4double RAD,BQD,BHX,BOC,BDC,BDD;
      G4double Z11,Z12,Z21,Z22;
      G4double FRH,FRO,FRD,FRDD;    
//C
//C     First zero the output B-field in case of abort
//C
      for(G4int i=0; i<3; i++) 
         {BFLD[i] = 0.;}
//C
//C     EXTRACT THE MULTIPOLE PARAMETERS FROM THE INPUT DATA ARRAY
//C
      LF1  = DATA[0];
      LU1  = DATA[1];
      LF2  = DATA[2];
      A    = DATA[9];
      B    = DATA[10];
      L    = DATA[11];
      RAD  = DATA[12];
      BQD  = DATA[13];
      BHX  = DATA[14];
      BOC  = DATA[15];
      BDC  = DATA[16];
      BDD  = DATA[17];
      Z11  = DATA[18];
      Z12  = DATA[19];
      Z21  = DATA[20];
      Z22  = DATA[21];
      FRH  = DATA[34];
      FRO  = DATA[35];
      FRD  = DATA[36];
      FRDD = DATA[37];
      DSH  = DATA[38];
      DSO  = DATA[39];
      DSD  = DATA[40];
      DSDD = DATA[41];
//C   DTF1= LF1/ VEL;
//C   DTF2= LF2/ VEL;
//C   DTU = LU1/ VEL;
      D = 2. * RAD;

      if (FRH == 0.0) FRH = 1.0;
      if (FRO == 0.0) FRO = 1.0;
      if (FRD == 0.0) FRD = 1.0;
      if (FRDD == 0.0) FRDD = 1.0;
      DH = FRH * D;
      DO = FRO * D; 
      DD = FRD * D;
      DDD = FRDD * D;
//C
//C     Extract the A-axis coordinates from the XPOS array
//C
      XA=XPOS[0];
      YA=XPOS[1];
      ZA=XPOS[2];
//C
//C     Calculate the B-axis coordinates (entrance VFB coordinates)
//C
      XB=-XA;
      YB=YA;
      ZB=A-ZA;
//C
//C     Calculate the C-axis coordinates (exit VFB coordinates)
//C
      XC=-XB;
      YC=YB;
      ZC=-ZB-L;
//C
//C
      BX = 0.;
      BY = 0.;
      BZ = 0.;
      BT = 0.;
      S = 0.;
//C
      if (ldiag) 
         {
          std::cout << std::string(50, '-') << std::endl;
          std::cout << " ENTER SUBROUTINE POLES" << std::endl;
          std::cout << " XA,YA,ZA= " << std::fixed << std::setprecision(3) << XA << " " << YA << " " << ZA << std::endl;
          }

      if (ldiag) 
         {
          std::cout << " XB,YB,ZB= " << std::fixed << std::setprecision(3) << XB << " " << YB << " " << ZB << std::endl;
          std::cout << " XC,YC,ZC= " << std::fixed << std::setprecision(3) << XC << " " << YC << " " << ZC << std::endl;
          std::cout << " Z11= " << std::fixed << std::setprecision(3) << Z11 << "  Z12= " << Z12 
                    << "  Z21= " << Z21 << "  Z22= " << Z22 << std::endl;
          std::cout << std::endl;
          }
//C
//C     Determine which zone of the multipole we are in by calling MITRAY_ZONE.
//c     IZONE=0  far entrance or exit (B=0)
//C     IZONE=1  entrance fringe field
//C     IZONE=2  uniform field region
//C     IZONE=3  exit fringe field
//C     IZONE=4  overlapping entrance and exit fringe fields for short magnet
//C              In this case, we calculate for the B-field at any point
//C              B=B(entrance) + B(exit) - B(uniform).  This simulates
//c              action of MIT-RAYTRACE in integrating *backwards* from end of
//c              the entrance fringe field, to the start of the exit fringe
//c              field, using the field of the uniform field region.
//C     IZONE=-1 error
//

  std::cout << "//////////////////////" << std::endl; 
  std::cout << "XPOS[0] " << XPOS[0] << std::endl; 
  std::cout << "XPOS[1] " << XPOS[1] << std::endl; 
  std::cout << "XPOS[2] " << XPOS[2] << std::endl; 
  std::cout << "ZB " << ZB << std::endl; 
  std::cout << "ZC " << ZC << std::endl; 
  std::cout << "Z11 " << Z11 << std::endl; 
  std::cout << "Z12 " << Z12 << std::endl; 
  std::cout << "Z22 " << Z22 << std::endl; 
  std::cout << "IZONE " << IZONE << std::endl; 
 
      mitray_zone(ZB,ZC,Z11,Z12,Z12,Z22,IZONE);
      
  std::cout << "XPOS[0] " << XPOS[0] << std::endl; 
  std::cout << "XPOS[1] " << XPOS[1] << std::endl; 
  std::cout << "XPOS[2] " << XPOS[2] << std::endl; 
  std::cout << "ZB " << ZB << std::endl; 
  std::cout << "ZC " << ZC << std::endl; 
  std::cout << "Z11 " << Z11 << std::endl; 
  std::cout << "Z12 " << Z12 << std::endl; 
  std::cout << "Z22 " << Z22 << std::endl; 
  std::cout << "IZONE " << IZONE << std::endl; 
//C
//C     For each zone, set variable IN accordingly, and call MITRAY_BPOLES
//C     to evaluate the magnetic field for that zone.
//C     Note that for IZONE=4 (overlapping entrance and exit fringe fields,
//C     for a very short multipole), the B-field that we want is
//C     B(total) = B(entrance fringe) + B(exit fringe) - B(uniform region)
//C

      if (IZONE == 0) 
         {
//C        
//C        ***********************************
//C        * FAR ENTRANCE OR EXIT ZONES, B=0 *
//C        ***********************************
//C
          if (ldiag) 
             {
              std::cout << std::string(50, '-') << std::endl;
              std::cout << " A SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
             }

          return;
          }
//C
//C
//C     Calculate the gradients for the various multipole components.
//C     The gradients evaluated below (with no minus signs) are correct when
//c     we are working in the C-axis system
//C     (uniform field and exit fringe field regions), but the signs of
//C     the quadrupole, octopole and dodecapole gradients must be changed
//C     when we work in the B-axis system (i.e. in the entrance fringe field
//C     region).
//C
      GRAD1C = BQD/RAD;
      GRAD2C = BHX/pow(RAD,2.0);
      GRAD3C = BOC/pow(RAD,3.0);
      GRAD4C = BDC/pow(RAD,4.0);
      GRAD5C = BDD/pow(RAD,5.0);
      
      if(IZONE==2 || IZONE==4)
        {
//C
//C        Come here for either uniform field region (IZONE=1) or
//C        overlapping entrance and exit fringe fields (IZONE=4).
//C
//C        *****************************************************
//C        * INTERIOR ("UNIFORM FIELD") ZONE FIELD CALCULATION *
//C        *****************************************************
//C
//C        Load the C-axis coordinates into array TC(i), which is what
//C        subroutine BPOLES is expecting for the uniform field region
//C        IN designates the region that we are in.
         IN=2;
         S=0.;
         TC[0]=XC;
         TC[1]=YC;
         TC[2]=ZC;
//C
//C        We are working in the C-axis system, so use C-axis gradients
         GRAD1=GRAD1C;
         GRAD2=GRAD2C;
         GRAD3=GRAD3C;
         GRAD4=GRAD4C;
         GRAD5=GRAD5C;

std::cout << "GRAD1 " << GRAD1 << std::endl;
std::cout << "GRAD2 " << GRAD2 << std::endl;
std::cout << "GRAD3 " << GRAD3 << std::endl;
std::cout << "GRAD4 " << GRAD4 << std::endl;
std::cout << "GRAD5 " << GRAD5 << std::endl;
         
//C
//C        CALL BPOLES TO CALCULATE B-FIELD IN C-AXIS SYSTEM
//C  
         mitray_bpoles();
         
         if (ldiag) 
            {
             std::cout << "UNIFORM FIELD REGION CALC." << std::endl;
             std::cout << " C SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << BX << " " << BY << " " << BZ << std::endl;
    	     std::cout << " A SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << BX << " " << BY << " " << BZ << std::endl;
             }

         if (IZONE == 2) 
            {
             BFLD[0] = BX; 
             BFLD[1] = BY;
             BFLD[2] = BZ;
             
             std::cout << "BFLD[0] " << BFLD[0] << std::endl;
             std::cout << "BFLD[1] " << BFLD[1] << std::endl;
             std::cout << "BFLD[2] " << BFLD[2] << std::endl;
             std::cout << "//////////////////////" << std::endl;
           
             return;
             }
         else 
             {
//c           IZONE=4  overlapping entrance & exit fringe fields
//c           We add the uniform field components BX,BY,BZ to the total
//c           field BFLD(i).
//c           Note minus sign, since B=B(entrance)+B(exit)-B(uniform)
  
              BFLD[0] = BFLD[0] - BX;
              BFLD[1] = BFLD[1] - BY;
              BFLD[2] = BFLD[2] - BZ;
              
              std::cout << "BFLD[0] " << BFLD[0] << std::endl;
              std::cout << "BFLD[1] " << BFLD[1] << std::endl;
              std::cout << "BFLD[2] " << BFLD[2] << std::endl;
              std::cout << "//////////////////////" << std::endl;
              }
         }

      if (IZONE == 1 || IZONE == 4) 
         {
//C
//C        Come here for either pure entrance fringe field (IZONE=1) or
//C        overlapping entrance and exit fringe fields (IZONE=4).
//C
//C        ******************************************
//C        * ENTRANCE FRINGE FIELD ZONE CALCULATION *
//C        ******************************************
//C
//C        Set IN=1 to designate entrance fringe field, and extract the
//C        fringing field coefficients
          IN = 1;
          C0 = DATA[22];
          C1 = DATA[23];
          C2 = DATA[24];
          C3 = DATA[25];
          C4 = DATA[26];
          C5 = DATA[27];

//C
//C       Load the B-axis coordinates into the array TC(i), which is
//C       what subroutine BPOLES is expecting
          TC[0] = XB;
          TC[1] = YB;
          TC[2] = ZB;

//C
//C       We are working in the B-axis system, so change the signs of
//C       the quadrupole, octopole and dodecapole gradients from those
//C       of the C-axis system.
          GRAD1 = -GRAD1C;
          GRAD2 = GRAD2C;
          GRAD3 = -GRAD3C;
          GRAD4 = GRAD4C;
          GRAD5 = -GRAD5C;

std::cout << "TTTTTTTTTTTTTTTT"<< std::endl;
std::cout << "C0 " << C0 << std::endl;
std::cout << "C1 " << C1 << std::endl;
std::cout << "C2 " << C2 << std::endl;
std::cout << "C3 " << C3 << std::endl;
std::cout << "C4 " << C4 << std::endl;
std::cout << "C5 " << C5 << std::endl;
std::cout << "TC[0] " << TC[0] << std::endl;
std::cout << "TC[1] " << TC[1] << std::endl;
std::cout << "TC[2] " << TC[2] << std::endl;
std::cout << "GRAD1 " << GRAD1 << std::endl;
std::cout << "GRAD2 " << GRAD2 << std::endl;
std::cout << "GRAD3 " << GRAD3 << std::endl;
std::cout << "GRAD4 " << GRAD4 << std::endl;
std::cout << "GRAD5 " << GRAD5 << std::endl;
std::cout << "TTTTTTTTTTTTTTTT"<< std::endl;

//C
//C
//C       Call BPOLES to calculate B-field in B-axis system
//C
          mitray_bpoles(); 

          if (ldiag) 
             {
              std::cout << "ENTRANCE FRINGE FIELD CALC." << std::endl;
              std::cout << " B SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << BX << " " << BY << " " << BZ << std::endl;
              std::cout << " A SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << -BX << " " << BY << " " << -BZ << std::endl;
              }
//C
//C       THE B-FIELD COMPONENTS BX,BY,BZ ARE IN THE B-AXIS SYSTEM; CONVERT TO
//C       THE A-AXIS SYSTEM.  STORE THEM IN OUTPUT ARRAY BFLD(i).
//C
          if (IZONE == 1) 
             {
              BFLD[0] = -BX;
              BFLD[1] = BY;
              BFLD[2] = -BZ;
              
              std::cout << "BFLD[0] " << BFLD[0] << std::endl;
              std::cout << "BFLD[1] " << BFLD[1] << std::endl;
              std::cout << "BFLD[2] " << BFLD[2] << std::endl;
              std::cout << "//////////////////////" << std::endl; 
        
              return;
              }
          else 
             {
//C          ! IZONE.EQ.4  overlapping entrance & exit fringe fields
//c          We add the entrance field components BX,BY,BZ to the total
//c          field BFLD(i), since B(total)=B(entrance)+B(exit)-B(uniform)
//c          Note the sign changes for BX and BZ since these are in the B-axis
//c          system and we need to change them to the A-axis system.
             BFLD[0] = BFLD[0] - BX;
             BFLD[1] = BFLD[1] + BY;
             BFLD[2] = BFLD[2] - BZ;
             
             std::cout << "BFLD[0] " << BFLD[0] << std::endl;
             std::cout << "BFLD[1] " << BFLD[1] << std::endl;
             std::cout << "BFLD[2] " << BFLD[2] << std::endl;
             std::cout << "//////////////////////" << std::endl; 
             }
         }

      if (IZONE == 3 || IZONE == 4) 
         {
//C
//C        Come here for either pure exit fringe field (IZONE=3) or
//C        overlapping entrance and exit fringe fields (IZONE=4).
//C
//C        **************************************
//C        * EXIT FRINGE FIELD ZONE CALCULATION *
//C        **************************************
//C
//C        Set IN=3 to designate exit fringe field, and extract fringing
//C        coefficients for the exit fringe field.
          IN = 3;
          C0 = DATA[28];
          C1 = DATA[29];
          C2 = DATA[30];
          C3 = DATA[31];
          C4 = DATA[32];
          C5 = DATA[33];
//C
//C        Load the C-axis coordinates into the array TC(i), which is
//C        what subroutine BPOLES is expecting
         TC[0] = XC;
         TC[1] = YC;
         TC[2] = ZC;
//C
//C        We are working in the C-axis system, so use C-axis gradients
         GRAD1 = GRAD1C;
         GRAD2 = GRAD2C;
         GRAD3 = GRAD3C;
         GRAD4 = GRAD4C;
         GRAD5 = GRAD5C;
//C
//C
//C        Call BPOLES to calculate B-field in C-axis system
//C
         mitray_bpoles();

         if (ldiag) 
            {
             std::cout << "EXIT FRINGE FIELD CALC." << std::endl;
             std::cout << " C SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << BX << " " << BY << " " << BZ << std::endl;
             std::cout << " A SYSTEM  BX,BY,BZ= " << std::fixed << std::setprecision(3) << BX << " " << BY << " " << BZ << std::endl;
             }
//C
//C        THE B-FIELD COMPONENTS BX,BY,BZ ARE IN THE C-AXIS SYSTEM; CONVERT TO
//C        THE A-AXIS SYSTEM (which is the same).  
//C        STORE THEM IN OUTPUT ARRAY BFLD(i).
//C
         if (IZONE == 3) 
            {
             BFLD[0] = BX;
             BFLD[1] = BY;
             BFLD[2] = BZ;
             
             std::cout << "BFLD[0] " << BFLD[0] << std::endl;
             std::cout << "BFLD[1] " << BFLD[1] << std::endl;
             std::cout << "BFLD[2] " << BFLD[2] << std::endl;
             std::cout << "//////////////////////" << std::endl;
        
             return;
             }
         else 
            {
//C          ! IZONE=4 overlapping entrance and exit fringe fields
//c          We add the exit field components BX,BY,BZ to the total
//c          field BFLD(i) since B(total)=B(entrance)+B(exit)-B(uniform)
             BFLD[0] = BFLD[0] + BX;
             BFLD[1] = BFLD[1] + BY;
             BFLD[2] = BFLD[2] + BZ;
             
             std::cout << "BFLD[0] " << BFLD[0] << std::endl;
             std::cout << "BFLD[1] " << BFLD[1] << std::endl;
             std::cout << "BFLD[2] " << BFLD[2] << std::endl;
             std::cout << "//////////////////////" << std::endl;
        
             return;
             }
         }

      if (IZONE == -1) 
         {
//C        UNKNOWN FIELD REGION, RETURN WITH B=0
          std::cout << "UNKNOWN MULTIPOLE FIELD REGION" << std::endl;
          std::cout << "!!! Abort current event !!!" << std::endl;
          G4RunManager::GetRunManager()->AbortEvent(); 
          
          //jstop = 1;
          //ieotri = 1;
      
          return;
          }
     }

void DRAGONEMField::mitray_bpoles() const
     {
      //C****
      //C**** CALCULATION OF MULTIPOLE(POLES) FIELD COMPONENTS
      //C****
      //C****
      //C****
      //C**** 2 - QUADRUPOLE  (GRAD1)
      //C**** 3 - HEXAPOLE    (GRAD2)
      //C**** 4 - OCTAPOLE    (GRAD3)
      //C**** 5 - DECAPOLE    (GRAD4)
      //C**** 6 - DODECAPOLE  (GRAD5)
      //C****
      //C****
      
      G4double X, Y, Z;
      G4double X2, X3, X4, X5, X6, X7;
      G4double Y2, Y3, Y4, Y5, Y6, Y7;
      G4double B2X, B2Y, B2Z;
      G4double B3X, B3Y, B3Z;
      G4double B4X, B4Y, B4Z;
      G4double B5X, B5Y, B5Z;
      G4double B6X, B6Y, B6Z;
      
      X = TC[0];
      Y = TC[1];
      Z = TC[2];
      X2 = X * X;
      X3 = X2 * X;
      X4 = X3 * X;
      X5 = X4 * X;
      X6 = X5 * X;
      X7 = X6 * X;
      Y2 = Y * Y;
      Y3 = Y2 * Y;
      Y4 = Y3 * Y;
      Y5 = Y4 * Y;
      Y6 = Y5 * Y;
      Y7 = Y6 * Y;
      
      std::cout << "VVVVVVVVVVVVVVVVVVV\n";
      std::cout << "IN " << IN << std::endl;
      std::cout << "X " << X << std::endl;
      std::cout << "Y " << Y << std::endl;
      std::cout << "Z " << Z << std::endl;
      std::cout << "X2 " << X2 << std::endl;
      std::cout << "X3 " << X3 << std::endl;
      std::cout << "X4 " << X4 << std::endl;
      std::cout << "X5 " << X5 << std::endl;
      std::cout << "X6 " << X6 << std::endl;
      std::cout << "X7 " << X7 << std::endl;
      std::cout << "Y2 " << Y2 << std::endl;
      std::cout << "Y3 " << Y3 << std::endl;
      std::cout << "Y4 " << Y4 << std::endl;
      std::cout << "Y5 " << Y5 << std::endl;
      std::cout << "Y6 " << Y6 << std::endl;
      std::cout << "Y7 " << Y7 << std::endl;
      std::cout << "VVVVVVVVVVVVVVVVVVV\n";
      
      if (IN == 2)
         {
          B2X = GRAD1 * Y;
          B2Y = GRAD1 * X;
          B3X = GRAD2 * 2.0 * X * Y;
          B3Y = GRAD2 * (X2 - Y2);
          B4X = GRAD3 * (3.0 * X2 * Y - Y3);
          B4Y = GRAD3 * (X3 - 3.0 * X * Y2);
          B5X = GRAD4 * 4.0 * (X3 * Y - X * Y3);
          B5Y = GRAD4 * (X4 - 6.0 * X2 * Y2 + Y4);
          B6X = GRAD5 * (5.0 * X4 * Y - 10.0 * X2 * Y3 + Y5);
          B6Y = GRAD5 * (X5 - 10.0 * X3 * Y2 + 5.0 * X * Y4);
          BX = B2X + B3X + B4X + B5X + B6X;
          BY = B2Y + B3Y + B4Y + B5Y + B6Y;
          BZ = 0.;
          BT = std::sqrt(BX * BX + BY * BY);
          
          std::cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n";
          std::cout << "IN " << IN << std::endl;
          std::cout << "BX " << BX << std::endl;
          std::cout << "BY " << BY << std::endl;
          std::cout << "BT " << BT << std::endl;
          std::cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n";
               
          return;
          }  
      else 
         {         
          S = Z/D; 
          std::cout << "GGGGGGGGGGGGANTESGGGGGGGGGGGGGGGGGGGGG\n";
          std::cout << "IN " << IN << std::endl;
          std::cout << "Z " << Z << std::endl;
          std::cout << "D " << D << std::endl;          
          std::cout << "G1 " << G1 << std::endl;
          std::cout << "G2 " << G2 << std::endl;
          std::cout << "G3 " << G3 << std::endl;
          std::cout << "G4 " << G4 << std::endl;
          std::cout << "G5 " << G5 << std::endl;
          std::cout << "G6 " << G6 << std::endl;
          std::cout << "DSH " << DSH << std::endl; 
          std::cout << "S " << S << std::endl; 
          std::cout << "RE " << RE << std::endl;
          std::cout << "GGGGGGGGGGANTESGGGGGGGGGGGGGGGGGGGGG\n";
        
          mitray_bpls(2, D, S, RE, G1, G2, G3, G4, G5, G6);  //AQUI ESTA EL PROBLEMA

          B2X = GRAD1 * (RE * Y - (G2 / 12.) * (3. * X2 * Y + Y3) +
                       (G4 / 384.) * (5. * X4 * Y + 6. * X2 * Y3 + Y5) -
                       (G6 / 23040.) * (7. * X6 * Y + 15. * X4 * Y3 + 9. * X2 * Y5 + Y7));
          B2Y = GRAD1 * (RE * X - (G2 / 12.) * (X3 + 3. * X * Y2) +
                       (G4 / 384.) * (X5 + 6. * X3 * Y2 + 5. * X * Y4) -
                       (G6 / 23040.) * (X7 + 9. * X5 * Y2 + 15. * X3 * Y4 + 7. * X * Y6));
          B2Z = GRAD1 * (G1 * X * Y - (G3 / 12.) * (X3 * Y + X * Y3) +
                       (G5 / 384.) * (X5 * Y + 2. * X3 * Y3 + X * Y5));
          
          std::cout << "GGGGGGGGGGDESPUESGGGGGGGGGGGGGGGGGGGGG\n"; 
          std::cout << "RE " << RE << std::endl;
          std::cout << "G1 " << G1 << std::endl;
          std::cout << "G2 " << G2 << std::endl;
          std::cout << "G3 " << G3 << std::endl;
          std::cout << "G4 " << G4 << std::endl;
          std::cout << "G5 " << G5 << std::endl;
          std::cout << "G6 " << G6 << std::endl;
          std::cout << "B2X " << B2X << std::endl;
          std::cout << "B2Y " << B2Y << std::endl;
          std::cout << "B2Z " << B2Z << std::endl;
          std::cout << "GGGGGGGGGGDESPUESGGGGGGGGGGGGGGGGGGGGG\n";

//C****
//C**** HEXAPOLE
//C****
          SS = Z/DH  + DSH;
          mitray_bpls( 3, DH, SS, RE, G1, G2, G3, G4, G5, G6 );
          B3X = GRAD2*( RE*2.*X*Y - (G2/48.)*(12.*X3*Y + 4.*X*Y3 ) );
          B3Y = GRAD2*( RE*(X2-Y2) - (G2/48.)*(3.*X4 + 6.*X2*Y2 - 5.*Y4 ) );
          B3Z = GRAD2*( G1*(X2*Y - Y3/3.) - (G3/48.)*(3.*X4*Y+2.*X2*Y3-Y5));
//C****
//C**** OCTAPOLE
//C****
          SS = Z/DO  + DSO;
          mitray_bpls( 4, DO, SS, RE, G1, G2, G3, G4, G5, G6 );
          B4X = GRAD3*( RE*(3.*X2*Y - Y3) - (G2/80.)*(20.*X4*Y - 4.*Y5 ) );
          B4Y = GRAD3*( RE*(X3 - 3.*X*Y2) - (G2/80.)*(4.*X5-20.*X*Y4 ) );
          B4Z = GRAD3*G1*(X3*Y - X*Y3 );
//C****
//C**** DECAPOLE
//C****
          SS = Z/DD  + DSD;
          mitray_bpls( 5, DD, SS, RE, G1, G2, G3, G4, G5, G6 );
          B5X = GRAD4*RE*(4.*X3*Y - 4.*X*Y3);
          B5Y = GRAD4*RE*(X4 - 6.*X2*Y2 + Y4 );
          B5Z = GRAD4*G1*(X4*Y - 2.*X2*Y3 + Y5/5. );
//C****
//C**** DODECAPOLE
//C****
          SS = Z/DDD + DSDD;
          mitray_bpls( 6, DDD,SS, RE, G1, G2, G3, G4, G5, G6 );
          B6X = GRAD5*RE*(5.*X4*Y - 10.*X2*Y3 + Y5 );
          B6Y = GRAD5*RE*(X5 - 10.*X3*Y2 + 5.*X*Y4 );
          B6Z = 0.;
//C****
//C**** TOTAL FIELD
//C****
          BX = B2X + B3X + B4X + B5X + B6X;
          BY = B2Y + B3Y + B4Y + B5Y + B6Y;
          BZ = B2Z + B3Z + B4Z + B5Z + B6Z;
          BT = std::sqrt( BX*BX + BY*BY + BZ*BZ );
           
          std::cout << "VVVVVVVVVVVVVVVVVVVVVVVV\n";
          std::cout << "IN " << IN << std::endl;
          std::cout << "BX " << BX << std::endl;
          std::cout << "BY " << BY << std::endl;
          std::cout << "BT " << BT << std::endl;
          std::cout << "VVVVVVVVVVVVVVVVVVVVVVVV\n";
          
          
          return;
          }
     }
//C
//C=======================================================================
//C


void DRAGONEMField::mitray_bpls(G4int IGP,G4double D,G4double S,G4double& RE, G4double& G1,G4double& G2,G4double& G3,G4double& G4,G4double& G5,G4double& G6) const
     {
      G4double S2,S3,S4,S5,CS;
      G4double CP1,CP2,CP3,CP4,CP5;
      G4double E,ERE,ERE1,ERE2,ERE3,ERE4,ERE5,ERE6;
      G4double CP12,CP13,CP14,CP22,CP15,CP16,CP23,CP32;
      
      std::cout << "--------------------\n";
      std::cout << "IGP " << IGP << std::endl;
      std::cout << "D " << D << std::endl;
      std::cout << "S " << S << std::endl;
      std::cout << "RE " << RE << std::endl;
      std::cout << "G1 " << G1 << std::endl;
      std::cout << "G2 " << G2 << std::endl;
      std::cout << "G3 " << G3 << std::endl;
      std::cout << "G4 " << G4 << std::endl;
      std::cout << "G5 " << G5 << std::endl;  
      std::cout << "G6 " << G6 << std::endl;      
      std::cout << "--------------------\n";

      S2 = S * S;
      S3 = S2 * S;
      S4 = S2 * S2;
      S5 = S4 * S;
      CS = C0 + C1 * S + C2 * S2 + C3 * S3 + C4 * S4 + C5 * S5;
      CP1 = (C1 + 2.0 * C2 * S + 3.0 * C3 * S2 + 4.0 * C4 * S3 + 5.0 * C5 * S4) / D;
      CP2 = (2.0 * C2 + 6.0 * C3 * S + 12.0 * C4 * S2 + 20.0 * C5 * S3) / (D * D);
      CP3 = (6.0 * C3 + 24.0 * C4 * S + 60.0 * C5 * S2) / (D * D * D);
      CP4 = (24.0 * C4 + 120.0 * C5 * S) / (D * D * D * D);

      CP5 = 120.0 * C5 / (D * D * D * D * D);
      
      std::cout << "ZZZZZZZZZZZZZZZZZZZZZZ\n";
      std::cout << "S " << S << std::endl;
      std::cout << "D " << D << std::endl;
      std::cout << "S2 " << S2 << std::endl;
      std::cout << "S3 " << S3 << std::endl;
      std::cout << "S4 " << S4 << std::endl;
      std::cout << "S5 " << S5 << std::endl;
      std::cout << "CS " << CS << std::endl;
      std::cout << "CP1 " << CP1 << std::endl;
      std::cout << "CP2 " << CP2 << std::endl;
      std::cout << "CP3 " << CP3 << std::endl;
      std::cout << "CP4 " << CP4 << std::endl;
      std::cout << "CP5 " << CP5 << std::endl;
      std::cout << "ZZZZZZZZZZZZZZZZZZZZZZ\n";
      

      if (std::abs(CS) > 70.0) {CS = std::copysign(70.0, CS);}
      E = std::exp(CS);
      RE = 1.0 / (1.0 + E);
      ERE = E * RE;
      ERE1 = ERE * RE;
      ERE2 = ERE * ERE1;
      ERE3 = ERE * ERE2;
      ERE4 = ERE * ERE3;
      
      ERE5 = ERE * ERE4;
      ERE6 = ERE * ERE5;

      CP12 = CP1 * CP1;
      CP13 = CP1 * CP12;
      CP14 = CP12 * CP12;
      CP22 = CP2 * CP2;
      
      CP15 = CP12 * CP13;
      CP16 = CP13 * CP13;
      CP23 = CP2 * CP22;
      CP32 = CP3 * CP3;
      
      std::cout << "ZZZZZZZZZZZZZZZZZZZZZZ\n";
      std::cout << "ERE " << ERE << std::endl;
      std::cout << "ERE1 " << ERE1 << std::endl;
      std::cout << "ERE2 " << ERE2 << std::endl;
      std::cout << "ERE3 " << ERE3 << std::endl;
      std::cout << "ERE4 " << ERE4 << std::endl;
      std::cout << "ERE5 " << ERE5 << std::endl;
      std::cout << "ERE6 " << ERE6 << std::endl;
      std::cout << "CP12 " << CP12 << std::endl;
      std::cout << "CP13 " << CP13 << std::endl;
      std::cout << "CP14 " << CP14 << std::endl;
      std::cout << "CP15 " << CP15 << std::endl;
      std::cout << "CP16 " << CP16 << std::endl;
      std::cout << "CP22 " << CP22 << std::endl;
      std::cout << "CP23 " << CP23 << std::endl;
      std::cout << "CP32 " << CP32 << std::endl;
      std::cout << "ZZZZZZZZZZZZZZZZZZZZZZ\n";

      if (IGP == 6) return;
      G1 = -CP1 * ERE1;

      if (IGP == 5) return;
      G2 = -(CP2 + CP12) * ERE1 + 2.0 * CP12 * ERE2;
      if (IGP == 4) return;
      G3 = -(CP3 + 3.0 * CP1 * CP2 + CP13) * ERE1 +
         6.0 * (CP1 * CP2 + CP13) * ERE2 - 6.0 * CP13 * ERE3;

      if (IGP == 3) return;
      G4 = -(CP4 + 4.0 * CP1 * CP3 + 3.0 * CP22 + 6.0 * CP12 * CP2 + CP14) * ERE1 +
         (8.0 * CP1 * CP3 + 36.0 * CP12 * CP2 + 6.0 * CP22 + 14.0 * CP14) * ERE2 -
         36.0 * (CP12 * CP2 + CP14) * ERE3 + 24.0 * CP14 * ERE4;

      if (IGP != 2) return;
      G5 = (-CP5 - 5.0 * CP1 * CP4 - 10.0 * CP2 * CP3 - 10.0 * CP12 * CP3 -
           15.*CP1*CP22 - 10.*CP13*CP2 - CP15)*ERE1 +
           (10.*CP1*CP4 + 20.*CP2*CP3 + 60.*CP12*CP3 + 90.*CP1*CP22 + 
           140.*CP13*CP2 + 30.*CP15)*ERE2 + (-60.*CP12*CP3 - 
           90.*CP1*CP22 - 360.*CP13*CP2 - 150.*CP15)*ERE3 + 
           (240.*CP13*CP2 + 240.*CP15)*ERE4 + (-120.*CP15)*ERE5;
      G6 = (-6.*CP1*CP5 - 15.*CP2*CP4 - 15.*CP12*CP4 - 10.*CP32 - 
           60.*CP1*CP2*CP3 - 20.*CP13*CP3 - 15.*CP23 - 45.*CP12*CP22 -
           15.*CP14*CP2 - CP16)*ERE1 + (12.*CP1*CP5 + 30.*CP2*CP4 + 
           90.* CP12*CP4 + 20.*CP32 + 360.*CP1*CP2*CP3 + 280.*CP13*CP3 +
           90.* CP23 + 630.*CP12*CP22 + 450.*CP14*CP2 + 62.*CP16)*ERE2 +
           (-90.*CP12*CP4 - 360.*CP1*CP2*CP3 - 720.*CP13*CP3 - 90.*CP23 - 
           1620.*CP12*CP22 - 2250.*CP14*CP2 - 540.*CP16)*ERE3 +
           (480.*CP13*CP3 + 1080.*CP12*CP22 + 3600.*CP14*CP2 +
           1560.* CP16)*ERE4 + (-1800.*CP14*CP2 - 1800.*CP16)*ERE5 +
           720.*CP16*ERE6;
       }





}
