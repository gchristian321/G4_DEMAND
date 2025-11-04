#include "DRAGONEMField.hh"


//C***********************************************************************
//C                         SOLENOID SUBROUTINES
//C***********************************************************************
//C

namespace DRAGON {

void DRAGONEMField::mitray_solnd(G4double* DATA, G4double* XPOS, G4double* BFLD) const
     {
//C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     Subroutine for solenoid, in GEANT implementation of MIT-RAYTRACE C
//C     adapted from:                                                    C
//C     Subroutine SOLND (NO, NP, T, TP ,NUM ) by S. Kowalski            C
//C     by Stanley Yen (TRIUMF)                                          C
//C                                                                      C
//C      TC(1) to  TC(6) =  (  X,  Y,  Z, VX, VY, VZ )                   C
//C     DTC(1) to DTC(6) =  ( VX, VY, VZ, VXDOT, VYDOT, VZDOT )          C
//C                                                                      C
//C     BF (positive) : Solenoid field in beam direction                 C
//C                                                                      C
//C     MODIFICATION HISTORY                                             C
//C     --------------------                                             C
//C                                                                      C
//C     Mar 16, 1998    S.Yen     Original adaptation from MIT-RAYTRACE  C
//C                                                                      C
//C     S. Yen (TRIUMF)  e-mail  STAN@TRIUMF.CA                          C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C

      G4double XA, YA, ZA;
      G4double A, B, L, D, Z11, Z22;
      G4double S, BT;
//C
//C     Extract the A-axes coordinates from XPOS array
//C
      XA = XPOS[0];
      YA = XPOS[1];
      ZA = XPOS[2];
//C
//C     EXTRACT THE SOLENOID PARAMETERS FROM THE ARRAY DATA(i)
//C
      A   = DATA[9];
      B   = DATA[10];
      L   = DATA[11];
      D   = DATA[12];
      BF  = DATA[13];
      Z11 = DATA[14];
      Z22 = DATA[15];
//C
      AL  = L/2.;
      RAD = D/2.;
      BX  = 0.;
      BY  = 0.;
      BZ  = 0.;
      BT  = 0.;
      S   = 0.;
//C
//C     Load position coordinates into TC(1),TC(2),TC(3) since this is
//C     what subroutine BSOL expects.  Note that these are relative to
//C     the geometric center of the solenoid.
//C
      TC[0] =  XA;
      TC[1] =  YA;
      TC[2] =  ZA - A - AL;
//C
//C     If outside the integration regions defined by Z11 before the 
//C     entrance edge of the solenoid, to Z22 beyond the exit edge of 
//C     the solenoid, then return without calculating anything

      if(TC[2] < (-(AL+Z11)) || TC[2] > (AL+Z22)) 
        {return;}
//C
//C     CALL BSOL TO CALCULATE B-FIELD COMPONENTS BX,BY,BZ, 
//C     RETURNED IN COMMON BLOCK
//C
      mitray_bsol();
//C
      BFLD[0] = BX;
      BFLD[1] = BY;
      BFLD[2] = BZ;
//C
       return;
      }
   
//C
//C=======================================================================
//C
void DRAGONEMField::mitray_bsol() const
     {
//C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     Routine valid for fields outside central zone of elemental       C
//C                              Solenoid                                C
//C                                                                      C
//C     BF    = field at center of infinite solenoid; curr. den. (NI/M)  C
//C                                                                      C
//C     M.W.GARRETTT  JOURNAL OF APP. PHYS. 34,(1963),P2567              C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C

      const G4double PI4 = 12.566370616;
      G4double X, Y, Z; 
      G4double R, COSA, COSB;
      G4double RADR,AAPR,AAMR,RCSQ;
      G4double ZZ, R1SQ, R1, RKSQ;
      G4double BZS1,BRS1,BZS2,BRS2;
      G4double VKS, VES, P;

      X = TC[0];
      Y = TC[1];
      Z = TC[2];

      R = std::sqrt( X*X + Y*Y );

      if( R  >  (RAD/1.0e4)  ) 
        {
         COSA = (AL-Z) / std::sqrt( RAD*RAD + std::pow((AL-Z),2));
         COSB =-(AL+Z) / std::sqrt( RAD*RAD + std::pow((AL+Z),2));

         BX = 0.;
         BY = 0.;
         BZ = BF*(COSA-COSB)/2.0;
         BT = std::fabs(BZ);
         
         return;
         }

      RADR = RAD+R;
      AAPR = 4.0*RAD/RADR;
      AAMR = (RAD-R)/(2.0*RAD);
      RCSQ = 4.0*RAD*R/(RADR*RADR);
//C
//C *** Solenoid left hand source
//C
      ZZ   = -(AL+Z);
      R1SQ = RADR*RADR + ZZ*ZZ;
      R1   = std::sqrt(R1SQ);
      RKSQ = 4.0*RAD*R/R1SQ;

      mitray_FB01AD(RKSQ, VKS, VES );
      mitray_FB03AD(RCSQ, RKSQ, P );

      BZS1 = AAPR*ZZ*(VKS+AAMR*(P-VKS) ) /R1;
      BRS1 = R1*(2.0*(VKS-VES) - RKSQ*VKS);
//C
//C *** Solenoid right hand source
//C
      ZZ   = AL-Z;
      R1SQ = RADR*RADR + ZZ*ZZ;
      R1   = std::sqrt(R1SQ);
      RKSQ = 4.0*RAD*R/R1SQ;

      mitray_FB01AD(RKSQ, VKS, VES );
      mitray_FB03AD(RCSQ, RKSQ, P );

      BZS2 = AAPR*ZZ*(VKS+AAMR*(P-VKS) ) /R1;
      BRS2 = R1*(2.0*(VKS-VES) - RKSQ*VKS);

      BZ = BF*( BZS2-BZS1 )/PI4;
      BR = BF*( BRS2-BRS1 )/(R*PI4);
      BX = BR * X /R;
      BY = BR *  Y/R;
      BT = std::sqrt(std::pow(BX,2) + std::pow(BY,2) + std::pow(BZ,2) );

      return;
      }
//C
//C=======================================================================
//C   

void DRAGONEMField::mitray_FB01AD(G4double C, G4double& VK, G4double& VE) const
     {
      const G4double XLG = 1.0e30;  // '7FFFFFFFFFFF'X 
      G4double D, E;

      D = 1.0 - C;
      if (D > 0.0) 
         {E = -std::log(D);}
//C
//C *** Harwell version of FB01AD
//C
      if (C >= 1.0) 
         {
          VE = 1.0;
          VK = XLG;
          
          return;
          }

      VE = E * (((((((((( 
          3.18591956555015718e-5 * D   + 0.989833284622538479e-3) * D  +
          0.643214658643830177e-2) * D + 0.16804023346363385e-1)  * D  +
          0.261450147003138789e-1) * D + 0.334789436657616262e-1) * D  +
          0.427178905473830956e-1) * D + 0.585936612555314917e-1) * D  +
          0.937499997212031407e-1) * D + 0.249999999999901772e0)  * D) +
          (((((((((
          0.149466217571813268e-3 * D  + 0.246850333046072273e-2) * D  +
          0.863844217360407443e-2) * D + 0.107706350398664555e-1) * D  +
          0.782040406095955417e-2) * D + 0.759509342255943228e-2) * D  +
          0.115695957452954022e-1) * D + 0.218318116761304816e-1) * D  +
          0.568051945675591566e-1) * D + 0.443147180560889526)    * D  +
          1.0;
//C
//C *** Routine modified to calculate VD and VE always
//C
      VK = E * (((((((((( 
          0.297002809665556121e-4  * D + 0.921554634963249846e-3) * D +
          0.597390429915542916e-2) * D + 0.155309416319772039e-1) * D +
          0.239319133231107901e-1) * D + 0.301248490128989303e-1) * D +
          0.373777397586236041e-1) * D + 0.48828041906862398e-1)  * D +
          0.703124997390383521e-1) * D + 0.124999999999908081e0)  * D + 
          0.5) +
          (((((((((( 
          0.139308785700664673e-3  * D + 0.229663489839695869e-2) * D +
          0.800300398064998537e-2) * D + 0.984892932217689377e-2) * D +
          0.684790928262450512e-2) * D + 0.617962744605331761e-2) * D +
          0.878980187455506468e-2) * D + 0.149380135326871652e-1) * D +
          0.308851462713051899e-1) * D + 0.965735902808562554e-1) * D +
          1.38629436111989062e0);
      }  
//C
//C=======================================================================
//C  

void DRAGONEMField::mitray_FB02AD(G4double CAYSQ, G4double SINP, G4double COSP, G4double& E, G4double& F) const
     {
      G4double PHI, H, A, N; 
      G4double SIG1, SIG2, SIG3, SIG4, SIN2, TERM, CRIT;
      G4double RECIP, FACT, H1, DEL1, DEL2, DEL3, DEL4; 
      G4double FACT1, FACTOR, CAYDSQ, FACTRO, FLOG1;
      G4double FACTN, FACTM;
      G4double T1, T2, CFI, CFJ, CFL, CFM, CFN, CFI1, CFJ1;

      PHI = std::atan(SINP / COSP);  

      if (CAYSQ * SINP * SINP - 0.5 >= 0) 
         {
          H = 1.0;
          A = PHI;
          N = 0;
          SIG1 = 0.0;
          SIG2 = 0.0;
          SIN2 = SINP * SINP;
          TERM = SINP * COSP * 0.5;
          CRIT = PHI;

          do {
              N = N + 1;
              RECIP = 1.0 / N;
              FACT = (N - 0.5) * RECIP;
              H1 = H;
              H = FACT * CAYSQ * H;
              A = FACT * A - TERM * RECIP;
              TERM = TERM * SIN2;
              CRIT = CRIT * SIN2;
              DEL1 = H * A;
              DEL2 = -0.5 * RECIP * CAYSQ * H1 * A;
              SIG1 = SIG1 + DEL1;
              SIG2 = SIG2 + DEL2;
              
              if (std::fabs(DEL1) - 4.0e-16 > 0) 
                 {
                  F = PHI + SIG1;
                  E = PHI + SIG2;
               
                  return;
                  }
              }
          while(std::fabs(CRIT) - std::fabs(A) <= 0);
          
          F = PHI + SIG1;
          E = PHI + SIG2;
               
          return;
          } 
      else 
         {
          CFI = 1.0;
          CFJ = 1.0;
          CFL = 0.0;
          CFM = 0.0;
          CFN = 0.0;
          SIG1 = 0.0;
          SIG2 = 0.0;
          SIG3 = 0.0;
          SIG4 = 0.0;

          N = 0;
          
          FACT1 = 1.0 - CAYSQ * SINP * SINP;
          FACTOR = 0.5 * COSP * std::sqrt(CAYSQ / FACT1);
          FACTRO = FACTOR + FACTOR;
          CAYDSQ = 1.0 - CAYSQ;

          do {
              N = N + 1;
              RECIP = 1.0 / N;
              FACTN = RECIP * (N - 0.5);
              FACTM = (N + 0.5) / (N + 1.0);
              FACTOR = FACTOR * FACT1;
              CFI1 = CFI;
              CFJ1 = CFJ;
              CFI = CFI * FACTN;
              CFJ = CFJ * FACTN * FACTN * CAYDSQ;
              CFL = CFL + 0.5 / (N * (N - 0.5));
              CFM = (CFM - FACTOR * RECIP * CFI) * FACTM * FACTM * CAYDSQ;
              CFN = (CFN - FACTOR * RECIP * CFI1) * FACTN * FACTM * CAYDSQ;
              DEL1 = CFM - CFJ * CFL;
              DEL2 = CFN - (FACTN * CFL - 0.25 * RECIP * RECIP) * CAYDSQ * CFJ1;
              DEL3 = CFJ;
              DEL4 = FACTM * CFJ;
              SIG1 = SIG1 + DEL1;
              SIG2 = SIG2 + DEL2;
              SIG3 = SIG3 + DEL3;
              SIG4 = SIG4 + DEL4;
              }
          while(fabs(DEL1) - 4.0e-16 <= 0);
               
          G4double CAYMOD = std::sqrt(CAYSQ);
          FLOG1 = std::log(4.0 / (std::sqrt(FACT1) + CAYMOD * COSP));
          T1 = (1.0 + SIG3) * FLOG1 + FACTRO * std::log(0.5 + 0.5 * CAYMOD * std::fabs(SINP));
          T2 = (0.5 + SIG4) * CAYDSQ * FLOG1 + 1.0 - FACTRO * (1.0 - CAYMOD * std::fabs(SINP));
          F = T1 + SIG1;
          E = T2 + SIG2;
          
          return;
          }
}
   
//C
//C=======================================================================
//C  

void DRAGONEMField::mitray_FB03AD(G4double GN, G4double CACA, G4double& P) const
     {
      G4double STH, CTH, CADA, CAPK, CAPE, E, F, PI;

      if (GN > 0) 
         {
          if (CACA >= 0) 
             {
              P = 1.5707963268 / std::sqrt(1.0 - GN);
              
              return;
              } 
          else 
             {
              STH = std::sqrt(-GN / (CACA - GN));
              CTH = std::sqrt(1.0 - STH * STH);
              CADA = 1.0 - CACA;

              mitray_FB01AD(CACA, CAPK, CAPE);
              mitray_FB02AD(CADA, STH, CTH, E, F);

              BR = CAPE * F - CAPK * (F - E);
              P = CAPK * CTH * CTH + STH * BR / std::sqrt(1.0 - GN);
              
              return;
              }
          } 
       else 
          {
           if (GN - CACA > 0) 
              {
               STH = std::sqrt(GN / CACA);
               CTH = std::sqrt(1.0 - STH * STH);

               mitray_FB01AD(CACA, CAPK, CAPE);
               mitray_FB02AD(CACA, STH, CTH, E, F);

               BR = CAPK * E - CAPE * F;
               P = CAPK + BR * STH / (CTH * std::sqrt(1.0 - GN));
               
               return;
               }
           else 
              {
               if (GN - CACA == 0)
                  {
                   mitray_FB01AD(CACA, CAPK, CAPE);
                   P = CAPE / (1.0 - CACA);
                   
                   return;
                   } 
               else
                  { 
                   CADA = 1.0 - CACA;
                   PI = 3.1415926536;
                   STH = std::sqrt((1.0 - GN) / CADA);
                   CTH = std::sqrt(1.0 - STH * STH);

                   mitray_FB01AD(CACA, CAPK, CAPE);
                   mitray_FB02AD(CADA, STH, CTH, E, F);

                   BR = PI / 2.0 + CAPK * (F - E) - CAPE * F;
                   P = CAPK + BR * std::sqrt(GN) / (CADA * STH * CTH);
                   
                   return;
                   }
               }
            } 
     }

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
      
 
}
