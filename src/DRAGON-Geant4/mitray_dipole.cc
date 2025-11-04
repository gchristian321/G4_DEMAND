#include "DRAGONEMField.hh"
#include "G4RunManager.hh"


namespace DRAGON {

/*
***********************************************************************
                          DIPOLE SUBROUTINES
************************************************************************
*/


void DRAGONEMField::mitray_dipo_AtoB(G4double XA, G4double YA, G4double ZA, G4double& XB, G4double& YB, G4double& ZB) const
     {
//C
//C	TRANSFORM FROM INITIAL ENTRANCE COORDINATES TO VFB COORD.
//C	I.E. FROM SYSTEM A TO SYSTEM B OF STANDARD LAYOUT, FOR A DIPOLE.
//C	ALPHA IS THE POLEFACE ROTATION ANGLE OF THE ENTRANCE FACE
//C
//C	INPUT:  XA, YA, ZA
//C	OUTPUT: XB, YB, ZB
//C	
      G4double COSA = std::cos(ALPHA/deg_rad);
      G4double SINA = std::sin(ALPHA/deg_rad);

      XB = (A-ZA)*SINA - XA*COSA;
      YB = YA;
      ZB = (A-ZA)*COSA + XA*SINA;
	
//C
//C	COORDINATE SYSTEM B IS DISPLACED BY XCR1 IN THE +XB DIRECTION
//C

      XB = XB - XCR1;

//C
//C	ALL COORDINATES ARE NOW IN TERMS OF COORDINATE SYSTEM B
//C
	
      return; 
      }

void DRAGONEMField::mitray_dipo_BtoC(G4double XB, G4double YB, G4double ZB, G4double& XC, 
    G4double& YC, G4double& ZC) const
    {
     //C
     //C	TRANFORM FROM SYSTEM B TO SYSTEM C COORDINATES OF A DIPOLE
     //C
	 
     G4double COPAB = std::cos( (PHI-ALPHA-BETA)/deg_rad );
     G4double SIPAB = std::sin( (PHI-ALPHA-BETA)/deg_rad );
     G4double COSPB = std::cos( (PHI/2. - BETA)/deg_rad );
     G4double SINPB = std::sin( (PHI/2. - BETA)/deg_rad );
     G4double SIP2  = std::sin( (PHI/2.)/deg_rad );
//C
//C	THE ORIGIN OF THE B SYSTEM IS DISPLACED FROM THE CENTRAL RAY
//C	BY AMOUNT XCR1 ALONG THE +XB AXIS
//C

     G4double XT = XB;
     G4double ZT = ZB;

//C
//C	TRANSFORM BACK TO B COORDINATE SYSTEM CENTERED ON CENTRAL RAY
//C
     XT = XT + XCR1;

//C
//C	NOW ROTATE/TRANSLATE TO GET COORDINATES IN TERMS OF C SYSTEM.
//C

     ZC = -ZT  *COPAB + XT  *SIPAB -2.*RB*SIP2*COSPB;
     XC = -ZT  *SIPAB - XT  *COPAB -2.*RB*SIP2*SINPB;

//C
//C	SHIFT C COORDINATE AXES FROM CENTRAL RAY BY AMOUNT XCR2 ALONG
//C	+XC AXIS
//C

     XC = XC - XCR2;
     YC = YB;

//C
//C	COORDINATES ARE NOW IN TERMS OF THE NEW, SHIFTED C SYSTEM.
//C
     return;
}

void DRAGONEMField::mitray_dipo_BBtoBA(G4double BXB, G4double BYB, G4double BZB, G4double& BXA,
    G4double& BYA, G4double& BZA) const
    {
//C
//C	SUBROUTINE TO TRANSFORM B FIELD FROM A TO B COORD SYSTEMS 
//C       OF A DIPOLE, I.E. FROM SYSTEM A TO SYSTEM B OF STANDARD LAYOUT.
//C	ALPHA IS THE POLEFACE ROTATION ANGLE OF THE ENTRANCE FACE
//C
//C	INPUT:  BXB, BYB, BZB 
//C	OUTPUT: BXA, BYA, BZA
//C

     G4double COSA = std::cos(-ALPHA/deg_rad);
     G4double SINA = std::sin(-ALPHA/deg_rad);

     BXA = -BZB*SINA - BXB*COSA;
     BYA = BYB;
     BZA = -BZB*COSA + BXB*SINA;

//C
//C	COORDINATE SYSTEM B IS DISPLACED BY XCR1 IN THE +XB DIRECTION
//C	ALL B-FIELD COMPONENTS ARE NOW IN TERMS OF COORDINATE SYSTEM A
//C

     return;
     }

void DRAGONEMField::mitray_dipo_BCtoBB(G4double BXC, G4double BYC, G4double BZC, G4double& BXB, 
    G4double& BYB, G4double& BZB) const
     {
//C
//C	TRANFORM B-FIELD COMPONENTS FROM SYSTEM C TO SYSTEM B COORDINATES
//C       OF A DIPOLE.
//C
//C	INPUT:  BXC, BYC, BZC   B-FIELD IN C-AXIS SYSTEM
//C	OUTPUT: BXB, BYB, BZB   B-FIELD IN B-AXIS SYSTEM
//C

      G4double COPAB = std::cos(-(PHI-ALPHA-BETA)/deg_rad);
      G4double SIPAB = std::sin(-(PHI-ALPHA-BETA)/deg_rad);
      G4double COSPB = std::cos(-(PHI/2.-BETA)/deg_rad);
      G4double SINPB = std::sin(-(PHI/2.-BETA)/deg_rad);
      G4double SIP2 =  std::sin(-(PHI/2.)/deg_rad);

//C
//C	NOW ROTATE TO GET COORDINATES IN TERMS OF C SYSTEM.
//C
      BZB = -BZC  *COPAB + BXC  *SIPAB;
      BXB = -BZC  *SIPAB - BXC  *COPAB;
      BYB = BYC;

//C
//C	B-FIELD COMPONENTS ARE NOW IN TERMS OF THE B SYSTEM.
//C
      return;
}


void DRAGONEMField::mitray_dipole(G4double* DATA, G4double* XPOS, G4double* BFLD) const
     {
//C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     Subroutine for dipole, in GEANT implementation of MIT-RAYTRACE   C
//C     adapted from:                                                    C
//C     Subroutine DIPOLE ( NO, NP, T, TP ,NUM ) by S. Kowalski          C
//C                                                                      C
//C      TC(1) to  TC(6) =  (  X,  Y,  Z, VX, VY, VZ )                   C
//C     DTC(1) to DTC(6) =  ( VX, VY, VZ, VXDOT, VYDOT, VZDOT )          C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C
       G4double LF1, LF2, LU1;
       G4String REGION;

//C
//C     SY  Addition Feb 10/98
//C     Extract coordinates of the field point,
//C     in the A-axis system of the device.
//C

       G4double XA=XPOS[0];
       G4double YA=XPOS[1];
       G4double ZA=XPOS[2];
	
       G4double BXA, BYA, BZA;
       G4double BXB, BYB, BZB;
       G4double BXC, BYC, BZC;
	
       REGION=" ";
       if (ldiag) 
          {
	   std::cout << std::endl;
	   std::cout << std::setw(50) << std::setfill('-') << "" << std::endl;
  	   std::cout << " ENTER SUBROUTINE DIPOLE" << std::endl;
  	   std::cout << " XA,YA,ZA=";
	   std::cout << std::fixed << std::setprecision(3);
	   std::cout << std::setw(10) << XA << " ";
	   std::cout << std::setw(10) << YA << " ";
	   std::cout << std::setw(10) << ZA << std::endl;
	   }
//C
//C     EXTRACT THE PARAMETERS FOR THE DIPOLE MAGNET
//C
       G4double IR = 0;
       LF1  = DATA[0];
       LU1  = DATA[1];
       LF2  = DATA[2];
       DG   = DATA[3];
       MTYP = DATA[4];
       IMAP = DATA[5];
       A    = DATA[10]; 
       B    = DATA[11];
       D    = DATA[12];
       RB   = DATA[13];
       BF   = DATA[14];
       PHI  = DATA[15];
       ALPHA= DATA[16];
       BETA = DATA[17];
       NDX  = DATA[18];
       BET1 = DATA[19];
       GAMA = DATA[20];
       DELT = DATA[21];
       G4double Z11  = DATA[24];
       G4double Z12  = DATA[25];
       G4double Z21  = DATA[26];
       G4double Z22  = DATA[27];
       G4double BR1  = DATA[40];
       G4double BR2  = DATA[41];
       XCR1 = DATA[42];
       XCR2 = DATA[43];
       G4double WDE  = DATA[48];
       G4double WDX  = DATA[49];
//C
//C    SY  Addition Feb 10/98
//C    These parameters define the bounds of the entrance fringe
//C    field, in the XB axis system.  The entrance fringe field
//C    is defined by XBMIN < XB < XBMAX  .AND.  Z12 < ZB < Z11
//C
       G4double XBMAX = WDE/2.;
       G4double XBMIN = -XBMAX;
//C
//C    Similarly, the bounds of the exit fringe field region are
//C    defined in the XC axis system.  The exit fringe field is 
//C    defined by XCMIN < XC < XCMAX .AND. Z21 < ZC < Z22
//C
       G4double XCMAX = WDX/2.;
       G4double XCMIN = -XCMAX;

       if(MTYP == 0) 
          MTYP = 1;
//C
//C     SY  Addition Feb 10/98
//C     We first zero the B-field components in case we need to abort
//C
       for(G4int i=0; i<3; i++) 
          {BFLD[i] = 0.;}

       BX = 0.;
       BY = 0.;
       BZ = 0.;
       BT = 0.;
       S = 0.;
       BR = BR1;

//C
//C    Transform from the A to the B coordinate system
       mitray_dipo_AtoB(XA,YA,ZA,XB,YB,ZB);

//C    Transform from the B to the C coordinate system
       mitray_dipo_BtoC(XB,YB,ZB,XC,YC,ZC);
//C    Print out coordinates if diagnostic mode
//C
       if(ldiag) 
         {
	  std::cout << " XB,YB,ZB= " << std::setw(10) << XB << std::setw(10) << YB 
	 	    << std::setw(10) << ZB << std::endl;

		// Output XC, YC, ZC
	  std::cout << " XC,YC,ZC= " << std::setw(10) << XC << std::setw(10) << YC 
		    << std::setw(10) << ZC << std::endl;

		// Output Z11, Z12, Z21, Z22, XBMIN, XBMAX, XCMIN, XCMAX
	  std::cout << " Z11= " << std::setw(10) << Z11 << "  Z12= " << std::setw(10) 
		    << Z12 << std::endl;
	  std::cout << " Z21= " << std::setw(10) << Z21 << "  Z22= " << std::setw(10) 
		    << Z22 << std::endl;
	  std::cout << " XBMIN= " << std::setw(10) << XBMIN << "  XBMAX= " << std::setw(10) 
		    << XBMAX << std::endl;
	  std::cout << " XCMIN= " << std::setw(10) << XCMIN << "  XCMAX= " << std::setw(10) 
		    << XCMAX << std::endl;
	  }
//C
//C   Now determine what region we are in -- 
//C   Before start of the entrance fringe field ("entrance far field", IN=-99 )
//C   entrance fringe field (IN=1), 
//C   "uniform" field region (IN=2), 
//C   exit fringe field region (IN=3),
//C   after the end of the exit fringe field ("exit far field", IN=+99)
//C   Because of the possibility of confusing entrance and exit regions, we
//C   first check for entrance/exit fringe field, then uniform field, then
//C   entrance/exit far field regions, to be sure that we get the most
//C   important regions first
//C

      if(ZB <= Z11 && ZB > Z12 && XB >= XBMIN &&XB <= XBMAX) 
        {
std::cout << "CASO 1" << std::endl;
//C
//C        *************************
//C        *                       *
//C        * ENTRANCE FRINGE FIELD *
//C        *                       *
//C        *************************
//C
//C        Entrance fringe field region, B-axis coordinates are used.
//C
         IR = 1;
         IN = 1;
         XC_OFFSET = RB*std::cos( ALPHA/deg_rad );
         ZC_OFFSET =-RB*std::sin( ALPHA/deg_rad );
//C
//C      Load the B-axis coordinates into TC(1),TC(2),TC(3), because this
//C      is where subroutine BDIP expects to find the coordinates.
//C

         TC[0] = XB;
         TC[1] = YB;
         TC[2] = ZB;
         
         C0   = DATA[28]; 
	 C1   = DATA[29];
	 C2   = DATA[30];
	 C3   = DATA[31];
	 C4   = DATA[32];
	 C5   = DATA[33];
	 DELS = DATA[44];
	 RCA  = DATA[46];
	 WDIP = DATA[48];
//C
//C      S2...S8 are the coefficients for the entrance face curvature
//C      SY Feb 20/98   Calculate powers of RB, use DATA(57)/RB4/RB3 
//C      instead of RB**7 to avoid exponent overflow in case of large RB
//C
         G4double RB2 = RB*RB;
         G4double RB3 = RB2*RB;
         G4double RB4 = RB3*RB;
         S2 = DATA[50] / RB    + RCA / 2.0;
	 S3 = DATA[51] / (RB * RB);
	 S4 = DATA[52] / (RB * RB * RB) + pow(RCA, 3) / 8.0;
	 S5 = DATA[53] / pow(RB, 4);
	 S6 = DATA[54] / (RB * RB * RB * RB) + pow(RCA, 5) / 16.0;
	 S7 = DATA[55] / pow(RB, 6);
	 S8 = DATA[56] / pow(RB, 7) + pow(RCA, 7) / 25.6;
//C
//C        CHECK IF WE HAVE A FLAT BOUNDARY
//C                 NSRF=0 FLAT
//C                     =1 CURVED
//C
std::cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV" << std::endl; 
std::cout << "S2 " << S2 << std::endl;
std::cout << "S3 " << S3 << std::endl;
std::cout << "S4 " << S4 << std::endl;
std::cout << "S5 " << S5 << std::endl;
std::cout << "S6 " << S6 << std::endl;
std::cout << "S7 " << S7 << std::endl;
std::cout << "S8 " << S8 << std::endl;

	 NSRF = 1;
	 if ((S2 == 0.0) && (S3 == 0.0) && (S4 == 0.0) && (S5 == 0.0) && 
	    (S6 == 0.0) && (S7 == 0.0) && (S8 == 0.0)) 
	    {NSRF = 0;}
//C
//C        Call BDIP to calculate the B-field components
//C
	 mitray_bdip();                     
//C
//C	 BX, BY, BZ ARE IN B-AXIS SYSTEM;  TRANSFORM TO A-AXIS SYSTEM
//C
 	 mitray_dipo_BBtoBA(BX, BY, BZ, BXA, BYA, BZA);

std::cout << "BX " << BX << std::endl;
std::cout << "BY " << BY << std::endl;
std::cout << "BZ " << BZ << std::endl;
std::cout << "BXA " << BXA << std::endl;
std::cout << "BYA " << BYA << std::endl;
std::cout << "BZA " << BZA << std::endl;
std::cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV" << std::endl;

	 if (ldiag) 
	    {
	     std::cout << "ENTRANCE FRINGE FIELD REGION" << std::endl;
	     std::cout << "B SYSTEM BX,BY,BZ=" << std::setw(10) << BX << std::setw(10) << BY << std::setw(10) << BZ << std::endl;
	     std::cout << "A SYSTEM BX,BY,BZ=" << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
	     }

	 BFLD[0] = BXA;
	 BFLD[1] = BYA;
	 BFLD[2] = BZA;

std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;

	 return;
	 }

//C
//C        -------------------------
//C
       else if (ZC > Z21 && ZC <= Z22 && XC >= XCMIN && XC <= XCMAX) 
               {
std::cout << "CASO 2" << std::endl;
//C
//C        *********************
//C        *                   *
//C        * EXIT FRINGE FIELD *
//C        *                   *
//C        *********************
//C
//C        Exit fringe field region, C-axis coordinates are used.
//C        Setup for second fringe field and integration
//C        IN=3 designates exit fringe field

		 IN = 3;
		 IR = 2;
		 XC_OFFSET = -RB * std::cos(BETA / deg_rad);
		 ZC_OFFSET = -RB * std::sin(BETA / deg_rad);
//C
//C        Load the C axis coordinates into array TC, because this is
//C        where subroutine BDIP expects to find the coordinates.
//C
                 TC[0] = XC;
		 TC[1] = YC;
		 TC[2] = ZC;

		 BR = BR2;
//C        C0,...,C5 are the expansion coefficients for the exit fringe field
		 C0 = DATA[34];  
		 C1 = DATA[35];
		 C2 = DATA[36];
		 C3 = DATA[37];
		 C4 = DATA[38];
		 C5 = DATA[39];
		 DELS = DATA[45];
//C        RCA is inverse radius of curvature of exit boundary
		 RCA = DATA[47];
		 WDIP = DATA[49];
//C
//C        S2...S8 are expansion coefficients for shape of exit
//C        face of dipole
//C        SY Feb 20, 1998  Calculate powers of RB, use DATA(64)/RB4/RB3 
//C        to avoid calculating RB**7, which could suffer exponent overflow
//C        in case of very large RB value, as in clamshell dipole.
//C
		 G4double RB2 = RB * RB;
		 G4double RB3 = RB2 * RB;
		 G4double RB4 = RB3 * RB;

		 S2 = DATA[57] / RB + RCA / 2.0;
		 S3 = DATA[58] / RB2;
		 S4 = DATA[59] / RB3 + std::pow(RCA, 3) / 8.0;
		 S5 = DATA[60] / RB4;
		 S6 = DATA[61] / (RB3 * RB2) + std::pow(RCA, 5) / 16.0;
		 S7 = DATA[62] / (RB3 * RB3);
		 S8 = DATA[63] / (RB4 * RB3) + std::pow(RCA, 7) / 25.6;
//C
//C        CHECK IF WE HAVE A FLAT BOUNDARY
//C                 NSRF=0 FLAT
//C                     =1 CURVED
//C
		 NSRF = 1;
		 if (S2 == 0.0 && S3 == 0.0 && S4 == 0.0 && S5 == 0.0 && S6 == 0.0 && S7 == 0.0 && S8 == 0.0) 
		    {NSRF = 0;}
//C
//C        Call BDIP to calculate magnetic field components
std::cout << "////////////////////////" << std::endl;
std::cout << "XC_OFFSET " << XC_OFFSET << std::endl;
std::cout << "ZC_OFFSET " << ZC_OFFSET << std::endl;
std::cout << "BR " << BR << std::endl;
std::cout << "S2 " << S2 << std::endl;
std::cout << "S3 " << S3 << std::endl;
std::cout << "S4 " << S4 << std::endl;
std::cout << "S5 " << S5 << std::endl;
std::cout << "S6 " << S6 << std::endl;
std::cout << "S7 " << S7 << std::endl;
std::cout << "S8 " << S8 << std::endl;
std::cout << "NSRF " << NSRF << std::endl;
std::cout << "TC[0] " << TC[0] << std::endl;
std::cout << "TC[1] " << TC[1] << std::endl;
std::cout << "TC[2] " << TC[2] << std::endl;
std::cout << "C0 " << C0 << std::endl;
std::cout << "C1 " << C1 << std::endl;
std::cout << "C2 " << C2 << std::endl;
std::cout << "C3 " << C3 << std::endl;
std::cout << "C4 " << C4 << std::endl;
std::cout << "C5 " << C5 << std::endl;
std::cout << "DELS " << DELS << std::endl;
std::cout << "RCA " << RCA << std::endl;
std::cout << "WDIP " << WDIP << std::endl;
//C
		 mitray_bdip();
//C
//C        BX, BY, BZ ARE IN C-AXIS SYSTEM;  FIRST TRANSFORM TO B-AXIS SYSTEM
//C
	 	 mitray_dipo_BCtoBB(BX, BY, BZ, BXB, BYB, BZB);
//C
//C        THEN TRANSFORM FROM B TO A-AXIS SYSTEM
//C
		 G4double BXA, BYA, BZA;
		
		 mitray_dipo_BBtoBA(BXB, BYB, BZB, BXA, BYA, BZA);

		 if (ldiag) 
		    {
		     std::cout << "EXIT FRINGE FIELD REGION" << std::endl;
		     std::cout << "C SYSTEM BX,BY,BZ= " << std::setw(10) << BX << std::setw(10) << BY << std::setw(10) << BZ << std::endl;
		     std::cout << "B SYSTEM BXB,BYB,BZB= " << std::setw(10) << BXB << std::setw(10) << BYB << std::setw(10) << BZB << std::endl;
		     std::cout << "A SYSTEM BXA,BYA,BZA= " << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
		     }

		 BFLD[0] = BXA;
		 BFLD[1] = BYA;
		 BFLD[2] = BZA;

std::cout << "BXA " << BXA << std::endl;
std::cout << "BYA " << BYA << std::endl;
std::cout << "BZA " << BZA << std::endl;
std::cout << "////////////////////////" << std::endl;	

std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;	
		
		 return;
	         }
//C
//C        -------------------------
//C
               else if (ZB <= Z12 && ZC <= Z21) 
                       {
std::cout << "CASO 3" << std::endl;                       
//C
//C        ************************
//C        *                      *
//C        * UNIFORM FIELD REGION *
//C        *                      *
//C        ************************
//C    
//C       UNIFORM FIELD REGION;  C-AXIS COORDINATES ARE USED
//C
		        S = 0.0;
		        IN = 2;
		        XC_OFFSET = -RB * cos(BETA / deg_rad);
		        ZC_OFFSET = -RB * sin(BETA / deg_rad);
//C
//C        Load the C axis coordinates into array TC, because this is where
//C        subroutine BDIP expects to find the coordinates.
		        TC[0] = XC; 
		        TC[1] = YC;
		        TC[2] = ZC;  

		        DELS = 0.0;

		        mitray_bdip();
//C
//C        BX, BY, BZ ARE IN C-AXIS SYSTEM.  TRANSFORM FIRST TO B-AXIS SYSTEM
//C
		        mitray_dipo_BCtoBB(BX, BY, BZ, BXB, BYB, BZB);
//C
//C        THEN TRANSFORM FIELD FROM B TO A-AXIS SYSTEM
//C
		       mitray_dipo_BBtoBA(BXB, BYB, BZB, BXA, BYA, BZA);

	 	       if (ldiag) 
	 	          {
			   std::cout << "UNIFORM FIELD REGION" << std::endl;
			   std::cout << "C SYSTEM BX,BY,BZ= " << std::setw(10) << BX << std::setw(10) << BY << std::setw(10) << BZ << std::endl;
			   std::cout << "B SYSTEM BXB,BYB,BZB= " << std::setw(10) << BXB << std::setw(10) << BYB << std::setw(10) << BZB << std::endl;
			   std::cout << "A SYSTEM BXA,BYA,BZA= " << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
		           }

		       BFLD[0] = BXA;  
		       BFLD[1] = BYA;  
	 	       BFLD[2] = BZA;  
	 	       
std::cout << "BXA " << BXA << std::endl;
std::cout << "BYA " << BYA << std::endl;
std::cout << "BZA " << BZA << std::endl;

std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;

		       return;
	               }
//C
//C        -------------------------
//C
	             else if (ZB > Z11 && XB >= XBMIN && XB <= XBMAX) 
	                     {
std::cout << "CASO 4" << std::endl;
//C
//C        **********************
//C        *                    *
//C        * ENTRANCE FAR FIELD *
//C        *                    *
//C        **********************
//C
		              IN = -99;
		              BXC = 0.0;
		              BYC = 0.0;
	 	              BZC = 0.0;
		              BXB = 0.0;
	 	              BYB = 0.0;
		              BZB = 0.0;
		              BXA = 0.0;
		              BYA = 0.0;
		              BZA = 0.0;

		              if (ldiag) 
		                 {
			          std::cout << "DIPOLE ENTRANCE FAR FIELD REGION" << std::endl;
			          std::cout << "C SYSTEM BXC,BYC,BZC= " << std::setw(10) << BXC << std::setw(10) << BYC << std::setw(10) << BZC << std::endl;
			          std::cout << "B SYSTEM BXB,BYB,BZB= " << std::setw(10) << BXB << std::setw(10) << BYB << std::setw(10) << BZB << std::endl;
			          std::cout << "A SYSTEM BXA,BYA,BZA= " << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
		                  }

		              BFLD[0] = BXA;  
		              BFLD[1] = BYA;  
		              BFLD[2] = BZA;  
		              
std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;
		
		              return;
	                      }
//C
//C        -------------------------
//C
	                    else if (ZC > Z22 && XC > XCMIN && XC <= XCMAX) 
	                            {
std::cout << "CASO 5" << std::endl;
//C
//C        ******************
//C        *                *
//C        * EXIT FAR FIELD *
//C        *                *
//C        ******************
//C
		                     IN = 99;
		                     IN = -99; 
		                     BXC = 0.0;
		                     BYC = 0.0;
		                     BZC = 0.0;
		                     BXB = 0.0;
		                     BYB = 0.0;
		                     BZB = 0.0;
		                     BXA = 0.0;
		                     BYA = 0.0;
		                     BZA = 0.0;

		                     if (ldiag) 
		                        {
			                 std::cout << "DIPOLE EXIT FAR FIELD REGION" << std::endl;
			                 std::cout << "C SYSTEM BXC,BYC,BZC= " << std::setw(10) << BXC << std::setw(10) << BYC << std::setw(10) << BZC << std::endl;
			                 std::cout << "B SYSTEM BXB,BYB,BZB= " << std::setw(10) << BXB << std::setw(10) << BYB << std::setw(10) << BZB << std::endl;
			                 std::cout << "A SYSTEM BXA,BYA,BZA= " << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
		                         }

		                     BFLD[0] = BXA;  
	 	                     BFLD[1] = BYA; 
		                     BFLD[2] = BZA; 
		                     
std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;
		
		                     return;
	                             } 
//C
//C        -------------------------
//C
	                           else 
	                              {
std::cout << "CASO 6" << std::endl;
//C        UNSPECIFIED FIELD REGION, RETURN WITH ZERO FIELD COMPONENTS
		                       std::cout << "UNKNOWN DIPOLE REGION" << std::endl;
		                       std::cout << "C SYSTEM BXC,BYC,BZC= " << std::setw(10) << BXC << std::setw(10) << BYC << std::setw(10) << BZC << std::endl;
		                       std::cout << "B SYSTEM BXB,BYB,BZB= " << std::setw(10) << BXB << std::setw(10) << BYB << std::setw(10) << BZB << std::endl;
		                       std::cout << "A SYSTEM BXA,BYA,BZA= " << std::setw(10) << BXA << std::setw(10) << BYA << std::setw(10) << BZA << std::endl;
	 	                       std::cout << "!!! Abort current event !!!" << std::endl;
	 	                       G4RunManager::GetRunManager()->AbortEvent();

		                       //jstop = 1;
		                       //ieotri = 1;
		                       BFLD[0] = 0.0; 
		                       BFLD[1] = 0.0; 
		                       BFLD[2] = 0.0; 
		                       
std::cout << "****************************************" << std::endl;
std::cout << "XPOS[0] " << XPOS[0] << std::endl;
std::cout << "XPOS[1] " << XPOS[1] << std::endl;
std::cout << "XPOS[2] " << XPOS[2] << std::endl;
std::cout << "BFLD[0] " << BFLD[0] << std::endl;
std::cout << "BFLD[1] " << BFLD[1] << std::endl;
std::cout << "BFLD[2] " << BFLD[2] << std::endl;
std::cout << "****************************************" << std::endl;
		
		                       return;
	                               }
     }

void DRAGONEMField::mitray_bdip() const
     {
//C****
//C****
//C**** MTYP=1  :    UNIFORM FIELD STANDARD APPROXIMATION
//C**** MTYP=2  :    UNIFORM FIELD MODIFIED ITERATIVE PROCEDURE
//C**** MTYP=3  :    NONUNIFORM FIELD STANDARD APPROXIMATION
//C**** MTYP=4  :    NONUNIFORM FIELD  B=BF/(1+N*DR/R)
//C**** MTYP=5  :    UNIFORM FIELD, CIRCULAR POLE OPTION
//C**** MTYP=6  :    PRETZEL MAGNET
//C****
//C**** THE RELATIONSHIP BETWEEN B0, ......... B12 AND B(I,J) RELATIVE TO
//C**** AXES (Z,X) IS GIVEN BY
//C****
//C****
//C****
//C**** B0  = B( 0, 0 )
//C**** B1  = B( 1, 0 )
//C**** B2  = B( 2, 0 )
//C**** B3  = B( 1, 1 )
//C**** B4  = B( 1,-1 )
//C**** B5  = B( 0, 1 )
//C**** B6  = B( 0, 2 )
//C**** B7  = B( 0,-1 )
//C**** B8  = B( 0,-2 )
//C**** B9  = B(-1, 0 )
//C**** B10 = B(-2, 0 )
//C**** B11 = B(-1, 1 )
//C**** B12 = B(-1,-1 )
//C****
//C****
//C**** MTYP = 1 , 2, 5
//C**** UNIFORM FIELD MAGNETS
//C****
std::cout << "MTYP " << MTYP << std::endl;
     if(MTYP == 1 || MTYP == 2 || MTYP == 5) 
       {
	if(IN == 1 || IN == 3) 
	  {
	   G4double X = TC[0];
	   G4double Y = TC[1];
	   G4double Z = TC[2];
//C****
//C**** MTYP=1,2,5 MAP ROUTINES/INTERPOLATE
//C****
      	   if( IMAP == 0 ) 
      	     {
	      mitray_bdpp( B0, Z, X );
      	      G4double S0 = S;
				
	      if( Y != 0) 
	        {
	         //IMAP = 0, Y != 0
		 if(MTYP != 2) 
		   {
//C****
//C****
//C**** MTYP = 1,5
//C**** NON-MIDPLANE FRINGING FIELD REGION
//C****
	 	    mitray_bdpp(B1, Z + DG, X);
		    mitray_bdpp(B2, Z + 2.0 * DG, X);
		    mitray_bdpp(B3, Z + DG, X + DG);
		    mitray_bdpp(B4, Z + DG, X - DG);
		    mitray_bdpp(B5, Z, X + DG);
		    mitray_bdpp(B6, Z, X + 2.0 * DG);
		    mitray_bdpp(B7, Z, X - DG);
		    mitray_bdpp(B8, Z, X - 2.0 * DG);
		    mitray_bdpp(B9, Z - DG, X);
		    mitray_bdpp(B10, Z - 2.0 * DG, X);
		    mitray_bdpp(B11, Z - DG, X + DG);
		    mitray_bdpp(B12, Z - DG, X - DG);
		    }
		 else 
		    {
//C****
//C**** MTYP = 2
//C**** NON-MIDPLANE FRINGING FIELD REGION
//C****
	 	     mitray_bdpp(B1, 1, 0, true);
		     mitray_bdpp(B2, 2, 0, true);
		     mitray_bdpp(B3, 1, 1, true);
		     mitray_bdpp(B4, 1, -1, true);
		     mitray_bdpp(B5, 0, 1, true);
		     mitray_bdpp(B6, 0, 2, true);
		     mitray_bdpp(B7, 0, -1, true);
		     mitray_bdpp(B8, 0, -2, true);
		     mitray_bdpp(B9, -1, 0, true);
		     mitray_bdpp(B10, -2, 0, true);
		     mitray_bdpp(B11, -1, 1, true);
		     mitray_bdpp(B12, -1, -1, true);
  		     }
//C****
//C**** CALCULATE BX, BY, AND BZ
//C****
		  S = S0;
		  G4double YG1 = Y / DG;
		  G4double YG2 = YG1 * YG1;
		  G4double YG3 = YG2 * YG1;
		  G4double YG4 = YG3 * YG1;
		  BX = YG1 * ((B5 - B7) * 2.0 / 3.0 - (B6 - B8) / 12.0) +
		       YG3 * ((B5 - B7) / 6.0 - (B6 - B8) / 12.0 -
	    	       (B3 + B11 - B4 - B12 - 2.0 * B5 + 2.0 * B7) / 12.0);

		  BY = B0 - YG2 * ((B1 + B9 + B5 + B7 - 4.0 * B0) * 2.0 / 3.0 -
		       (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0) +
		       YG4 * (-(B1 + B9 + B5 + B7 - 4.0 * B0) / 6.0 +
		       (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0 +
		       (B3 + B11 + B4 + B12 - 2.0 * B1 - 2.0 * B9 -
		       2.0 * B5 - 2.0 * B7 + 4.0 * B0) / 12.0);

		  BZ = YG1 * ((B1 - B9) * 2.0 / 3.0 - (B2 - B10) / 12.0) +
		       YG3 * ((B1 - B9) / 6.0 - (B2 - B10) / 12.0 -
		       (B3 + B4 - B11 - B12 - 2.0 * B1 + 2.0 * B9) / 12.0);

		  BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
		  
		  return;
		  } 
	     else 
	        {
//C****
//C**** CONSTANT FIELD REGION
//C****
		 BX = 0.;
		 BY = B0;
		 BZ = 0.;
		 BT = B0;
		
		 return;
		 }
	     }
	else 
	   {
	    //IMAP != 0
 	    mitray_bdmp( B0, Z, X );
	    G4double S0 = 0.;
				
	    if( Y != 0. ) 
	      {
	       mitray_bdmp(B1, Z + DG, X);
	       mitray_bdmp(B2, Z + 2.0 * DG, X);
	       mitray_bdmp(B3, Z + DG, X + DG);
	       mitray_bdmp(B4, Z + DG, X - DG);
	       mitray_bdmp(B5, Z, X + DG);
	       mitray_bdmp(B6, Z, X + 2.0 * DG);
	       mitray_bdmp(B7, Z, X - DG);
	       mitray_bdmp(B8, Z, X - 2.0 * DG);
	       mitray_bdmp(B9, Z - DG, X);
	       mitray_bdmp(B10, Z - 2.0 * DG, X);
	       mitray_bdmp(B11, Z - DG, X + DG);
	       mitray_bdmp(B12, Z - DG, X - DG);
//C****
//C**** CALCULATE BX, BY, AND BZ
//C****
	       S = S0;
	       G4double YG1 = Y / DG;
	       G4double YG2 = YG1 * YG1;
	       G4double YG3 = YG2 * YG1;
	       G4double YG4 = YG3 * YG1;
	       BX = YG1 * ((B5 - B7) * 2.0 / 3.0 - (B6 - B8) / 12.0) +
		    YG3 * ((B5 - B7) / 6.0 - (B6 - B8) / 12.0 -
		    (B3 + B11 - B4 - B12 - 2.0 * B5 + 2.0 * B7) / 12.0);

	       BY = B0 - YG2 * ((B1 + B9 + B5 + B7 - 4.0 * B0) * 2.0 / 3.0 -
		    (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0) +
	  	    YG4 * (-(B1 + B9 + B5 + B7 - 4.0 * B0) / 6.0 +
		    (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0 +
		    (B3 + B11 + B4 + B12 - 2.0 * B1 - 2.0 * B9 -
		    2.0 * B5 - 2.0 * B7 + 4.0 * B0) / 12.0);

	       BZ = YG1 * ((B1 - B9) * 2.0 / 3.0 - (B2 - B10) / 12.0) +
		    YG3 * ((B1 - B9) / 6.0 - (B2 - B10) / 12.0 -
		    (B3 + B4 - B11 - B12 - 2.0 * B1 + 2.0 * B9) / 12.0);

	       BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
		
	       return;
	       }
	    else 
	       {
		//IMAP != 0, Y = 0
		BX = 0.;
		BY = B0;
		BZ = 0.;
		BT = B0;
		
		return;
		}
	   }
         }

       else if(IN == 2) 
              {
		BX = 0.;
      		BY = BF;
      		BZ = 0.;
      		BT = BF;
      		std::cerr << "!!! Abort current run !!!" << std:: endl;
		G4RunManager::GetRunManager()->AbortEvent();
		}
		
             else if(IN == 4) 
		    {
//C****
//C**** CONSTANT FIELD REGION
//C****
		     BX = 0.0;
		     BY = BR;
		     BZ = 0.0;
		     BT = BR;
			
		     return;
		     }
		
		   else 
		      {
		       //IN != 1, 2, 3, 4
		       std::cerr << "0 ERROR - GO TO -  IN BFUN   IN= " << std::setw(5) << IN << std::endl;
		       std::cerr << "!!! Abort current run !!!" << std::endl;
		       BX = 0.;
      		       BY = BF;
      		       BZ = 0.;
      		       BT = BF;
      		       std::cerr << "!!! Abort current run !!!" << std:: endl;
	  	       G4RunManager::GetRunManager()->AbortEvent();
		       }
       }
     else if(MTYP == 3 || MTYP == 4) 
            {
	     mitray_ndip();
	     
	     return;
	     }
          else if(MTYP == 6) 
                 {
	          mitray_bpretz();
	
  	          return;
	          }
		else 
		   {
		    //MTYP != 1, 2, 3, 4, 5, 6
		    std::cerr << "**error** in MITRAY_BDIP" << std::endl;
		    std::cerr << "          Illegal value MTYP=" << MTYP << std::endl;
		    std::cerr << "!!! Abort current run !!!" << std::endl;
		    G4RunManager::GetRunManager()->AbortEvent();
	            }
     }

void DRAGONEMField::mitray_bdpp(G4double& BFLD, G4double Z, G4double X, G4bool BDPPXflag) const
     {
//C****
//C****
//C****
//C**** MTYP=1  :    UNIFORM FIELD STANDARD APPROXIMATION
//C**** MTYP=2  :    UNIFORM FIELD MODIFIED ITERATIVE PROCEDURE
//C****              MORE ACCURATE 3'RD AND HIGHER ORDER CURVATURES
//C**** MTYP=5  :    UNIFORM FIELD, CIRCULAR POLE OPTION
//C****
//C****
	
      if(!BDPPXflag) 
        {
	 //Execute as usual
	 if(MTYP == 1) 
	   {
//C****
//C**** MTYP=1  :    UNIFORM FIELD STANDARD APPROXIMATION
//C****
	    S = ( Z - xmitray_zefb(X) )/D + DELS;
	    }
	 else if(MTYP == 2) 
	        {
//C****
//C**** MTYP=2  :    UNIFORM FIELD, ITERATIVE CALCULATION
//C****
		 mitray_sdip(X, Z);
		 }
 	 else if(MTYP == 5) 
 	        {
//C****
//C**** MTYP=5  :    UNIFORM FIELD, CIRCULAR POLE OPTION
//C****
		 if(std::fabs(RCA >= 1e-8)) 
		   {G4double A = 1./RCA;}
		 else 
		    {S = Z/D + DELS;}
		 }
	       else 
	          {return;}
	 }
       else 
	  {mitray_sdip(Z, X, true);}
	
     G4double CS = C0 + S * (C1 + S * (C2 + S * (C3 + S * (C4 + S * C5))));
     if (std::fabs(CS) > 70.0) 
        {CS = std::copysign(70.0, CS);}

     G4double E = std::exp(CS);
     G4double P0 = 1.0 + E;
     G4double DB = BF - BR;
     BFLD = BR + DB / P0;
	
     return;
     }

void DRAGONEMField::mitray_ndip() const
     {
//C****
//C****
//C**** MTYP = 3 OR 4
//C**** THIS VERSION OF BFUN IS MAINLY FOR NONUNIFORM FIELD MAGNETS
//C**** THE CENTRAL FIELD REGION IS REPRESENTED TO 3'RD ORDER ON-AND-
//C**** OFF THE MIDPLANE BY ANALYTIC EXPRESSIONS. SEE SLAC NO. 75
//C**** FRINGE FIELD REGIONS REPRESENTED BY FERMI TYPE FALL-OFF
//C**** ALONG WITH RADIAL FALL-OFF
//C**** COMPONENTS OF 'B' IN FRINGE REGION EVALUATED BY NUMERICAL METHODS
//C****
//C****
//C**** THE RELATIONSHIP BETWEEN B0, ......... B12 AND B(I,J) RELATIVE TO
//C**** AXES (Z,X) IS GIVEN BY
//C****
//C****
//C**** B0  = B( 0, 0 )
//C**** B1  = B( 1, 0 )
//C**** B2  = B( 2, 0 )
//C**** B3  = B( 1, 1 )
//C**** B4  = B( 1,-1 )
//C**** B5  = B( 0, 1 )
//C**** B6  = B( 0, 2 )
//C**** B7  = B( 0,-1 )
//C**** B8  = B( 0,-2 )
//C**** B9  = B(-1, 0 )
//C**** B10 = B(-2, 0 )
//C**** B11 = B(-1, 1 )
//C**** B12 = B(-1,-1 )
//C****
//C****
      G4double DR1, DR2, DR3, DR4, DR5, DR6, DR7, DR8, DR9, DR10, DR11, DR12;
	
      G4double X = TC[0];
      G4double Y = TC[1];
      G4double Z = TC[2];

      G4double DX = X - XC_OFFSET;
      G4double DZ = Z - ZC_OFFSET;

      G4double RP = std::sqrt(DX * DX + DZ * DZ);
      G4double DR = RP - RB;
std::cout << "Y " << Y << std::endl;
std::cout << "DR " << DR << std::endl;
std::cout << "IN " << IN << std::endl;
std::cout << "IMAP " << IMAP << std::endl;

      if(IN == 1 || IN == 3) 
        {
//C****
//C**** FRINGING FIELD ZONES
//C****
//C**** CHECK IF FIELD MAP CALCULATED
//C****
	 if(IMAP == 0) 
	   {
	    //MTYP=3,4 STANDARD ROUTINES
            G4double ZFB = xmitray_zefb(X);   

	    if (Z > ZFB) 
	       {DR = std::sqrt(DX * DX + std::pow(ZFB - ZC_OFFSET, 2)) - RB;}

std::cout << "ANTES" << std::endl;
std::cout << "X " << X << std::endl;
std::cout << "ZFB " << ZFB << std::endl;
std::cout << "DR " << DR << std::endl;
std::cout << "B0 " << B0 << std::endl;
std::cout << "Z " << Z << std::endl;
                                            
	    mitray_ndpp(B0, Z, X, DR);      
std::cout << "DESPUES" << std::endl;
std::cout << "X " << X << std::endl;
std::cout << "ZFB " << ZFB << std::endl;
std::cout << "DR " << DR << std::endl;
std::cout << "B0 " << B0 << std::endl;
std::cout << "Z " << Z << std::endl;

	    if(Y == 0.) 
	      {
	       BX = 0.;
	       BY = B0;
	       BZ = 0.;
	       BT = B0;
	       
	       return;
	       }
	     else 
	        {
		 //NON MID-PLANE FRINGING FIELD REGION
		 if(Z > ZFB) 
		   {
		    DR1 = DR;
		    DR2 = DR;
		    DR9 = DR;
		    DR10 = DR;

		    G4double XP = X + DG;
		    G4double ZFB = xmitray_zefb(XP);
		    G4double DX = XP - XC_OFFSET;
		    DR3 = std::sqrt(DX * DX + std::pow(ZFB - ZC_OFFSET, 2)) - RB;
		    DR5 = DR3;
		    DR11 = DR3;

		    XP = X - DG;
		    ZFB = xmitray_zefb(XP);
		    DX = XP - XC_OFFSET;
		    DR4 = std::sqrt(DX * DX + std::pow(ZFB - ZC_OFFSET, 2)) - RB;
		    DR7 = DR4;
		    DR12 = DR4;

		    XP = X + 2.0 * DG;
		    ZFB = xmitray_zefb(XP);
		    DX = XP - XC_OFFSET;
		    DR6 = std::sqrt(DX * DX + std::pow(ZFB - ZC_OFFSET, 2)) - RB;

		    XP = X - 2.0 * DG;
		    ZFB = xmitray_zefb(XP);
		    DX = XP - XC_OFFSET;
		    DR8 = std::sqrt(DX * DX + std::pow(ZFB - ZC_OFFSET, 2)) - RB;
		    }
		 else 
		    {
		     DR1  = std::sqrt(DX * DX + std::pow(DZ + DG, 2)) - RB;
		     DR2  = std::sqrt(DX * DX + std::pow(DZ + 2.0 * DG, 2)) - RB;
		     DR3  = std::sqrt(std::pow(DX + DG, 2) + std::pow(DZ + DG, 2)) - RB;
		     DR4  = std::sqrt(std::pow(DX - DG, 2) + std::pow(DZ + DG, 2)) - RB;
		     DR5  = std::sqrt(std::pow(DX + DG, 2) + DZ * DZ) - RB;
		     DR6  = std::sqrt(std::pow(DX + 2.0 * DG, 2) + DZ * DZ) - RB;
		     DR7  = std::sqrt(std::pow(DX - DG, 2) + DZ * DZ) - RB;
		     DR8  = std::sqrt(std::pow(DX - 2.0 * DG, 2) + DZ * DZ) - RB;
		     DR9  = std::sqrt(DX * DX + std::pow(DZ - DG, 2)) - RB;
		     DR10 = std::sqrt(DX * DX + std::pow(DZ - 2.0 * DG, 2)) - RB;
		     DR11 = std::sqrt(std::pow(DX + DG, 2) + std::pow(DZ - DG, 2)) - RB;
		     DR12 = std::sqrt(std::pow(DX - DG, 2) + std::pow(DZ - DG, 2)) - RB;
std::cout << "DR1 " << DR1 << std::endl;
std::cout << "DR2 " << DR2 << std::endl;
std::cout << "DR3 " << DR3 << std::endl;
std::cout << "DR4 " << DR4 << std::endl;
std::cout << "DR5 " << DR5 << std::endl;
std::cout << "DR6 " << DR6 << std::endl;
std::cout << "DR7 " << DR7 << std::endl;
std::cout << "DR8 " << DR8 << std::endl;
std::cout << "DR9 " << DR9 << std::endl;
std::cout << "DR10 " << DR10 << std::endl;
std::cout << "DR11 " << DR11 << std::endl;
std::cout << "DR12 " << DR12 << std::endl;
		     }
				
//C**** CALL NDPP ( B1 , Z + DG, X  , DR1 )
//C**** CALL NDPP ( B2 , Z + 2.*DG, X  , DR2 )
//C**** CALL NDPP ( B3 , Z + DG, X + DG  , DR3 )
//C**** CALL NDPP ( B4 , Z + DG, X - DG  , DR4 )
//C**** CALL NDPP ( B5 , Z , X + DG , DR5 )
//C**** CALL NDPP ( B6 , Z , X + 2.*DG  , DR6 )
//C**** CALL NDPP ( B7 , Z , X - DG , DR7 )
//C**** CALL NDPP ( B8 , Z , X - 2.*DG  , DR8 )
//C**** CALL NDPP ( B9 , Z - DG, X  , DR9 )
//C**** CALL NDPP ( B10, Z - 2.*DG, X, DR10 )
//C**** CALL NDPP ( B11, Z - DG, X + DG  , DR11 )
//C**** CALL NDPP ( B12, Z - DG, X - DG  , DR12 )
//C****
//C****

std::cout << "ANTES" << std::endl;
std::cout << "B0 " << B0 << std::endl;
std::cout << "B1 " << B1 << std::endl;
std::cout << "B2 " << B2 << std::endl;
std::cout << "B3 " << B3 << std::endl;
std::cout << "B4 " << B4 << std::endl;
std::cout << "B5 " << B5 << std::endl;
std::cout << "B6 " << B6 << std::endl;
std::cout << "B7 " << B7 << std::endl;
std::cout << "B8 " << B8 << std::endl;
std::cout << "B9 " << B9 << std::endl;
std::cout << "B10 " << B10 << std::endl;
std::cout << "B11 " << B11 << std::endl;
std::cout << "B12 " << B12 << std::endl;
		 mitray_ndpp(B1, 1, 0, DR1, true);
		 mitray_ndpp(B2, 2, 0, DR2, true);
		 mitray_ndpp(B3, 1, 1, DR3, true);
		 mitray_ndpp(B4, 1, -1, DR4, true);
		 mitray_ndpp(B5, 0, 1, DR5, true);
		 mitray_ndpp(B6, 0, 2, DR6, true);
		 mitray_ndpp(B7, 0, -1, DR7, true);
		 mitray_ndpp(B8, 0, -2, DR8, true);
		 mitray_ndpp(B9, -1, 0, DR9, true); 
std::cout << "-----------------------------------------" << std::endl;
std::cout << "B10 " << B10 << std::endl;
		 mitray_ndpp(B10, -2, 0, DR10, true);		 //ESTE COMO EJEMPLO
std::cout << "B10 " << B10 << std::endl;
std::cout << "-----------------------------------------" << std::endl;
		 mitray_ndpp(B11, -1, 1, DR11, true);
		 mitray_ndpp(B12, -1, -1, DR12, true);
		 
std::cout << "DESPUES" << std::endl;
std::cout << "B0 " << B0 << std::endl;
std::cout << "B1 " << B1 << std::endl;
std::cout << "B2 " << B2 << std::endl;
std::cout << "B3 " << B3 << std::endl;
std::cout << "B4 " << B4 << std::endl;
std::cout << "B5 " << B5 << std::endl;
std::cout << "B6 " << B6 << std::endl;
std::cout << "B7 " << B7 << std::endl;
std::cout << "B8 " << B8 << std::endl;
std::cout << "B9 " << B9 << std::endl;
std::cout << "B10 " << B10 << std::endl;
std::cout << "B11 " << B11 << std::endl;
std::cout << "B12 " << B12 << std::endl;


//C****
//C**** OFF-MIDPLANE FIELD COMPONENTS BX, BY, AND BZ
//C****
		 G4double YG1 = Y / DG;
		 G4double YG2 = YG1 * YG1;
		 G4double YG3 = YG2 * YG1;
		 G4double YG4 = YG3 * YG1;

                 BX = YG1 * ((B5 - B7) * 2.0 / 3.0 - (B6 - B8) / 12.0) +
		      YG3 * ((B5 - B7) / 6.0 - (B6 - B8) / 12.0 -
		      (B3 + B11 - B4 - B12 - 2.0 * B5 + 2.0 * B7) / 12.0);

		 BY = B0 - YG2 * ((B1 + B9 + B5 + B7 - 4.0 * B0) * 2.0 / 3.0 -
		      (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0) +
		      YG4 * (-(B1 + B9 + B5 + B7 - 4.0 * B0) / 6.0 +
		      (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0 +
		      (B3 + B11 + B4 + B12 - 2.0 * B1 - 2.0 * B9 -
		      2.0 * B5 - 2.0 * B7 + 4.0 * B0) / 12.0);

		 BZ = YG1 * ((B1 - B9) * 2.0 / 3.0 - (B2 - B10) / 12.0) +
		      YG3 * ((B1 - B9) / 6.0 - (B2 - B10) / 12.0 -
		      (B3 + B4 - B11 - B12 - 2.0 * B1 + 2.0 * B9) / 12.0);

		 BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);

std::cout << "YG1 " << YG1 << std::endl;
std::cout << "YG2 " << YG2 << std::endl;
std::cout << "YG3 " << YG3 << std::endl;
std::cout << "YG4 " << YG4 << std::endl;
std::cout << "BX " << BX << std::endl;
std::cout << "BY " << BY << std::endl;
std::cout << "BZ " << BZ << std::endl;
std::cout << "BT " << BT << std::endl;

		 return;
		 }
	    }
	  else 
	     {
//C****
//C**** MTYP=3,4 MAP ROUTINES/INTERPOLATE
//C****
//C****
	      mitray_bdmp(B0, Z, X);
	      if(Y == 0) 
	        {
		 BX = 0.;
		 BY = B0;
		 BZ = 0.;
		 BT = B0;
		
		 return;
		 }
	      else 
	         {
		  mitray_bdmp(B1, Z + DG, X);
		  mitray_bdmp(B2, Z + 2.0 * DG, X);
		  mitray_bdmp(B3, Z + DG, X + DG);
		  mitray_bdmp(B4, Z + DG, X - DG);
		  mitray_bdmp(B5, Z, X + DG);
		  mitray_bdmp(B6, Z, X + 2.0 * DG);
		  mitray_bdmp(B7, Z, X - DG);
		  mitray_bdmp(B8, Z, X - 2.0 * DG);
		  mitray_bdmp(B9, Z - DG, X);
		  mitray_bdmp(B10, Z - 2.0 * DG, X);
		  mitray_bdmp(B11, Z - DG, X + DG);
		  mitray_bdmp(B12, Z - DG, X - DG);
				
		  G4double YG1 = Y / DG;
		  G4double YG2 = YG1 * YG1;
		  G4double YG3 = YG2 * YG1;
		  G4double YG4 = YG3 * YG1;

		  BX = YG1 * ((B5 - B7) * 2.0 / 3.0 - (B6 - B8) / 12.0) +
		       YG3 * ((B5 - B7) / 6.0 - (B6 - B8) / 12.0 -
		       (B3 + B11 - B4 - B12 - 2.0 * B5 + 2.0 * B7) / 12.0);

		  BY = B0 - YG2 * ((B1 + B9 + B5 + B7 - 4.0 * B0) * 2.0 / 3.0 -
		       (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0) +
		       YG4 * (-(B1 + B9 + B5 + B7 - 4.0 * B0) / 6.0 +
		       (B2 + B10 + B6 + B8 - 4.0 * B0) / 24.0 +
		       (B3 + B11 + B4 + B12 - 2.0 * B1 - 2.0 * B9 -
		       2.0 * B5 - 2.0 * B7 + 4.0 * B0) / 12.0);

		  BZ = YG1 * ((B1 - B9) * 2.0 / 3.0 - (B2 - B10) / 12.0) +
		       YG3 * ((B1 - B9) / 6.0 - (B2 - B10) / 12.0 -
		       (B3 + B4 - B11 - B12 - 2.0 * B1 + 2.0 * B9) / 12.0);

		  BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
				
		  return;
		  }
	     }
	 }
       else if(IN == 2) 
              {
	       G4double DRR1 = DR/RB;
      	       G4double DRR2 = DRR1*DRR1;
      	       G4double DRR3 = DRR2*DRR1;
      	       G4double DRR4 = DRR3*DRR1;

		if(Y == 0) 
		 {
//C****
//C**** MID-PLANE UNIFORM FIELD REGION
//C****
		  BX = 0.0;
		  BY = 0.0;

		  if(MTYP == 3) 
		    {BY = BF * (1.0 - NDX * DRR1 + BET1 * DRR2 + GAMA * DRR3 + DELT * DRR4);}
		  if(MTYP == 4) 
		    {BY = BF / (1.0 + NDX * DRR1);}

		  BZ = 0.0;
		  BT = BY;

		  return;
		  }
		else 
		   {
//C****
//C**** NON MID-PLANE UNIFORM FIELD REGION
//C****
		    G4double YR1 = Y/RB;
		    G4double YR2 = YR1*YR1;
		    G4double YR3 = YR2*YR1;
		    G4double YR4 = YR3*YR1;
		    G4double RR1 = RB/RP;
		    G4double RR2 = RR1*RR1;
		    G4double RR3 = RR2*RR1;
			
		    if(MTYP == 3) 
		      {
		       G4double BRR = BF * ((-NDX + 2.0 * BET1 * DRR1 + 3.0 * GAMA * DRR2 + 4.0 * DELT * DRR3) * YR1 -
				      (NDX * RR2 + 2.0 * BET1 * RR1 * (1.0 - RR1 * DRR1) +
				      3.0 * GAMA * (2.0 + 2.0 * RR1 * DRR1 - RR2 * DRR2) +
				      4.0 * DELT * (6.0 * DRR1 + 3.0 * RR1 * DRR2 - RR2 * DRR3)) * YR3 / 6.0);

		       BY = BF * (1.0 - NDX * DRR1 + BET1 * DRR2 + GAMA * DRR3 + DELT * DRR4 -
			    0.5 * YR2 * (-NDX * RR1 + 2.0 * BET1 * (1.0 + RR1 * DRR1) +
			    3.0 * GAMA * DRR1 * (2.0 + RR1 * DRR1) +
			    4.0 * DELT * DRR2 * (3.0 + RR1 * DRR1)) +
		  	    YR4 * (-NDX * RR3 + 2.0 * BET1 * (RR3 * DRR1 - RR2) +
			    3.0 * GAMA * (4.0 * RR1 - 2.0 * RR2 * DRR1 + RR3 * DRR2) +
			    4.0 * DELT * (6.0 + 12.0 * RR1 * DRR1 - 3.0 * RR2 * DRR2 + RR3 * DRR3)) / 24.0);
		       BX = BRR * DX / RP;
			    BZ = BRR * DZ / RP;
			    BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
				
		       return;
		       }
		     else if(MTYP == 4) 
		            {
			     G4double DNR1 = 1.0 + NDX * DRR1;
			     G4double DNR2 = DNR1 * DNR1;
			     G4double DNR3 = DNR2 * DNR1;
			     G4double DNR4 = DNR3 * DNR1;
			     G4double DNR5 = DNR4 * DNR1;

			     G4double BRR = BF * NDX * (-YR1 / DNR2 + YR3 * (6.0 * NDX * NDX / DNR4 -
				 	    2.0 * NDX * RR1 / DNR3 - RR2 / DNR2) / 6.0);

			     BY = BF * (1.0 / DNR1 + 0.5 * YR2 * NDX * (-2.0 * NDX / DNR3 + RR1 / DNR2) +
				  YR4 * NDX * (24.0 * NDX * NDX * NDX / DNR5 - 
				  12.0 * NDX * NDX * RR1 / DNR4 - 
				   2.0 * NDX * RR2 / DNR3 - RR3 / DNR2) / 24.0);

			     BX = BRR * DX / RP;
			     BZ = BRR * DZ / RP;
			     BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
				
			     return;
			     } 
			   else 
			      {
			       std::cerr << std::setw(3) << IN << " ERROR -GO TO -  IN BFUN   IN=" 
				 	 << IN << "   MTYP=" << std::setw(4) << MTYP << std::endl;
			       std::exit(EXIT_FAILURE); 
			       }
		   }
	      }
  	    else if(IN == 4) 
  	           {
		    BX = 0;
	  	    BY = BR;
		    BZ = BR;
		    BT = BR;
		
		    return;
	            }
 	         else 
 	            {
		     std::cerr << std::setw(3) << IN << " ERROR -GO TO -  IN BFUN   IN=" 
			       << IN << "   MTYP=" << std::setw(4) << MTYP << std::endl;
		     std::exit(EXIT_FAILURE); 
	             }
}


void DRAGONEMField::mitray_ndpp(G4double& BFLD, G4double Z, G4double X, G4double DR, G4bool NDPPXflag) const
     {
      if(!NDPPXflag) 
        {
         std::cout << "MITRAY_NDPP" << std::endl;
         mitray_sdip(X, Z);    //****************REVISA QUE PASA AQUI con S DENTRO DE mitray_sdip***************************
         }
      
      if(NDPPXflag) 
        {
         std::cout << "MITRAY_NDPPX" << std::endl;
         mitray_sdip(Z, X, true);
         }
 
      G4double DRR1 = DR / RB;
      G4double DRR2 = DRR1 * DRR1;
      G4double DRR3 = DRR2 * DRR1;
      G4double DRR4 = DRR3 * DRR1;
                                             
std::cout << "DRR1 " << DRR1 << std::endl;
std::cout << "DRR2 " << DRR2 << std::endl;
std::cout << "DRR3 " << DRR3 << std::endl;
std::cout << "DRR4 " << DRR4 << std::endl;

      G4double CS = C0 + S * (C1 + S * (C2 + S * (C3 + S * (C4 + S * C5))));
      
      std::cout << "C0 " << C0 << std::endl;
      std::cout << "C1 " << C1 << std::endl;
      std::cout << "C2 " << C2 << std::endl;
      std::cout << "C3 " << C3 << std::endl;
      std::cout << "C4 " << C4 << std::endl;
      std::cout << "C5 " << C5 << std::endl;  
      std::cout << "S " << S << std::endl; 
      std::cout << "CS " << CS << std::endl;
      
      if (std::abs(CS) > 70.0) 
         {CS = std::copysign(70.0, CS);}
      G4double E = std::exp(CS);
      G4double P0 = 1.0 + E;
      G4double DB = BF - BR;
      BFLD = 0.0;

      if(MTYP == 3) 
        {BFLD = BR + (1.0 - NDX * DRR1 + BET1 * DRR2 + GAMA * DRR3 + DELT * DRR4) * DB / P0;}
      else if(MTYP == 4) 
             {BFLD = BR + (1.0 / (1.0 + NDX * DRR1)) * DB / P0;}
             
std::cout << "MTYP " << MTYP << std::endl; 
std::cout << "E " << E << std::endl;
std::cout << "CS " << CS << std::endl;
std::cout << "P0 " << P0 << std::endl;
std::cout << "DB " << DB << std::endl;
std::cout << "BFLD " << BFLD << std::endl;

//C****
//C**** WRITE(6,100) X, Y, Z,  DR, S, BFLD
//C*100 FORMAT( 1P6D15.4 )
//C****
      return;
}

void DRAGONEMField::mitray_bpretz() const
     {
//C****
//C****
//C**** MTYP=6
//C****
//C****
//C**** PRETZEL MAGNET FIELD COMPONENTS
//C**** DG = SMALL NEGATIVE NUMBER
//C****
//C****
      G4double G1 = BF / D;
      G4double Y = TC[1];
      G4double Z = TC[2]; 

      if (Z <= DG) 
         {
	  G4double BY0 = G1 * std::abs(Z);
	  G4double BY1 = BY0 * NDX / Z;
	  G4double BY2 = BY1 * (NDX - 1.0) / Z;
	  G4double BY3 = BY2 * (NDX - 2.0) / Z;
	  G4double BY4 = BY3 * (NDX - 3.0) / Z;

	  BX = 0.0;
	  BY = BY0 - Y * Y * BY2 / 2.0 + Y * Y * Y * Y * BY4 / 24.0;
	  BZ = Y * BY1 - Y * Y * Y * BY3 / 6.0;
	  BT = std::sqrt(BX * BX + BY * BY + BZ * BZ);
		
	  return;
	  } 
      else 
         {
	  BX = 0.0;
	  BY = 0.0;
	  BZ = 0.0;
	
	  return;
	  }
}

//Aqui se asignan los valores a S (EL PROBLEMA ES AQUI)  
// Aqui se esta cambiando S mas veces de lo que deberia (solo una vez por llamada de mitray_ndpp
// para definir los valores del grupo B0-B12)
void DRAGONEMField::mitray_sdip(G4double X, G4double Z, G4bool SIJflag) const
     {
//C****
//C****
//C****
//C**** MTYP=2,3,4  :
//C****
//C****
//C**** FIELD POINT (X,Z)
//C****
//C****
//C**** CHECK TO SEE IF BOUNDARY IS FLAT
//C****
	if(SIJflag) 
	  {
std::cout << "DENTRO DE MITRAY_SIJ" << std::endl;
	   G4double A = (Z * DCS + X * DSN) * DG;
	   G4double DSD = -DCS * (xmitray_zefb(XO+A*DCS) - ZO - A*DSN);
	   S = SS +  ( ( X*DCS - Z*DSN )*DG + DSD )/D;
	   
	   return;
	   }

std::cout << "DENTRO DE MITRAY_SDIP" << std::endl;	
std::cout << "WWWWWWWWWWWWWWWWWWWWWWWW" << std::endl; 
std::cout << "NSRF " << NSRF << std::endl;
	if(NSRF == 0) 
	  {
	   S = Z / D + DELS;
	   SS = S;
	   DCS = 1.0;
	   DSN = 0.0;

	   ZO = xmitray_zefb(X);

std::cout << "D " << D << std::endl;
std::cout << "Z " << Z << std::endl;
std::cout << "DELS " << DELS << std::endl;
std::cout << "S " << S << std::endl;
std::cout << "SS " << SS << std::endl;
std::cout << "ZO " << ZO << std::endl;	
std::cout << "WWWWWWWWWWWWWWWWWWWWWWWW" << std::endl;   
	   return;
  	   }
	else 
	   {
//C****
//C**** FIND POINT ON EFFECTIVE FIELD BOUNDARY THROUGH FIELD POINT
//C**** PARALLEL TO Z-AXIS
//C****
	    G4double XXP, ZZP;
	    G4double ZP = xmitray_zefb(X);

//C****
//C**** INTERVAL OF SEARCH, AZ
//C****
	    G4double AZ = (Z - ZP) / 5.0;
	    G4double ZSIGN = std::copysign(1.0, AZ);
	    G4double AZMAX = std::sqrt(X * X + Z * Z) / 5.0;
	    if (AZ > AZMAX) 
	       {AZ = AZMAX;}

            AZ = std::abs(AZ);
	    G4double XP = X - 5 * AZ;
	    G4double IXP = 1;
	    G4double DP = 1.0e15;

	    for(G4int I = 1; I <= 11; I++) 
	       {
		ZP = xmitray_zefb(XP);
		XXP = X-XP;
		ZZP = Z-ZP;
		G4double DD = XXP*XXP + ZZP*ZZP;
		if(DD < DP) 
		  {
		   IXP = I;
		   DP = DD;
		   }
		XP += AZ;
		}
//C****
//C****  DIVIDE INTERVAL AND REPEAT FOR MORE EXACT
//C****  SHORTEST DISTANCE.
//C****
	    G4double X1 = X + AZ * (IXP - 6); 
	    AZ = AZ / 5.0;           
	    XP = X1 - 5*AZ;
	    IXP = 1;
	    DP = 1.0e15;               	

	    for(G4int I=1; I<=11; I++) 
	       {
		ZP = xmitray_zefb(XP);
		XXP = X-XP;
		G4double ZZP = Z-ZP;
		G4double DD = XXP*XXP + ZZP*ZZP;
		if(DD < DP) 
		  {
		   IXP = I;
		   DP = DD;
	 	   }
		XP += AZ;
		}

	    XO = X1 + AZ * (IXP - 6); 
	    ZO = xmitray_zefb(XO); 
	    G4double XPO = X - XO; 
	    G4double ZPO = Z - ZO;
	    G4double RO = XPO * XPO + ZPO * ZPO;

//C****
//C**** INTERPOLATE FOR MORE ACCURATE LOCATION
//C****
	    if(IXP != 1 && IXP != 11) 
	      {
	       XP  = XO + AZ;
      	       ZP  = xmitray_zefb(XP);
      	       XXP = X-XP;
      	       ZZP = Z-ZP;
      	       G4double R1  = XXP*XXP + ZZP*ZZP;

//C****
//C**** CALCULATE POINT ON THE OTHER SIDE
//C****
 	       G4double XPM = XO - AZ;
	       G4double ZPM = xmitray_zefb(XPM);
	       XXP = X-XPM;
	       ZZP = Z-ZPM;
	       G4double R2  = XXP*XXP + ZZP*ZZP;

	       if(R1 > R2) 
	         {
//C****
//C**** SWAP POINTS
//C****
 		  XP  = XO;
		  ZP  = ZO;
		  R1  = RO;
		  XO  = XPM;
		  ZO  = ZPM;
		  RO  = R2;
		  }

	       G4double X12 = XP-XO;
	       G4double Z12 = ZP-ZO;
	       G4double CC  = X12*X12 + Z12*Z12;
	       XO  = XO + (CC+RO-R1)*AZ/(2*CC);
	       ZO  = xmitray_zefb(XO);
	       G4double XPO = X - XO;
 	       G4double ZPO = Z - ZO;
	       RO = XPO*XPO + ZPO*ZPO;
	       }

	    if(RO < 1.0e-15) 
	      {RO = 1.0e-15;}
	    if(RO > 1.0e15)
	      {RO = 1.0e15;}

	    G4double DZDXO = xmitray_dzdx(XO);
	    G4double COSTH = std::sqrt(1.0 / (1.0 + DZDXO * DZDXO));
	    G4double DELTAX = std::sqrt(RO) * COSTH / 4.0;
//C****
//C****
//C**** WRITE(6,100) X, Z, XO, ZO, COSTH, DELTAX
//C****
//C**** PREPARE TO CALCULATE A PAIR OF EQUALLY SPACED IN X
//C**** DISTANCES ON EITHER SIDE OF RO
//C****
            G4double RINV4 = 1.0/(RO*RO);
//C****
//C**** CALCULATE REPRESENTATIVE DISTANCE
//C****
 	    G4double CX = XO - 2*DELTAX;

	    for(G4int J=1; J<=6; J++) 
	       {
		if(J != 3) 
		  {
		   ZP = xmitray_zefb(CX);
		   G4double XDI = X - CX;
		   G4double ZDI = Z - ZP;
		   G4double RR = XDI*XDI + ZDI*ZDI;
		   
		   if(RR < 1.0e-15) 
		     {RR - 1.0e-15;}
		   if(RR > 1.0e15) 
		     {RR = 1.0e15;}
		   RINV4 = RINV4 + 1.0 / (RR*RR);
		   }
		CX = CX + DELTAX;
		}

	    G4double DP2 = std::sqrt(1.0/RINV4);
	    DP = std::sqrt(DP2);
  	
  	    G4double S = 1.41875 * ZSIGN * DP/D + DELS;
//C****
//C**** Parameters for off midplane calculation
//C****
	    SS = S;
	    G4double DELTA = std::atan(xmitray_dzdx(XO));
	    DCS = std::cos(DELTA);
	    DSN = std::sin(DELTA);
//C****
//C*100 FORMAT( 1P6D15.4 )
//C**** WRITE(6,100) X, Z, DELS, S
//C****
  	    return;
	    }
}

G4double DRAGONEMField::xmitray_zefb(G4double XP) const
         {
          G4double XP2 = XP*XP;
          G4double XP3 = XP2*XP;
          G4double XP4 = XP3 * XP;
          G4double ZEFB= -(S2*XP2 + S3*XP3 + S4*XP4 + S5*XP4*XP + S6*XP4*XP2 + S7*XP4*XP3 
	 	         + S8*XP4*XP4);
    
          return ZEFB;
          }

G4double DRAGONEMField::xmitray_dzdx(G4double XP) const
         {
	  G4double XP2 = XP*XP;
          G4double XP3 = XP2*XP;
          G4double XP4 = XP3 * XP;
          G4double DZDX= -(2.*S2*XP + 3.*S3*XP2+ 4.*S4*XP3 + 5.*S5*XP4 +6.*S6*XP4*XP 
	  	         + 7.*S7*XP4*XP2 + 8.*S8*XP4*XP3);
          return DZDX;
          }

G4double DRAGONEMField::xmitray_dzdx2(G4double XP) const
         {
	  G4double XP2 = XP*XP;
          G4double XP3 = XP2*XP;
          G4double XP4 = XP3 * XP;
	  G4double DZDX2 = -(2.*S2 + 6.*S3*XP+ 12.*S4*XP2 + 20.*S5*XP3 + 30.*S6*XP4 
		           + 42.*S7*XP4*XP + 56.*S8*XP4*XP2);
	
	  return DZDX2;
          }

void DRAGONEMField::mitray_bdmp(G4double& BZZ, G4double Z, G4double X) const
     {
      G4double FZ[3];
      G4double DX = IX + X/DG;
      G4double DZ = IZ + Z/DG;
//C**** NXP = DX
//C**** NZQ = DZ
      G4double NXP = DX + 0.5;
      G4double NZQ = DZ + 0.5;
      G4double PX = DX - NXP;
      G4double QZ = DZ - NZQ;
//C****
//C****
//C**** 6-POINT BIVARIATE INTERPOLATION 'ABRAMOWITZ'
//C****
//C****
//C****      BZZ =  BF*( ( QZ*(QZ-1.) * BZMAP( NXP, NZQ-1, IR, IMAP ) +
//C****     1        PX*(PX-1.) * BZMAP( NXP-1, NZQ, IR, IMAP ) +
//C****     2        PX*(PX-2.*QZ+1.) * BZMAP( NXP+1, NZQ, IR, IMAP ) +
//C****     3        QZ*(QZ-2.*PX+1.)*BZMAP( NXP, NZQ+1, IR, IMAP )  )/2.
//C****     4        PX*QZ * BZMAP( NXP+1, NZQ+1, IR, IMAP ) +
//C****     5        (1.+PX*QZ-PX*PX-QZ*QZ) * BZMAP( NXP, NZQ, IR, IMAP )
//C****
//C****
      G4double QZ2 = QZ * QZ;
      G4double QZ3 = QZ2 * QZ;
      G4double QZ4 = QZ3 * QZ;

      for(G4int I=1; I<=3; I++) 
         {
	  G4int NXX = NXP - 2 + I;
  	  G4int NZQ_int = static_cast<G4int>(NZQ);
	  G4double BM2 = BZMAP[NXX - 1][NZQ_int - 3][IR - 1][IMAP - 1];
	  G4double BM1 = BZMAP[NXX - 1][NZQ_int - 2][IR - 1][IMAP - 1];
 	  G4double B00 = BZMAP[NXX - 1][NZQ_int - 1][IR - 1][IMAP - 1];
 	  G4double BP1 = BZMAP[NXX - 1][NZQ_int][IR - 1][IMAP - 1];
 	  G4double BP2 = BZMAP[NXX - 1][NZQ_int + 1][IR - 1][IMAP - 1];


	  G4double A1 = ((BP1 - BM1) * 8 - BP2 + BM2) / 12.0;
	  G4double A2 = ((BP1 + BM1) * 16 - BP2 - BM2 - 30 * B00) / 24.0;
	  G4double A3 = ((BM1 - BP1) * 2 + BP2 - BM2) / 12.0;
	  G4double A4 = (-4 * (BP1 + BM1) + BP2 + BM2 + 6 * B00) / 24.0;

	  FZ[I-1] = B00 + A1 * QZ + A2 * QZ2 + A3 * QZ3 + A4 * QZ4;
	  }

      G4double C1 = (FZ[2] - FZ[0]) / 2.0;
      G4double C2 = (FZ[2] + FZ[0] - 2.0 * FZ[1]) / 2.0;
      BZZ = BF * (FZ[1] + C1 * PX + C2 * PX * PX);
	
      return;
}

void DRAGONEMField::mitray_fmap() 
     {
//C****
//C****
//C**** Control Section for calculating Dipole Field Maps on a
//C**** Rectangular grid.
//C****
//C****
//C****
      G4int PRNT = 0;
//C****
//C****
//C**** BZMAP(NX,NZ,IR,IMAP)
//C****    NX         : X-POSITION INDEX
//C****    NZ         : Z-POSITION INDEX
//C****    IR=1       : ENTRANCE FRINGE FIELD
//C****    IR=2       : EXIT FRINGE FIELD
//C****    IMAP=(1-5) : IDENTIFIES FIELD MAP
//C****
//C****
//C****
//C**** CLEAR FIELD MAP AND INDEX ARRAYS
//C****
     for(G4int I=0; I<5; I++) 
        {JMAP[I] = 0;}
	
     for(G4int I1=0; I1<101; I1++) 
        {
	 for(G4int I2=0; I2<101; I2++) 
	    {
	     for(G4int I3=0; I3<2; I3++) 
	        {
		 for(G4int I4=0; I4<5; I4++) 
		    {BZMAP[I1][I2][I3][I4] = 0.;}
		 }
	     }
	 }
//C****
//C**** CYCLE THROUGH ELEMENTS TO FIND DIPOLES WHICH NEED FIELD MAPS
//C**** CALCULATED.
//C****
     for(G4int I=0; I<no; I++) 
        {
	 if(idata[I-1] != 2) 
	    continue;
//C****
//C**** CHECK FOR MAP INDEX AND WHETHER DIPOLE NEEDS A MAP TO BE
//C**** CALCULATED
//C****
	 MTYP = data[4][I];
	 IMAP = data[5][I];
	 
	 if(IMAP == 0) 
	   continue;
//C****
//C**** CHECK FOR VALID MTYP'S
//C****
	 if(MTYP <= 5) 
	   {
//C****
//C**** CHECK TO SEE IF MAP WITH THIS INDEX IMAP=(1-5) HAS ALLREADY
//C**** BEEN CALCULATED
//C**** CHECK IMAP INDEX LIMITS
//C****
	    if(IMAP <= 4) 
	      {
	       if(JMAP[IMAP-1] != 0) 
	         {continue;}
//C****
//C**** IDENTIFY DIPOLE MAGNETIC ELEMENT USED TO CALCULATE FIELD MAP
//C**** FOR INDEX IMAP.
//C****
               else 
                  {
		   JMAP[IMAP-1] = I;
//C****
//C**** CALCULATE FIELD MAP FOR INDEX IMAP
//C****
		   if(PRNT == 0) 
		     {
		      std::cout << "1" << std::endl; 
		      std::cout << " FIELD MAP PARAMETERS " << std::endl;
		      std::cout << std::setw(5) << " ";
		      std::cout << "IMAP   NO   MTYP   IR   NXLO  NXHI  NZLO  NZHI  " << std::endl;
		      PRNT = 1;
		      mitray_dmap(I);
		      continue;
		      }
		   }
	       }
             else 
                {
		 std::cerr << std::endl;
		 std::cerr << " *** FATAL ERROR **** - ELEMENT NO= " << std::setw(5) << I << std::endl;
		 std::cerr << "    EXCEEDS MAXIMUM FIELD MAP INDEX" << std::endl;
		 std::cerr << "  IMAP= " << std::setw(5) << IMAP << std::endl;
		 }
	    }
          else 
             {
	      std::cout << std::endl
			<< " ***WARNING*** FIELD MAPS NOT IMPLEMENTED FOR THIS"
			<< "\n MTYP:  NO=" << std::setw(5) << I 
			<< "  MTYP=" << std::setw(5) << MTYP 
			<< "  IMAP=" << std::setw(5) << IMAP 
			<< std::endl;
//C****
//C**** RESET INVALID FIELD-MAP
//C****
	      data[5][I] = 0.;
	      continue;
	      }
	}
     return;
}

void DRAGONEMField::mitray_dmap(G4int II) 
     {
//C****
//C****
//C**** CALCULATE FIELD MAPS
//C****
//C****
      G4double LF1, LF2;

      G4double NXLMT[5][2][2], NZLMT[5][2][2];

      G4int NXMAX = 101;
      G4int NZMAX = 101;
      G4int NZ1MAX = 60;
      G4int NZ2MAX = 40;

      IX = NXMAX/2 + 1;
      IZ = NZ2MAX + 1; 

      LF1  = data[0][II - 1];  
      LF2  = data[2][II - 1];  
      DG   = data[3][II - 1]; 
      MTYP      = data[4][II - 1]; 
      IMAP      = data[5][II - 1];  
      D    = data[12][II - 1]; 
      RB   = data[13][II - 1]; 
      BF   = data[14][II - 1]; 
      G4double ALPHA= data[16][II - 1]; 
      G4double BETA = data[17][II - 1]; 
      NDX  = data[18][II - 1]; 
      BET1 = data[19][II - 1]; 
      GAMA = data[20][II - 1]; 
      DELT = data[21][II - 1]; 
      G4double Z11  = data[24][II - 1]; 
      G4double Z12  = data[25][II - 1]; 
      G4double Z21  = data[26][II - 1]; 
      G4double Z22  = data[27][II - 1];
      G4double BR1  = data[40][II - 1]; 
      G4double BR2  = data[41][II - 1]; 
      G4double WDE  = data[48][II - 1]; 
      G4double WDX  = data[49][II - 1]; 

      if (MTYP == 0) 
         {MTYP = 1;}

      if (WDE == 0) 
         {
	  WDE = 5.0 * D;
 	  data[48][II - 1] = WDE;
	  }
      if (WDX == 0) 
         {
	  WDX = 5.0 * D;
	  data[49][II - 1] = WDX;
	  }

      G4double DX1 = (WDE + 2.0 * std::abs(Z11) * std::tan(std::abs(ALPHA / 57.29578))) / (NXMAX - 7);
      G4double DX2 = (WDX + 2.0 * std::abs(Z22) * std::tan(std::abs(BETA / 57.29578))) / (NXMAX - 7);
      G4double DZ11 = (LF1 + std::abs(Z11)) / (NZ1MAX - 3);
      G4double DZ12 = (LF1 + std::abs(Z12)) / (NZ2MAX - 3);
      G4double DZ21 = (LF2 + std::abs(Z21)) / (NZ2MAX - 3);
      G4double DZ22 = (LF2 + std::abs(Z22)) / (NZ1MAX - 3);
      G4double DGI = DG;
      
      if(DX1 > DG) DG = DX1;
      if(DX2 > DG) DG = DX2;
      if(DZ11 > DG) DG = DZ11;
      if(DZ12 > DG) DG = DZ12;
      if(DZ21 > DG) DG = DZ21;
      if(DZ22 > DG) DG = DZ22;
      if(DG != DGI) 
        {
	 data[3][II-1] = DG;
	 std::cerr << "***WARNING*** INPUT DG CHANGED TO STAY WITHIN ARRAY LIMITS : "
        	   << "DG(Input)=" << std::fixed << std::setprecision(3) << std::setw(10) << DGI
          	   << "   DG(Calc.)=" << std::setw(10) << DG << std::endl;
	 }
//C****
//C**** IR=1
//C****
//C****
      G4int IFLAG = 0;
      IR = 1;
      G4double NDX1 = (WDE + 2. * std::abs(Z11) * std::tan(std::abs(ALPHA / 57.29578))) / (2. * DG);
      G4double NXLO = IX - NDX1 - 3;
      G4double NXHI = IX + NDX1 + 3;
      G4double NZLO = IZ - 3 + (Z12 - LF1) / DG;
      G4double NZHI = IZ + 3 + (Z11 + LF1) / DG;

      NXLMT[IMAP-1][IR-1][0] = NXLO;
      NXLMT[IMAP-1][IR-1][1] = NXHI;
      NZLMT[IMAP-1][IR-1][0] = NZLO;
      NZLMT[IMAP-1][IR-1][1] = NZHI;
//C****
//C**** CHECK IF INDEX .LT. 1 ; PRINT WARNING
//C**** CHECK IF NX .GT. NXMAX  ; PRINT WARNING
//C**** CHECK IF NZ .GT. NZMAX  ; PRINT WARNING
//C****
      if (NXLO <= 0) 
         {
	  NXLMT[IMAP-1][IR-1][0] = 1;
  	  IFLAG = 1;
	  }

      if (NXHI > NXMAX) 
         {
	  NXLMT[IMAP-1][IR-1][1] = NXMAX;
	  IFLAG = 1;
	  }

      if (NZLO <= 0) 
         {
	  NZLMT[IMAP-1][IR-1][0] = 1;
	  IFLAG = 1;
	  }

      if (NZHI > NZMAX) 
         {
	  NZLMT[IMAP-1][IR-1][1] = NZMAX;
	  IFLAG = 1;
	  }

      if(IFLAG != 0) 
        {
	 std::cout << "***WARNING*** MAP INDICES EXCEED LIMITS -RESET-\n"
                   << " NO=" << std::setw(4) << II
                   << " IMAP=" << std::setw(4) << IMAP
                   << " IR=" << std::setw(4) << IR << "\n"
                   << " NXLO=" << std::setw(4) << NXLO
                   << " NXHI=" << std::setw(4) << NXHI
                   << " NZLO=" << std::setw(4) << NZLO
                   << " NZHI=" << std::setw(4) << NZHI << std::endl;
	 }
//C****
//C**** IR=2
//C****
      IFLAG = 0;
      IR = 2;
      NDX1 = (WDX + 2.0 * std::abs(Z22) * std::tan(std::abs(BETA / 57.29578))) / (2.0 * DG);
      NXLO = IX - NDX1 - 3;
      NXHI = IX + NDX1 + 3;
      NZLO = IZ - 3 + (Z21 - LF2) / DG;
      NZHI = IZ + 3 + (Z22 + LF2) / DG;
      NXLMT[IMAP-1][IR-1][0] = NXLO;
      NXLMT[IMAP-1][IR-1][1] = NXHI;
      NZLMT[IMAP-1][IR-1][0] = NZLO;
      NZLMT[IMAP-1][IR-1][1] = NZHI;
//C****
//C**** CHECK IF INDEX .LT. 1 ; PRINT WARNING
//C**** CHECK IF NX .GT. NXMAX  ; PRINT WARNING
//C**** CHECK IF NZ .GT. NZMAX  ; PRINT WARNING
//C****
      if (NXLO <= 0) 
         {
	  NXLMT[IMAP-1][IR-1][0] = 1;
	  IFLAG = 1;
	  }
	
      if (NXHI > NXMAX) 
         {
	  NXLMT[IMAP-1][IR-1][1] = NXMAX;
 	  IFLAG = 1;
	  }
	
      if (NZLO <= 0) 
         {
	  NZLMT[IMAP-1][IR-1][0] = 1;
	  IFLAG = 1;
	  }
	
      if (NZHI > NZMAX) 
         {
	  NZLMT[IMAP-1][IR-1][1] = NZMAX;
	  IFLAG = 1;
	  }
	
      if(IFLAG != 0) 
        {
  	 std::cout << "***WARNING*** MAP INDICES EXCEED LIMITS -RESET-\n"
                   << " NO=" << std::setw(4) << II
                   << " IMAP=" << std::setw(4) << IMAP
                   << " IR=" << std::setw(4) << IR << "\n"
                   << " NXLO=" << std::setw(4) << NXLO
                   << " NXHI=" << std::setw(4) << NXHI
                   << " NZLO=" << std::setw(4) << NZLO
                   << " NZHI=" << std::setw(4) << NZHI << std::endl;
	 }
//C****
//C**** CALCULATE MAPS FOR ENTRANCE AND EXIT FRINGE FIELDS
//C****
      for(IR=1; IR<=2; IR++) 
         {
	  if(IR == 1) 
	    {
             //C****
             //C**** SETUP ENTRANCE FRINGE FIELD PARAMETERS
             //C****
	     XC = RB * std::cos(ALPHA / 57.29578);
 	     ZC = -RB * std::sin(ALPHA / 57.29578);
	     BR = BR1;
			
	     C0 = data[28][II - 1];
	     C1 = data[29][II - 1];
	     C2 = data[30][II - 1];
	     C3 = data[31][II - 1];
	     C4 = data[32][II - 1];
	     C5 = data[33][II - 1];
	     DELS = data[44][II - 1];
	     RCA = data[46][II - 1];
	     S2 = data[50][II - 1] / RB + RCA / 2.0;
	     S3 = data[51][II - 1] / std::pow(RB, 2);
	     S4 = data[52][II - 1] / std::pow(RB, 3) + std::pow(RCA, 3) / 8.0;
	     S5 = data[53][II - 1] / std::pow(RB, 4);
	     S6 = data[54][II - 1] / std::pow(RB, 5) + std::pow(RCA, 5) / 16.0;
	     S7 = data[55][II - 1] / std::pow(RB, 6);
	     S8 = data[56][II - 1] / std::pow(RB, 7) + std::pow(RCA, 7) / 25.6;
//C****
//C**** CHECK IF WE HAVE A FLAT BOUNDARY
//C****       NSRF=0 FLAT
//C****           =1 CURVED
//C****
	     NSRF = 1;
	     if ((S2 == 0.0) && (S3 == 0.0) && (S4 == 0.0) &&
		(S5 == 0.0) && (S6 == 0.0) && (S7 == 0.0) &&
		(S8 == 0.0)) 
		{NSRF = 0;}
	     }
	   if(IR == 2) 
	     {
//C****
//C**** SETUP EXIT FRINGE FIELD PARAMETERS
//C****
//C****
	      XC = -RB * std::cos(BETA / 57.29578);
	      ZC = -RB * std::sin(BETA / 57.29578);
	      BR = BR2;
	      C0 = data[34][II - 1];
	      C1 = data[35][II - 1];
	      C2 = data[36][II - 1];
	      C3 = data[37][II - 1];
	      C4 = data[38][II - 1];
	      C5 = data[39][II - 1];
	      DELS = data[45][II - 1];
	      RCA = data[47][II - 1];
	      S2 = data[57][II - 1] / RB + RCA / 2.0;
	      S3 = data[58][II - 1] / pow(RB, 2);
	      S4 = data[59][II - 1] / pow(RB, 3) + pow(RCA, 3) / 8.0;
	      S5 = data[60][II - 1] / pow(RB, 4);
	      S6 = data[61][II - 1] / pow(RB, 5) + pow(RCA, 5) / 16.0;
	      S7 = data[62][II - 1] / pow(RB, 6);
	      S8 = data[63][II - 1] / pow(RB, 7) + pow(RCA, 7) / 25.6;
//C****
//C**** CHECK IF WE HAVE A FLAT BOUNDARY
//C****       NSRF=0 FLAT
//C****           =1 CURVED
//C****
	      NSRF = 1;
	      if ((S2 == 0.0) && (S3 == 0.0) && (S4 == 0.0) &&
		 (S5 == 0.0) && (S6 == 0.0) && (S7 == 0.0) &&
		 (S8 == 0.0)) 
		 {NSRF = 0;}
	      }
		
	    NXLO = NXLMT[IMAP - 1][IR - 1][0];
	    NXHI = NXLMT[IMAP - 1][IR - 1][1];
	    NZLO = NZLMT[IMAP - 1][IR - 1][0];
	    NZHI = NZLMT[IMAP - 1][IR - 1][1];
//C****
//C**** MTYP = 3, 4
//C****
	    if(MTYP == 3 || MTYP || 4) 
	      {
	       for(G4int I=NXLO; I<=NXHI; I++) 
	          {
		   for(G4int J=NZLO; J<=NZHI; J++) 
		      {
		       G4double X = (I-IX) * DG;
		       G4double Z = (J-IZ) * DG;
		       G4double DX = X - XC;
		       G4double DZ = Z - ZC;
		       G4double ZFB = xmitray_zefb(X);
			
		       if(Z > ZFB) 
		         DZ = ZFB - ZC;
		       G4double DR = std::sqrt(DX*DX + DZ*DZ) - RB;
		       mitray_ndpp(B0, Z, X, DR);
		       BZMAP[I-1][J-1][IR-1][IMAP-1] = B0/BF;
		       }
		   }
	       }
//C****
//C**** MTYP = 0, 1, 2, 5
//C****
            if(MTYP = 0 || MTYP == 1 || MTYP == 2 || MTYP == 5) 
              {
	       for(G4int I=NXLO; I<=NXHI; I++) 
	          {
		   for(G4int J=NZLO; J<=NZHI; J++) 
		      {
		       G4double X = (I-IX) * DG;
		       G4double Z = (J-IZ) * DG;
		       mitray_bdpp(B0, Z, X);
		       BZMAP[I-1][J-1][IR-1][IMAP-1] = B0/BF;
		       }
		   }
		}

	    std::cout << std::setw(3) << "" 
		      << std::setw(6) << IMAP
		      << std::setw(6) << II
		      << std::setw(6) << MTYP
		      << std::setw(6) << IR
		      << std::setw(6) << NXLO
		      << std::setw(6) << NXHI
		      << std::setw(6) << NZLO
		      << std::setw(6) << NZHI
		      << std::endl;
	 }
//C****
//C**** PRINT MAPS
//C***
      G4int NP = 1000;
      if(NP <= 100) 
        {
	 for(IR=1; IR<=2; IR++) 
	    {
	     if(IR == 1) 
	       {
		std::cout << '1' << std::endl
          		  << "   ENTRANCE FRINGING FIELD MAP :  IMAP=" << std::setw(3) << IMAP << std::endl;
		}
			
	     if(IR == 2) 
	       {
		std::cout << '1' << std::endl
          		  << "   EXIT FRINGING FIELD MAP :  IMAP=" << std::setw(3) << IMAP << std::endl;
	        }
	     NXLO = NXLMT[IMAP-1][IR-1][0];
 	     NXHI = NXLMT[IMAP-1][IR-1][1];
	     NZLO = NZLMT[IMAP-1][IR-1][0];
	     NZHI = NZLMT[IMAP-1][IR-1][1];
			
	     for(G4int I=NXLO; I<=NXHI; I++) 
	        {
		 G4int J1 = I;
		 G4int J2 = I+14;
		 
		 if(J1 > NXLO) 
		   {std::cout << "1";}
				
		 if(J2 > NXHI) 
		   J2 = NXHI;

		 std::cout << "    NX";
		 for (G4int J = J1; J <= J2; J++) 
		     {std::cout << std::setw(8) << J;}
		 std::cout << std::endl;

		 std::cout << "     X";
		 for (G4int J = J1; J <= J2; J++) 
		     {std::cout << std::setw(8) << std::fixed << std::setprecision(3) << (J - 51) * DG;} 
		 std::cout << std::endl;
 		 std::cout << "  NZ     Z" << std::endl;
		 
		 G4int L_int = static_cast<G4int>(L);
		 for(G4int KSK = NZLO; KSK <= NZHI; KSK++) 
		    {
		     std::cout << std::setw(4) << L 
			       << std::fixed << std::setprecision(3) << (L - IZ) * DG; 
		     for (G4int J = J1; J <= J2; ++J) 
		         {std::cout << std::fixed << std::setprecision(4) << BZMAP[J-1][L_int - 1][IR - 1][IMAP - 1];}
		     std::cout << std::endl;			
		     }
		  }
	   }
       } 
     return;
}


}
