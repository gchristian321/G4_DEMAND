#include "DRAGONEMField.hh"     //!geant


namespace DRAGON {

void DRAGONEMField::mitray_edipol(G4double* DATA, G4double* XPOS, G4double* EFLD) const {
    /*
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    C                                                                      C
    C     Subroutine for electric dipole, in GEANT implementation of       C
    C     MIT-RAYTRACE                                                     C
    C     adapted from:                                                    C
    C     Subroutine EDIPL (NO, NP, T, TP, NUM)  by S. Kowalski            C
    C                                                                      C
    C      TC(1) to  TC(6) =  (  X,  Y,  Z, VX, VY, VZ )                   C
    C     DTC(1) to DTC(6) =  ( VX, VY, VZ, VXDOT, VYDOT, VZDOT )          C
    C     T = TIME                                                         C
    C                                                                      C
    C     Input:  DATA(i)    array containing parameters of the            C
    C                        electric dipole                               C
    C             XPOS(i),   i=1,3  contain the X,Y,Z coordinates of the   C
    C                        field point, in the A-axis coordinate system  C
    C                        of the electric dipole.                       C
    C                                                                      C
    C     Output: EFLD(i),   i=1,3  contain the electric field components  C
    C                        Ex, Ey, Ez at the specified field point, in   C
    C                        the A-axis coordinate system of the E-dipole. C
    C                                                                      C
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  */

    G4double LF1, LF2, LU1, EFF;
    G4String REGION;
    G4double EXA, EYA, EZA;
    G4double Z11, Z12, Z21, Z22;
    
//C     S.Yen addition  July 31, 1999
//C     Extract coordinates of the field point in the
//C     A-axis coordinate system of the electric dipole

    G4double XA = XPOS[0];
    G4double YA = XPOS[1];
    G4double ZA = XPOS[2];
    REGION = " ";
    
    if (ldiag) {
       std::cout << std::endl;
        std::cout << std::setw(50) << std::setfill('-') << "" << std::endl;
        std::cout << " ENTER SUBROUTINE EDIPOLE" << std::endl;
        std::cout << " XA,YA,ZA=";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << std::setw(10) << XA << " ";
        std::cout << std::setw(10) << YA << " ";
        std::cout << std::setw(10) << ZA << std::endl;
        }

//C
//C     Extract the parameters for the electric dipole
//C

    LF1  = DATA[0];
    LU1  = DATA[1];
    LF2  = DATA[2];
    DG   = DATA[3];
    A    = DATA[10];
    B    = DATA[11];
    D    = DATA[12];
    RB   = DATA[13];
    EFF  = DATA[14];
    PHI  = DATA[15];
    EC2  = DATA[16];
    EC4  = DATA[17];
    WE   = DATA[18];
    WC   = DATA[19];
    Z11  = DATA[24];
    Z12  = DATA[25];
    Z21  = DATA[26];
    Z22  = DATA[27];
    
//C
//C     SY Addition Sept 20/99
//C     We must define the bounds of the entrance fringe field
//C     in the XB axis system.  The entrance fringe field region is
//C     defined by ( XBMIN < XB < XBMAX )  &  ( Z12 < ZB < Z11 )
//C     If the field point is outside the gap width between the 
//C     electrodes, then we assume that the field is zero.
//C

    G4double XBMAX = (0.5*D);
    G4double XBMIN = -XBMAX;
    
//C
//C     Similarly the bounds of the exit fringe field are defined in
//C     the XC axis system.  The exit fringe field is defined by
//C     ( XCMIN < XC < XCMAX ) & ( Z21 < ZC < Z22 )
//C

    G4double XCMAX = XBMAX;
    G4double XCMIN = XBMIN;
    
//C
//C     SY Addition Sept 20/99
//C     First zero the E-field components in case we need to abort
//C

    for(G4int i=0; i<3; i++) {
        EFLD[i] = 0.;
    }

    if(WE == 0.) WE = 1000. * RB;
    
    BX = 0.;
    BY = 0.;
    BZ = 0.;
    EX = 0.;
    EY = 0.;
    EZ = 0.;
    ET = 0.;
    S = 0.;
 
//C
//C     TRANSFORM FROM INITIAL ENTRANCE (A-SYSTEM) COORDINATES TO ENTRANCE
//C     EFB (B-SYSTEM) COORDINATES.
//C

    G4double XB =  -XA;
    G4double YB = YA;
    G4double ZB = (A-ZA);
 
//C
//C     TRANSFORM FROM THE B-SYSTEM COORDINATES TO THE EXIT EFB (C-SYSTEM)
//C     COORDINATES.
//C
//C     (These are the same equations used for B-->C tranformation for
//C      the magnetic dipole, but with ALPHA=BETA=0)
//C
    G4double COPAB =std::cos((PHI)/deg_rad);
    G4double SIPAB =std::sin((PHI)/deg_rad);
    G4double COSPB =std::cos((PHI/2.)/deg_rad);
    G4double SINPB =std::sin((PHI/2.)/deg_rad);
    G4double SIP2 =std::sin((PHI/2.)/deg_rad);
    G4double ZC = -ZB  *COPAB +  XB  *SIPAB -2.*RB*SIP2*COSPB;
    G4double XC = -ZB  *SIPAB -  XB  *COPAB -2.*RB*SIP2*SINPB;
    G4double YC = YB;
//C
//C     Print out coordinates if in diagostic mode      
//C
    if (ldiag) {
        std::cout << " XB,YB,ZB=" << std::fixed << std::setprecision(3) 
              << XB << " " << YB << " " << ZB << std::endl;
    
        std::cout << " XC,YC,ZC=" << std::fixed << std::setprecision(3) 
              << XC << " " << YC << " " << ZC << std::endl;
    
        std::cout << " Z11=" << std::fixed << std::setprecision(3) << Z11 
              << "  Z12=" << Z12 
              << "  Z21=" << Z21 
              << "  Z22=" << Z22 << std::endl;
    
        std::cout << " XBMIN=" << std::fixed << std::setprecision(3) << XBMIN 
              << "  XBMAX=" << XBMAX 
              << "  XCMIN=" << XCMIN 
              << "  XCMAX=" << XCMAX << std::endl;
    
        std::cout << std::endl;
    }
    
//C
//C     Now determine which region we are in.  Choices are:
//C     Before start of entrance fringe field ("entrance far field") IN=-99
//C     entrance fringe field region  (IN=1)
//C     "uniform" field region (IN=2)
//C     exit fringe field region (IN=3)
//C     after the end of the exit fringe field ("exit far field") IN=+99
//C     Because of the possibility of confusing entrance and exit regions,
//C     we first check for entrance/exit fringe field, then uniform field,
//C     then entrance/exit far field regions, to make sure that we get the
//C     most important regions first.
//C
//C**** IN DESIGNATES MAGNET REGIONS FOR BFUN
//C

std::cout << "YYYYYYYYYYYYYY" << std::endl;
std::cout << "ZB " << ZB << "\n";
std::cout << "ZC " << ZC << "\n";
std::cout << "Z11 " << Z11 << "\n";
std::cout << "Z12 " << Z12 << "\n";
std::cout << "Z21 " << Z21 << "\n";
std::cout << "Z22 " << Z22 << "\n";
std::cout << "XB " << XB << "\n";
std::cout << "XC " << XC << "\n";
std::cout << "XBMIN " << XBMIN << "\n";
std::cout << "XCMIN " << XCMIN << "\n";
std::cout << "XBMAX " << XBMAX << "\n";
std::cout << "YYYYYYYYYYYYYY" << std::endl;

    if (ZB <= Z11 && ZB > Z12 && XB >= XBMIN && XB <= XBMAX) {
std::cout << "SIIIIIIIII1" << "\n";
//C
//C        *************************
//C        *                       *
//C        * ENTRANCE FRINGE FIELD *
//C        *                       *
//C        *************************
//C
//C        Entrance fringe field region, B-axis coordinates are used.
//C
        IN = 1;
        XC_OFFSET = RB;
        ZC_OFFSET = 0.0;
        EF = EFF;

        //Get Enge coefficients Cn for entrance fringe field.

        C0 = DATA[34];
        C1 = DATA[35];
        C2 = DATA[36];
        C3 = DATA[37];
        C4 = DATA[38];
        C5 = DATA[39];
//C
//C       Load the B-axis coordinates into TC(1), TC(2), TC(3), because this 
//C       is where the subroutine EDIP expects to find the coordinates
//C
        TC[0] = XB;
        TC[1] = YB;
        TC[2] = ZB;
//C
//C       Call subroutine EDIP to calculate the E-field components
//C
        mitray_edip();
//C
//C       The E-field components EX, EY, EZ have been passed back
//C       via common block /MITRAY11/, and are in the B-axis system.
//C       Transform the E-field components back to the A-axis system.
//C       Actually there is no change, since the A-axis and B-axis systems
//C       are related by a simple translation for the electrostatic dipole,
//C       with no rotation.
//C
        EXA = EX;
        EYA = EY;
        EZA = EZ;
//C
//C       Print out diagnostics if required
//C
        if (ldiag) {
            std::cout << "ENTRANCE FRINGE FIELD REGION" << std::endl;
            std::cout << "B SYSTEM " << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
    
            std::cout << "ENTRANCE FRINGE FIELD REGION" << std::endl;
            std::cout << "B SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
        }
//C
//C       Load EXA, EYA, EZA into array EFLD and exit.
//C
        EFLD[0] = EXA;
        EFLD[1] = EYA;
        EFLD[2] = EZA;
        
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
//            std::cout << "EFLD[0] " << EFLD[0] << "\n";
//            std::cout << "EFLD[1] " << EFLD[1] << "\n";
//            std::cout << "EFLD[2] " << EFLD[2] << "\n";
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
               
        return;
        }
//C
//    -------------------------------------
//C
 else if(ZC > Z21 & ZC <= Z22 && XC >= XCMIN && XC <= XCMAX) {
std::cout << "SIIIIIIIII2" << "\n";
//C
//C        *********************
//C        *                   *
//C        * EXIT FRINGE FIELD *
//C        *                   *
//C        *********************
//C
//C        SET INDICATOR IN=3 for exit fringe field
//C
        IN = 3;
        XC_OFFSET=-RB;     //! ADDED NOV 23/99
        ZC_OFFSET=0.;      //! ADDED NOV 23/99
        EF=-EFF;           //! ADDED NOV 23/99
//C
//C       Get Enge coefficients Cn for the exit fringe field      
//C
        C0 = DATA[34];
        C1 = DATA[35];
        C2 = DATA[36];
        C3 = DATA[37];
        C4 = DATA[38];
        C5 = DATA[39];
//C     
//C        LOAD THE C AXIS COORDINATES INTO ARRAY TC, BECAUSE THIS IS
//C        WHERE SUBROUTINE BDIP EXPECTS TO FIND THEM.
//C
        TC[0] = XC;
        TC[1] = YC;
        TC[2] = ZC;
        
//C        Call subroutine EDIP to calculate the E-field components
//C        in the C-axis system
//C
        mitray_edip();
//C
//C        The electric field components Ex, Ey, Ez are in the C-axis
//C        system.  Transform them to the A-axis system.
//C
       mitray_edip_ECtoEA(EX,EY,EZ,EXA,EYA,EZA);
//C
//C       Print out diagnostics if required
//C
        if (ldiag) {
            std::cout << "EXIT FRINGE FIELD REGION" << std::endl;
            std::cout << "C SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
            std::cout << "EXIT FRINGE FIELD REGION" << std::endl;
            std::cout << "C SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
        }
//C
//C       Load EXA, EYA, EZA into array EFLD and exit.
//C
        EFLD[0] = EXA;
        EFLD[1] = EYA;
        EFLD[2] = EZA;
              
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
//            std::cout << "EFLD[0] " << EFLD[0] << "\n";
//            std::cout << "EFLD[1] " << EFLD[1] << "\n";
//            std::cout << "EFLD[2] " << EFLD[2] << "\n";
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
        
        return;
        }
//C
//C       --------------------------------------
//C
 else if(ZB <= Z12 && ZC <= Z21) {
std::cout << "SIIIIIIIII3" << "\n";
//C
//C        ************************
//C        *                      *
//C        * UNIFORM FIELD REGION *
//C        *                      *
//C        ************************
//C      
//C        Set indicator IN=2 for uniform field region
//C
        IN = 2;
        XC_OFFSET = -RB;
        ZC_OFFSET = 0.0;
        EF = -EFF;
        S = 0.;
              
//C        LOAD THE C-AXIS COORDINATES INTO ARRAY TC BECAUSE THIS
//C        IS WHERE SUBROUTINE EDIP EXPECTS TO FIND THEM
//C
        TC[0] = XC;
        TC[1] = YC;
        TC[2] = ZC;
//C
//C        CALL SUBROUTINE EDIP TO CALCULATE THE E-FIELD COMPONENTS IN
//C        THE C-AXIS SYSTEM
//C
	mitray_edip();
	
//C        EX, EY, EZ ARE THE E-FIELD COMPONENTS IN THE C-SYSTEM.
//C        TRANSFORM THEM TO THE A-SYSTEM
//C

     mitray_edip_ECtoEA(EX, EY, EZ, EXA, EYA, EZA);
//C
//C        Print out diagnostics if required
//C
        if (ldiag) {
            std::cout << "UNIFORM FIELD REGION" << std::endl;
            std::cout << "C SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
                << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
                << EXA << " " << EYA << " " << EZA << std::endl;
    
            std::cout << "UNIFORM FIELD REGION" << std::endl;
            std::cout << "C SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
                << EX << " " << EY << " " << EZ << std::endl;
            std::cout << "A SYSTEM  EX,EY,EZ=" << std::fixed << std::setprecision(3) 
                << EXA << " " << EYA << " " << EZA << std::endl;
        }
//C
//C        Load EXA, EYA, EZA into array EFLD and exit.
//C
        EFLD[0] = EXA;
        EFLD[1] = EYA;
        EFLD[2] = EZA;
        
//                    std::cout << "LLLLLLLLLLLLLLL" << "\n";
 //           std::cout << "EFLD[0] " << EFLD[0] << "\n";
 //           std::cout << "EFLD[1] " << EFLD[1] << "\n";
 //           std::cout << "EFLD[2] " << EFLD[2] << "\n";
 //           std::cout << "LLLLLLLLLLLLLLL" << "\n";
        
        return;
      
        }
//C
//C     -----------------------------------------
//C
    else if(ZB > Z11 && XB >= XBMIN && XB <= XBMAX) {
std::cout << "SIIIIIIIII4" << "\n";
//C
//C        **********************
//C        *                    *
//C        * ENTRANCE FAR FIELD *
//C        *                    *
//C        **********************
//C        
//C        SET INDICATOR IN=-99
//C 
        IN = -99;
        EXA = 0.;
        EYA = 0.;
        EZA = 0.;
//C
//C       PRINT DIAGNOSTICS IF REQUIRED
//C
        if (ldiag) {
            std::cout << "DIPOLE ENTRANCE FAR FIELD REGION" << std::endl;
            std::cout << "A SYSTEM  EXA,EYA,EZA=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
            std::cout << "DIPOLE ENTRANCE FAR FIELD REGION" << std::endl;
            std::cout << "A SYSTEM  EXA,EYA,EZA=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
        }
//C
//C        LOAD E-FIELD COMPONENTS INTO ARRAY EFLD AND EXIT
//C
        EFLD[0] = EXA;
        EFLD[1] = EYA;
        EFLD[2] = EZA;
        
//                    std::cout << "LLLLLLLLLLLLLLL" << "\n";
//            std::cout << "EFLD[0] " << EFLD[0] << "\n";
//            std::cout << "EFLD[1] " << EFLD[1] << "\n";
//            std::cout << "EFLD[2] " << EFLD[2] << "\n";
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
        
        return;
   }
//C
//C        -----------------------------------------
//C

 else if(ZC > Z22 && XC > XCMIN && XC <= XCMAX) {
std::cout << "SIIIIIIIII5" << "\n";
//C
//C           ******************
//C           *                *
//C           * EXIT FAR FIELD *
//C           *                *
//C           ******************
//C
//C        SET INDICATOR IN=+99
//C
        IN = +99;
        EXA = 0.;
        EYA = 0.;
        EZA = 0.;
//C
//C        PRINT DIAGNOSTICS IF REQUIRED
//C
        if (ldiag) {
            std::cout << "DIPOLE EXIT FAR FIELD REGION" << std::endl;
            std::cout << "A SYSTEM  EXA,EYA,EZA=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;

            std::cout << "DIPOLE EXIT FAR FIELD REGION" << std::endl;
            std::cout << "A SYSTEM  EXA,EYA,EZA=" << std::fixed << std::setprecision(3) 
              << EXA << " " << EYA << " " << EZA << std::endl;
        }
//C
//C        LOAD E-FIELD COMPONENTS INTO ARRAY EFLD AND EXIT
//C
        EFLD[0] = EXA;
        EFLD[1] = EYA;
        EFLD[2] = EZA;
        
//                    std::cout << "LLLLLLLLLLLLLLL" << "\n";
 //           std::cout << "EFLD[0] " << EFLD[0] << "\n";
 //           std::cout << "EFLD[1] " << EFLD[1] << "\n";
 //           std::cout << "EFLD[2] " << EFLD[2] << "\n";
 //           std::cout << "LLLLLLLLLLLLLLL" << "\n";

        return;
    }

//C
//C      ----------------------------------------------------
//C
 else {
//C        UNSPECIFIED FIELD REGION, RETURN WITH ZERO FIELD COMPONENTS
//C
  if(ldiag) {
            std::cout << "UNKNOWN ELECTRIC DIPOLE REGION" << std::endl;

            std::cout << "A SYSTEM XA,YA,ZA=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XA << " " 
                << std::setw(12) << YA << " " 
                << std::setw(12) << ZA << std::endl;

            std::cout << "B SYSTEM XB,YB,ZB=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XB << " " 
                << std::setw(12) << YB << " " 
                << std::setw(12) << ZB << std::endl;

            std::cout << "C SYSTEM XC,YC,ZC=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XC << " " 
                << std::setw(12) << YC << " " 
                << std::setw(12) << ZC << std::endl;

            std::cout << "RETURN E-FIELD EXA=0, EYA=0, EZA=0" << std::endl;

            std::cout << "UNKNOWN ELECTRIC DIPOLE REGION" << std::endl;

            std::cout << "A SYSTEM XA,YA,ZA=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XA << " " 
                << std::setw(12) << YA << " " 
                << std::setw(12) << ZA << std::endl;

            std::cout << "B SYSTEM XB,YB,ZB=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XB << " " 
                << std::setw(12) << YB << " " 
                << std::setw(12) << ZB << std::endl;

            std::cout << "C SYSTEM XC,YC,ZC=" << std::fixed << std::setprecision(4)
                << std::setw(12) << XC << " " 
                << std::setw(12) << YC << " " 
                << std::setw(12) << ZC << std::endl;

            std::cout << "RETURN E-FIELD EXA=0, EYA=0, EZA=0" << std::endl;

            //SET ALL E-FIELD COMPONENTS TO ZERO
            
            EFLD[0] = 0;
            EFLD[1] = 0;
            EFLD[2] = 0;
            
//            std::cout << "LLLLLLLLLLLLLLL" << "\n";
 //           std::cout << "EFLD[0] " << EFLD[0] << "\n";
 //           std::cout << "EFLD[1] " << EFLD[1] << "\n";
 //           std::cout << "EFLD[2] " << EFLD[2] << "\n";
 //           std::cout << "LLLLLLLLLLLLLLL" << "\n";
             
            return;
        }

        return;
    }

//=========================================================================

}

void DRAGONEMField::mitray_edip() const {


//C****
//C**** CALCULATES E-FIELD COMPONENTS FOR A CYLINDRICAL
//C**** ELECTROSTATIC DEFLECTOR
//C****

    G4double K;

//C****
//C****
//C       ADDED BY SY OCT 30/99  INITIAL VALUES OF E-FIELD        
//C
    EX = 0.;
    EY = 0.;
    EZ = 0.;

    G4double X = TC[0];
    G4double Y = TC[1];
    G4double Z = TC[2];
    G4double DX = X - XC_OFFSET;
    G4double RP2 = DX * DX + Z * Z;

std::cout << "TC[0] " << TC[0] << "\n";
std::cout << "TC[1] " << TC[1] << "\n";
std::cout << "TC[2] " << TC[2] << "\n";
std::cout << "DX " << DX << "\n";
std::cout << "RP2 " << RP2 << "\n";

    //TRACED UP TO HERE: NOW NEED TO PRINT VARIABLES TO TEST THE FOLLOWING LINES
    if(IN == 1 || IN ==3) {
//C****
//C****   FRINGE FIELD REGION    IN=1 OR IN=3
//C****
        G4double RP = std::sqrt(RP2);
        
        if( std::abs(X) < std::abs(XC_OFFSET) ) {
            G4double DR = RP-RB;
            G4double SINT = Z/RP;
            G4double COST = std::abs(XC_OFFSET-X)/RP;
            G4double THETA = std::asin(SINT);
            G4double S = THETA * RB / D + EC2 * Y * Y / (WE * WE) + EC4 * std::pow((Y / WE), 4);
            
            std::cout << "DR " << DR << "\n";
            std::cout << "SINT " << SINT << "\n";
            std::cout << "COST " << COST << "\n";
            std::cout << "THETA " << THETA << "\n";
            std::cout << "S " << S << "\n";
                        
            std::cout << "D " << D << "\n";
            std::cout << "S " << S << "\n";
            std::cout << "RE " << RE << "\n";
            std::cout << "G1 " << G1 << "\n";
            std::cout << "G2 " << G1 << "\n";
            std::cout << "G3 " << G3 << "\n";
            std::cout << "G4 " << G4 << "\n";
            std::cout << "G5 " << G5 << "\n";
            std::cout << "G6 " << G6 << "\n";
                        
            mitray_edpp( D, S, RE, G1, G2, G3, G4, G5, G6 );
            
            std::cout << "QQQQQQQQQQQQQQQQQQQQ" << "\n";
            std::cout << "D " << D << "\n";
            std::cout << "S " << S << "\n";
            std::cout << "RE " << RE << "\n";
            std::cout << "G1 " << G1 << "\n";
            std::cout << "G2 " << G1 << "\n";
            std::cout << "G3 " << G3 << "\n";
            std::cout << "G4 " << G4 << "\n";
            std::cout << "G5 " << G5 << "\n";
            std::cout << "G6 " << G6 << "\n";
            std::cout << "QQQQQQQQQQQQQQQQQQQQ" << "\n";
            
            G4double DRR  = DR/RB;
            G4double DRR2 = DRR*DRR;
            G4double DRR3 = DRR2*DRR;
            G4double DRR4 = DRR3*DRR;
            G4double EFR = EF * RB * ( RE - DRR2 * G2 / 2.0
                            + DRR3 * G2 / 2.0
                            + DRR4 * (G4 - 11.0 * G2) / 24.0 ) / RP;
            G4double EFT = EF * RB * ( DRR * G1 - DRR2 * G1 / 2.0
                            + DRR3 * (2.0 * G1 - G3) / 6.0
                            - DRR4 * (G1 - G3) / 4.0 ) / RP;
            EX = EFR*COST - EFT*SINT;
            EZ = EFR*SINT + EFT*COST;
            ET = std::sqrt( EX * EX + EY * EY + EZ * EZ);
            
            std::cout << "DRR " << DRR << "\n";
            std::cout << "DRR2 " << DRR2 << "\n";
            std::cout << "DRR3 " << DRR3 << "\n";
            std::cout << "DRR4 " << DRR4 << "\n";
            std::cout << "EFR " << EFR << "\n";
            std::cout << "EFT " << EFT << "\n";
            std::cout << "EX " << EX << "\n";
            std::cout << "EZ " << EZ << "\n";
            std::cout << "ET " << ET << "\n";
            
            if(IN == 1) {
                 EZ = -EZ;
            }

            return;
        } 
        
        else {
            std::cerr << " ** ERROR ** IN SUBROUTINE MITRAY_EDIP" << std::endl;
            std::cerr << "               X = " << X << ",    XC_OFFSET = " 
                << XC_OFFSET << std::endl;
            return;
        }
    }
    
    else if(IN == 2) {
//C****
//C****   UNIFORM FIELD REGION   IN=2
//C****
        EX = EF * RB * DX / RP2;
        EY = 0.;
        EZ = EF * RB * Z / RP2;
        ET = std::sqrt(EX * EX + EZ * EZ);
      
        return;
    } 
    
    else {
//C       ADDED BY SY  OCT 30/99 RETURN WITH ZERO FIELD IN CASE OF ILLEGAL
//C       VALUE OF 'IN'
        std::cout << " ** ERROR ** -GO TO-  IN SUBROUTINE MITRAY_EDIP" << std::endl;
        std::cout << "               INVALID VALUE  IN = " << IN << std::endl;
        return;
    }   

}


void DRAGONEMField::mitray_edpp(G4double D, G4double S, G4double& RE, G4double& G1, G4double& G2, G4double& G3, G4double& G4, G4double& G5, G4double& G6) const {

//C****
//C**** CALCULATE S; DETERMINE E-FIELD IN FRINGE REGIONS
//C****
    G4double K;

    G4double S2 = S * S;
    G4double S3 = S2 * S;
    G4double S4 = S2 * S2;
    G4double S5 = S4 * S;
    G4double CS = C0 + C1 * S + C2 * S2 + C3 * S3 + C4 * S4 + C5 * S5;
    G4double RBD = RB / D;
    G4double CP1 = (C1 + 2.0 * C2 * S + 3.0 * C3 * S2 + 4.0 * C4 * S3 + 5.0 * C5 * S4) * RBD;
    G4double CP2 = (2.0 * C2 + 6.0 * C3 * S + 12.0 * C4 * S2 + 20.0 * C5 * S3) * RBD * RBD;
    G4double CP3 = (6.0 * C3 + 24.0 * C4 * S + 60.0 * C5 * S2) * std::pow(RBD, 3);
    G4double CP4 = (24.0 * C4 + 120.0 * C5 * S)*pow(RBD,4.);

    std::cout << "-------------------" << "\n";
    std::cout << "S2 " << S2 << "\n";
    std::cout << "S3 " << S3 << "\n";
    std::cout << "S4 " << S4 << "\n";
    std::cout << "S5 " << S5 << "\n";
    std::cout << "CS " << CS << "\n";
    std::cout << "RBD " << RBD << "\n";
    std::cout << "CP1 " << CP1 << "\n";
    std::cout << "CP2 " << CP2 << "\n";
    std::cout << "CP3 " << CP3 << "\n";
    std::cout << "CP4 " << CP4 << "\n";
    std::cout << "-------------------" << "\n";
  
    if(std::abs(CS) > 70.) CS = copysign(70, CS);
    
    G4double E = std::exp(CS);
    RE = 1./(1.+E);
    G4double ERE = E*RE;
    G4double ERE1 = ERE*RE;
    G4double ERE2 = ERE*ERE1;
    G4double ERE3 = ERE*ERE2;
    G4double ERE4 = ERE*ERE3;

    G4double CP12 = CP1*CP1;
    G4double CP13 = CP1*CP12;
    G4double CP14 = CP12*CP12;
    G4double CP22 = CP2*CP2;

    std::cout << "KKKKKKKKKKKKKK" << "\n";
    std::cout << "CS " << CS << "\n";
    std::cout << "E " << E << "\n";
    std::cout << "RE " << RE << "\n";
    std::cout << "ERE " << ERE << "\n";
    std::cout << "ERE1 " << ERE1 << "\n";
    std::cout << "ERE2 " << ERE2 << "\n";
    std::cout << "ERE3 " << ERE3 << "\n";
    std::cout << "ERE4 " << ERE4 << "\n";
    std::cout << "CP12 " << CP12 << "\n";
    std::cout << "CP13 " << CP13 << "\n";
    std::cout << "CP14 " << CP14 << "\n";
    std::cout << "CP22 " << CP22 << "\n";   
    
    G1 = -CP1*ERE1;
    G2 = -( CP2+CP12 )*ERE1 + 2.*CP12 * ERE2;
    G3 = -(CP3 + 3.*CP1*CP2 + CP13) * ERE1 + 6.*(CP1*CP2 + CP13)*ERE2 - 6.*CP13*ERE3;
    G4 = -(CP4 + 4.*CP1*CP3 + 3.*CP22 + 6.*CP12*CP2 + CP14)*ERE1 +
        (8.*CP1*CP3 + 36.*CP12*CP2 + 6.*CP22 + 14.*CP14)*ERE2 -
        36.*(CP12*CP2 + CP14)*ERE3 + 24.*CP14*ERE4;
        
    std::cout << "G1 " << G1 << "\n";
    std::cout << "G2 " << G2 << "\n";
    std::cout << "G3 " << G3 << "\n";
    std::cout << "G4 " << G4 << "\n";
    std::cout << "KKKKKKKKKKKKKK" << "\n";

    return;

}



void DRAGONEMField::mitray_edip_ECtoEA(G4double EXC, G4double EYC, G4double EZC, G4double& EXA, G4double& EYA, G4double& EZA) const {
//C
//C     TRANSFORM E-FIELD COMPONENTS FROM SYSTEM C TO SYSTEM A COORDINATES
//C     OF AN ELECTOSTATIC DIPOLE.
//C
//C     INPUT:  EXC, EYC, EZC     E-FIELD COMPONENTS IN C-AXIS SYSTEM
//C     OUTPUT: EXA, EYA, EZA     E-FIELD COMPONENTS IN A-AXIS SYSTEM
//C

//C
//C     We use here the same formulae as for the magnetic dipole, but
//C     with ALPHA=BETA=0
//C
    G4double COPAB = std::cos(-PHI/deg_rad);
    G4double SIPAB = std::sin(-PHI/deg_rad);

//C
//C     Next 3 lines not needed
//C
//C     COSPB = DCOS(-PHI/2. /57.29577951D0)
//C     SINPB = DSIN(-PHI/2. /57.29577951D0)
//C     SIP2 = DSIN(-PHI/2. /57.29577951D0)
//C
//C    NOW ROTATE TO GET E-FIELD COMPONENTS IN B-SYSTEM, IN TERMS OF C-SYSTEM
//C

    G4double EZB = -EZC * COPAB + EXC* SIPAB;
    G4double EXB = -EZC * SIPAB - EXC * COPAB;
    G4double EYB = EYC;

    //A-axis and B-axis systems are related by simple translation and
    //no rotation for the electrostatic dipole, so the E-field components
    //are the same in both B-axis and A- axis systems.

    EXA=EXB;
    EYA=EYB;
    EZA=EZB;
    
//    std::cout << "LLLLLLLLLLLLLLL" << "\n";
//    std::cout << "EXA " << EXA << "\n";
//    std::cout << "EYA " << EYA << "\n";
//    std::cout << "EZA " << EZA << "\n";
//    std::cout << "LLLLLLLLLLLLLLL" << "\n";
    

    return;

}





}
