#ifndef mitray_setup_h
#define mitray_setup_h 1

#include "G4Types.hh"
#include "G4String.hh"



//C *** MITRAY COMMON BLOCK

      static const G4int nmax = 300, mmax = 75;

//C   MMAX is maximum number of parameters for element
//C   NMAX is maximum number of elements in beam transport system

      static G4int no;
      static G4int idata[nmax];

      static G4String ititle[nmax]; 

       static G4double data[mmax][nmax];


//C   NO        I*4               Number of elements in beam transport 
//C                                 system 

//C   IDATA(i)  I*4(nmax)         Type code for ith element  (2=DIPO, 
//C                                  3=EINZ,  7=EDIP,  8=VELS,  9=POLE, 
//C                                 10=MULT, 11=SHRT, 12=DRIF, 13=COLL,
//C                                 14=SOLE, 15=LENS, 16=ACCE)
//C
//C   ITITLE(i) CHAR*12(nmax)     User-specified name for ith element
//C                                 (e.g. B1, QUAD2, ARGUS, GEORGE)
//C
//C   DATA(j,i) REAL*8(mmax,nmax) Parameters for the ith element
//C                                 (field strengths, radii, etc.)
//C                                 Up to 75 parameters, for each of
//C                                 200 elements.
//C                                     

#endif
