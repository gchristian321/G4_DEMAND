#ifndef geom_mpole_h
#define geom_mpole_h 1

#include "G4Types.hh"

//C. *** Geometry COMMON BLOCK for all POLES magnets in the beam-line

      static const G4int max_mpole = 20;

      extern G4int nmpole;                          //! # of mpoles

      extern G4double pos_mpole[max_mpole][3];      //! position  of mpole in beam-line

      extern G4double dx_mpole[3][max_mpole];       //! displacement of A-frame

      extern G4double efblength_mpole[max_mpole];   //! mpole gap width
      extern G4double r_mpole[max_mpole];           //! mpole aperture radius

      extern G4double z11_mpole[max_mpole];         //! RAYTRACE - z11
      extern G4double z22_mpole[max_mpole];         //! RAYTRACE - z22

      extern G4int jcol_mpole[2][max_mpole];        //! collimator type
      extern G4double xcol_mpole[2][max_mpole];     //! collimator x-position
      extern G4double ycol_mpole[2][max_mpole];     //! collimator y-position
      extern G4double dxcol_mpole[2][max_mpole];    //! collimator x-size
      extern G4double dycol_mpole[2][max_mpole];    //! collimator y_size


#endif
