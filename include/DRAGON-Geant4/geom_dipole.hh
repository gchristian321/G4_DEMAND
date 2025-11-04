#ifndef geom_dipole_h
#define geom_dipole_h 1

#include "G4Types.hh"

//C. *** Geometry COMMON BLOCK for all DIPOLE magnets in the beam-line

//    The radius of curvature is POSITIVE 
//        if bend from entrance to exit is in a clockwise direction
//    The radius of curvature is NEGATIVE
//        if bend from entrance to exit is in a counter clockwise direction

      static const G4int max_dipole = 9;

      extern G4int ndipole;                        //! # of dipoles

      extern G4double pos_dipole[3][max_dipole];    //! position  of dipole in beam-line

      extern G4int irot_dipole[max_dipole];        //! rotation matrix index of A-frame
      extern G4double dx_dipole[3][max_dipole];     //! displacement of A-frame 

      extern G4double gap_dipole[max_dipole];      //! dipole gap width
      extern G4double phi_dipole[max_dipole];      //! dipole bend angle
      extern G4double r_dipole[max_dipole];        //! dipole radius of curvature
      extern G4double dr_dipole[max_dipole];       //! dipole delta_r of radius

      extern G4double alpha_dipole[max_dipole];    //! RAYTRACE - alpha
      extern G4double beta_dipole[max_dipole];    //! RAYTRACE - beta

      extern G4double z11_dipole[max_dipole];      //! RAYTRACE - z11
      extern G4double z22_dipole[max_dipole];      //! RAYTRACE - z22

      extern G4int jcol_dipole[2][max_dipole];      //! collimator type
      extern G4double xcol_dipole[2][max_dipole];   //! collimator x-position
      extern G4double ycol_dipole[2][max_dipole];   //! collimator y-position
      extern G4double dxcol_dipole[2][max_dipole];  //! collimator x-size
      extern G4double dycol_dipole[2][max_dipole];  //! collimator y_size


#endif
