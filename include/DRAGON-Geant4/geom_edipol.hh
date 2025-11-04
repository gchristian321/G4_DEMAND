#ifndef geom_edipol_h
#define geom_edipol_h 1

#include "G4Types.hh"


//C. *** Geometry COMMON BLOCK for all electrostatic DEFLECTORS in the beam-line

//C.    The radius of curvature is POSITIVE 
//C.        if bend from entrance to exit is in a clockwise direction
//C.    The radius of curvature is NEGATIVE
//C.        if bend from entrance to exit is in a counter clockwise direction

      static const G4int max_edipol = 9;
 
      extern G4int nedipol;                         //! # of electrostatic deflectors

      extern G4double pos_edipol[3][max_edipol];    //! position  of deflector in beam-line

      extern G4int irot_edipol[max_edipol];         //! rotation matrix index of A-frame
      extern G4double dx_edipol[3][max_edipol];     //! displacement of A-frame 

      extern G4double gap_edipol[max_edipol];       //! deflector gap width
      extern G4double phi_edipol[max_edipol];       //! deflector bend angle
      extern G4double r_edipol[max_edipol];         //! deflector radius of curvature
      extern G4double dr_edipol[max_edipol];        //! deflector delta_r of radius

      extern G4double z11_edipol[max_edipol];       //! RAYTRACE - z11
      extern G4double z22_edipol[max_edipol];       //! RAYTRACE - z22

      extern G4int jcol_edipol[2][max_edipol];      //! collimator type
      extern G4double xcol_edipol[2][max_edipol];   //! collimator x-position
      extern G4double ycol_edipol[2][max_edipol];   //! collimator y-position
      extern G4double dxcol_edipol[2][max_edipol];  //! collimator x-size
      extern G4double dycol_edipol[2][max_edipol];  //! collimator y_size




#endif
