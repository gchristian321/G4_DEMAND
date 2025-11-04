#ifndef geom_sole_h
#define geom_sole_h 1

#include "G4Types.hh"


//C. *** Geometry COMMON BLOCK for SOLENOID magnets in the beam-line

      static const G4int max_sole = 9;
   
      extern G4int nsole;                          //! # of solenoids

      extern G4double pos_sole[3][max_sole];       //! position of solenoid in beam-line

      extern G4double dx_sole[3][max_sole];        //! displacement of A-frame

      extern G4double efblength_sole[max_sole];    //! solenoid length
      extern G4double r_sole[max_sole];            //! solenoid radius

      extern G4double z11_sole[max_sole];          //! RAYTRACE - z11
      extern G4double z22_sole[max_sole];          //! RAYTRACE - z22

      extern G4int jcol_sole[2][max_sole];         //! collimator type
      extern G4double xcol_sole[2][max_sole];      //! collimator x-position
      extern G4double ycol_sole[2][max_sole];      //! collimator y-position
      extern G4double dxcol_sole[2][max_sole];     //! collimator x-size
      extern G4double dycol_sole[2][max_sole];     //! collimator y_size



#endif
