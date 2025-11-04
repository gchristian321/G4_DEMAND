#ifndef geant3functions_h
#define geant3functions_h 1

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

namespace DRAGON {

//Geant3 helper functions transferred into Geant4
//Made by Ben Marlow, DRAGON Co-op Student on 10 December, 2024

//===========================================================

//VECTOR FUNCTIONS

//===========================================================
static const G4int max_photon = 10;

void ucopy(G4double* array1, G4double* array2, G4int dimension);

void ucopy(G4String* names1, G4String* names2, G4int dimension);

void ucopy(G4int* array1, G4int* array2, G4int dimension);

void ucopy(G4double* array1, G4double array2[][max_photon], G4int rows, G4int col_index);

void vzero(G4double* array, G4int dimension);

void vzero(G4int* array, G4int dimension);

void vfill(G4double* array, G4int dimension, G4double value);

void vfill(G4int* array, G4int dimension, G4int value);

void vfill(G4int* array, G4int dimension, G4int value);

void vscale(G4double* output_array, G4double scale, G4double* input_array, G4int dimension);

void vunit(G4double* input_array, G4int dimension);

void vadd(G4double* array1, G4double* array2, G4double* array3, G4int dimension);

G4double vdotn(const G4double vec1[3][10], const G4double vec2[3][10], G4int col);

G4double Vmod(G4double* vector, G4int n);

void ublank(G4String* array, G4int dimension);

void sortzv(G4double* a, G4int* index, G4int n, G4int mode, G4int nway = 0, G4int nsort = 0);

//===========================================================

//MATRIX FUNCTIONS

//===========================================================

void gitran(G4double x[3], G4double dx[3], G4int irot, G4double xnew[3], G4double* theta = nullptr, G4double* phi = nullptr);

void gdtom(G4double xd[3], G4double xm[3], G4int iflag, G4String devname);

void gmtod(G4double xm[3], G4double xd[3], G4int iflag, G4String devname);

}

#endif


