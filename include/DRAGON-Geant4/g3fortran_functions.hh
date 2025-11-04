#ifndef g3fortran_functions_h
#define g3fortran_functions_h 1

#include "G4Types.hh"

void vzero(G4double* array, G4int dimension);
void ucopy(G4double* array1, G4double* array2, G4int dimension);
void gtrmul(G4double DX1[3],G4double RMAT1[10],G4double DX2[3],G4int IROT,G4double DXNEW[3],G4double RMATN[10]);
void grmul(G4double RMAT[10],G4int IROT,G4double RMATN[10]);


#endif


