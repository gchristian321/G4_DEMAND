#ifndef NuclearMasses_h
#define NuclearMasses_h  1

#include <cmath>
#include "globals.hh"
//

namespace g4gen {

/// Utility class to lookup nuclear masses from the 2012 AME.
class NuclearMasses 
{
private:
  
  // Default constructor - this class should exist only once!
  NuclearMasses();

public:

  // Destructor (generated)
  ~NuclearMasses() { };

  // Following values migrate to AME12
  enum  {nEntries = 3353,MaxA = 295, ZMax = 120}; 

public:

  // Operation: GetMassExcess
  //   values imported from The Ame2003 atomic mass evaluation (II)  
  static G4double GetMassExcess(G4int Z, G4int A); 

  // Operation: GetAtomicMass .. in Geant4 Energy units!
  //      Atomic_Mass =  MassExcess + A*amu_c2
  static G4double GetAtomicMass(G4int Z, G4int A);

  // Operation: GetNuclearMass
  //      Nuclear_Mass = Atomic_Mass - electronMass
  static G4double GetNuclearMass(G4int Z, G4int A);

  // Operation: GetBindingEnergy
  static G4double GetBindingEnergy(G4int Z, G4int A);

  // Operation: GetBetaDecayEnergy
  static G4double GetBetaDecayEnergy(G4int Z, G4int A);

  // Is the nucleus (A,Z) in table?
  static G4bool IsInTable(G4int Z, G4int A);

  static G4int MaxZ(G4int A);
  static G4int MinZ(G4int A);

	// GET FROM STRING //

	static void GetZAFromSymbol(const G4String& symbol, G4int* Z, G4int * A);
	
  // Operation: GetMassExcess
  //   values imported from The Ame2003 atomic mass evaluation (II)  
  static G4double GetMassExcess(const G4String& symbol); 

  // Operation: GetAtomicMass .. in Geant4 Energy units!
  //      Atomic_Mass =  MassExcess + A*amu_c2
  static G4double GetAtomicMass(const G4String& symbol);

  // Operation: GetNuclearMass
  //      Nuclear_Mass = Atomic_Mass - electronMass
  static G4double GetNuclearMass(const G4String& symbol);

  // Operation: GetBindingEnergy
  static G4double GetBindingEnergy(const G4String& symbol);

  // Operation: GetBetaDecayEnergy
  static G4double GetBetaDecayEnergy(const G4String& symbol);

  // Is the nucleus (A,Z) in table?
  static G4bool IsInTable(const G4String& symbol);

	

private:

  // Operation: GetIndex
  static G4int GetIndex(G4int Z, G4int A);
  

  // Data Members for Class Attributes
  //----------------------------------  

  // The following arrays are static for allow inicialization.
  // The inicialization is Done in NuclearMasses.cc

  // Mass Excess
  static const G4double MassExcess[nEntries];
  
  
  // Beta Decay Energy
  static const G4double BetaEnergy[nEntries];

    
  // Table of Z (number of protons) and A (number of nucleons)
  //        indexArray[0][ ] --> Z
  //        indexArray[1][ ] --> A
  static const G4int indexArray[2][nEntries];

  // Reduced Table of A for shorter index search.
  //         The index in this table coincide with A-1
  //         For each A value shortTable[A-1] has the index of the 1st occurrence in
  //         the indexArray[][]
  static const G4int shortTable[MaxA+1];

  // electrom mass
  static G4ThreadLocal G4double electronMass[ZMax];
  static G4ThreadLocal G4bool   isIntialized;

};

}

#endif






