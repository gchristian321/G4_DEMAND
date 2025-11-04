//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr01/include/PhysicsListMessenger.hh
/// \brief Definition of the PhysicsListMessenger class
//
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsListMessenger
//
// Created: 31.01.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef DRAGONPhysicsListMessenger_h
#define DRAGONPhysicsListMessenger_h 1


#include "G4UImessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UIcmdWithAString;
class G4UIcmdWithADouble;

namespace DRAGON
{
class DRAGONPhysicsList;

class DRAGONPhysicsListMessenger: public G4UImessenger
{
public:
  DRAGONPhysicsListMessenger(DRAGONPhysicsList*);
  ~DRAGONPhysicsListMessenger() override;
    
  void SetNewValue(G4UIcommand*, G4String) override;
  
private:
  
  DRAGONPhysicsList* fPhys;
   
  G4UIcmdWithAString* fFKINCmd = nullptr;
  G4UIcmdWithAString* fMINPUTCmd = nullptr;
  G4UIcmdWithAString* fTUNECmd = nullptr;
  G4UIcmdWithADouble* fBEAMCmd = nullptr;
  G4UIcmdWithAString* fMTUNCmd = nullptr;
  G4UIcmdWithAString* fMCUTCmd = nullptr;
  
};
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

