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
/// \file SAXSPhysicsListMessenger.cc
/// \brief Definition of the SAXSPhysicsListMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DRAGONEMFieldMessenger.hh"
#include "DRAGONEMField.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace DRAGON {

DRAGONEMFieldMessenger::DRAGONEMFieldMessenger(DRAGONEMField* field)
:fEMField(field)
{
  fEMFieldDir = new G4UIdirectory("/DRAGON/emf/");
  fEMFieldDir->SetGuidance("Control the electromagnetic field");
          
  //OJO esto hay que eliminarlo (no es un parametro que deba definir el usuario)
  fSCALCmd = new G4UIcmdWithAString("/DRAGON/emf/scal", this); 
  fSCALCmd->SetParameterName("scal", false); 
  fSCALCmd->SetGuidance("Set the bscale and escale parameters");
  fSCALCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStepperCMD = new G4UIcmdWithAnInteger("/DRAGON/emf/setStepperType", this);
  fStepperCMD->SetGuidance("Select stepper type for field");
  fStepperCMD->SetParameterName("choice", true);
  fStepperCMD->SetDefaultValue(4);
  fStepperCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMinStepCMD = new G4UIcmdWithADouble("/DRAGON/emf/setMinStep", this);
  fMinStepCMD->SetGuidance("Define minimal step");
  fMinStepCMD->SetParameterName("min step", true);
  fMinStepCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDeltaChordCMD = new G4UIcmdWithADouble("/DRAGON/emf/setDeltaChord", this);
  fDeltaChordCMD->SetGuidance("Define delta chord");
  fDeltaChordCMD->SetParameterName("delta chord", true);
  fDeltaChordCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDeltaOneStepCMD = new G4UIcmdWithADouble("/DRAGON/emf/setDeltaOneStep", this);
  fDeltaOneStepCMD->SetGuidance("Define delta one step");
  fDeltaOneStepCMD->SetParameterName("delta one step", true);
  fDeltaOneStepCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDeltaIntersectionCMD = new G4UIcmdWithADouble("/DRAGON/emf/setDeltaIntersection", this);
  fDeltaIntersectionCMD->SetGuidance("Define delta intersection");
  fDeltaIntersectionCMD->SetParameterName("delta intersection", true);
  fDeltaIntersectionCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEpsMinCMD = new G4UIcmdWithADouble("/DRAGON/emf/setEpsMin", this);
  fEpsMinCMD->SetGuidance("Define eps min");
  fEpsMinCMD->SetParameterName("eps min", true);
  fEpsMinCMD->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEpsMaxCMD = new G4UIcmdWithADouble("/DRAGON/emf/setEpsMax", this);
  fEpsMaxCMD->SetGuidance("Define eps max");
  fEpsMaxCMD->SetParameterName("eps max", true);
  fEpsMaxCMD->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONEMFieldMessenger::~DRAGONEMFieldMessenger()
{
  delete fEMFieldDir;
  delete fSCALCmd;
  delete fStepperCMD;
  delete fMinStepCMD;
  delete fDeltaChordCMD;
  delete fDeltaOneStepCMD;
  delete fDeltaIntersectionCMD;
  delete fEpsMinCMD;
  delete fEpsMaxCMD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  DRAGONEMFieldMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fSCALCmd) 
      {fEMField->SetSCAL(newValue);}     
  else if (command == fStepperCMD) {
    fEMField->SetStepperType(fStepperCMD->GetNewIntValue(newValue));
  }
  else if (command == fMinStepCMD) {
    fEMField->SetMinStep(fMinStepCMD->GetNewDoubleValue(newValue));
  }
  else if (command == fDeltaChordCMD) {
    fEMField->SetDeltaChord(fDeltaChordCMD->GetNewDoubleValue(newValue));
  }
  else if (command == fDeltaOneStepCMD) {
    fEMField->SetDeltaOneStep(fDeltaOneStepCMD->GetNewDoubleValue(newValue));
  }
 else  if (command == fDeltaIntersectionCMD) {
    fEMField->SetDeltaIntersection(fDeltaIntersectionCMD->GetNewDoubleValue(newValue));
  }
  else if (command == fEpsMinCMD) {
    fEMField->SetEpsMin(fEpsMinCMD->GetNewDoubleValue(newValue));
  }
  else if (command == fEpsMaxCMD) {
    fEMField->SetEpsMax(fEpsMaxCMD->GetNewDoubleValue(newValue));
  }
  else 
      G4cerr << "Error reading user command: " << newValue << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
