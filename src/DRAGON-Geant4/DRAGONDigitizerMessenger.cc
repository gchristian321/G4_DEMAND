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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DRAGONDigitizerMessenger.hh"
#include "DRAGONDigitizer.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace DRAGON {

DRAGONDigitizerMessenger::DRAGONDigitizerMessenger(DRAGONDigitizer* DRAGONdigitizer)
:fDRAGONDigitizer(DRAGONdigitizer)
{
  fDigiDir = new G4UIdirectory("/DRAGON/digi/");
  fDigiDir->SetGuidance("Control the digitalization step");
          
  fANALCmd = new G4UIcmdWithADouble("/DRAGON/digi/anal", this);  
  fANALCmd->SetGuidance("Set the E_threshold parameter");
  fANALCmd->SetParameterName("E_threshold", false);
  fANALCmd->SetDefaultValue(2.0);
  fANALCmd->SetRange("E_threshold>=0");
  fANALCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTHLDCmd = new G4UIcmdWithAString("/DRAGON/digi/thld", this);  
  fTHLDCmd->SetGuidance("Set the tot_thrshld and pmt_thrshld parameters");
  fTHLDCmd->SetParameterName("thld", false);
  fTHLDCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONDigitizerMessenger::~DRAGONDigitizerMessenger()
{
  delete fDigiDir;
  delete fANALCmd;
  delete fTHLDCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  DRAGONDigitizerMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fANALCmd) 
      {fDRAGONDigitizer->SetANAL(fANALCmd->GetNewDoubleValue(newValue));}
  else if(command == fTHLDCmd)
      fDRAGONDigitizer->SetTHLD(newValue); 
  else 
     G4cerr << "Error reading user command: " << newValue << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
