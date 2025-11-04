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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DRAGONRunActionMessenger.hh"
#include "DRAGONRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


namespace DRAGON
{

DRAGONRunActionMessenger::DRAGONRunActionMessenger(DRAGONRunAction* runAction)
  : G4UImessenger(), fRunAction(runAction)
{
	fRunDir = new G4UIdirectory("/DRAGON/run/");
	fRunDir->SetGuidance("Run Action control");
	
	fREVSCmd = new G4UIcmdWithAnInteger("/DRAGON/run/revs", this);
	fREVSCmd->SetGuidance("Set irevs parameter");
	fREVSCmd->SetParameterName("irevs", true);
	fREVSCmd->SetDefaultValue(0);
	fREVSCmd->SetRange("irevs>=0");
	fREVSCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fIswitCmd = new G4UIcmdWithAString("/DRAGON/run/swit", this);
	fIswitCmd->SetGuidance("Set iswit array parameters (10)");
	fIswitCmd->SetParameterName("iswit", true);
	fIswitCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONRunActionMessenger::~DRAGONRunActionMessenger()
{
	delete fREVSCmd;
	delete fIswitCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fREVSCmd)
	  {fRunAction->SetIrevs(fREVSCmd->GetNewIntValue(newValue));}
        else if(command == fIswitCmd)
               {fRunAction->SetIswit(newValue);}
        	else 
	   G4cerr << "Error reading user command: " << newValue << G4endl;
}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
