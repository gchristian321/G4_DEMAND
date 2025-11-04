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
//  Author: F. Poignant, floriane.poignant@gmail.com
//

#include "DRAGONPrimaryGeneratorActionMessenger.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"



namespace DRAGON
{

DRAGONPrimaryGeneratorActionMessenger::DRAGONPrimaryGeneratorActionMessenger(DRAGONPrimaryGeneratorAction* primary)
  :G4UImessenger(),fPrimary(primary)
{
//C *** Full Monte Simulation
//C.
//C.-->   Define *** USER COMMANDS ***
//C.
	fGPSDir = new G4UIdirectory("/DRAGON/gps/");
	fGPSDir->SetGuidance("GPS control");
//C.
//C *** Define the reaction number(I), beam charge(real), recoil charge(real)
//C *** Using a negative value of the reaction number loads that reaction and
//C *** tune, but passes the beam rather than recoils.
	fFKINCmd = new G4UIcmdWithAString("/DRAGON/gps/fkin", this);
	fFKINCmd->SetGuidance("Set lkine, fkine[0] and fkine[1] parameters");
        fFKINCmd->SetParameterName("fkin", false);
	fFKINCmd->AvailableForStates(G4State_Idle,G4State_PreInit);

	fKINECmd = new G4UIcmdWithAString("/DRAGON/gps/kine", this);
	fKINECmd->SetGuidance("Set ikine and pkine, based on the value of iswit[3] and iswit[4]");
        fKINECmd->SetParameterName("kine", true);
	fKINECmd->AvailableForStates(G4State_Idle,G4State_PreInit);
//C.
//C *** Gamma Detector
//C.
	fGKINCmd = new G4UIcmdWithAString("/DRAGON/gps/gkin", this);
	fGKINCmd->SetGuidance("Set mkine and gkine[10] parameters");
        fGKINCmd->SetParameterName("gkin", true);
	fGKINCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fEGAMCmd = new G4UIcmdWithAString("/DRAGON/gps/egam", this);
	fEGAMCmd->SetGuidance("Set egamma[10] parameters");
        fEGAMCmd->SetParameterName("egam", false);
	fEGAMCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fRAYFILECmd = new G4UIcmdWithAString("/DRAGON/gps/rayfile", this);
	fRAYFILECmd->SetGuidance("Set RAYFILE file name");
        fRAYFILECmd->SetParameterName("rayfile", true);
	fRAYFILECmd->AvailableForStates(G4State_Idle,G4State_PreInit);
       }

DRAGONPrimaryGeneratorActionMessenger::~DRAGONPrimaryGeneratorActionMessenger()
{
    delete fGPSDir; 
    delete fFKINCmd;
    delete fKINECmd;
    delete fGKINCmd;
    delete fEGAMCmd;
    delete fRAYFILECmd;
}

void DRAGONPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if(command ==  fKINECmd) {
		if(fPrimary->GetDRAGONRunAction()->GetIswit()[3] == 0 && fPrimary->GetDRAGONRunAction()->GetIswit()[4] == 0) {
			fPrimary->SetKINE(newValue, 0);
		}
		else if(fPrimary->GetDRAGONRunAction()->GetIswit()[3] == 0 && fPrimary->GetDRAGONRunAction()->GetIswit()[4] == 1) {
			fPrimary->SetKINE(newValue, 1);
		}
	}
	else if(command ==  fGKINCmd) {
		fPrimary->SetGKIN(newValue);
	}
	else if (command == fEGAMCmd)
	{fPrimary->SetEGAM(newValue);}
	else if (command == fRAYFILECmd)
	{fPrimary->SetRAYFILE(newValue);}
	else 
	   G4cerr << "Error reading user command: " << newValue << G4endl;
}

}
