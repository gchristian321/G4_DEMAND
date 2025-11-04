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
/// \file hadronic/Hadr01/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
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
// 

#include "DRAGONPhysicsList.hh"
#include "DRAGONPhysicsListMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4RunManager.hh"
#include "DRAGONHistoManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace DRAGON
{

DRAGONPhysicsListMessenger::DRAGONPhysicsListMessenger(DRAGONPhysicsList* Phys)
  :G4UImessenger(),fPhys(Phys)
{
//C.
//C *** Full Monte Simulation
//C.
//C *** Define the reaction number(I), beam charge(real), recoil charge(real)
//C *** Using a negative value of the reaction number loads that reaction and
//C *** tune, but passes the beam rather than recoils.
	fFKINCmd = new G4UIcmdWithAString("/DRAGON/phys/fkin", this);
	fFKINCmd->SetGuidance("Set lkine, fkine[0] and fkine[1]");
	fFKINCmd->AvailableForStates(G4State_Idle,G4State_PreInit);  
	
	fMINPUTCmd = new G4UIcmdWithAString("/DRAGON/phys/input", this);
	fMINPUTCmd->SetGuidance("Set input namelist for the reaction");
	fMINPUTCmd->AvailableForStates(G4State_Idle,G4State_PreInit); 
//C.
//C *** Define the beam energy and tune scale for non-resonant reactions	
    fBEAMCmd = new G4UIcmdWithADouble("/DRAGON/phys/beam", this);
	fBEAMCmd->SetGuidance("Set beam energy (MeV)");
	fBEAMCmd->SetParameterName("beamenerg", false);
	fBEAMCmd->SetDefaultValue(0);
	fBEAMCmd->SetRange("beamenerg>0");
	fBEAMCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
//C.
//C*** Define the reference tune specifying reference energy(MeV), atomic#,charge
        fTUNECmd = new G4UIcmdWithAString("/DRAGON/phys/tune", this);
        fTUNECmd->SetParameterName("tune", true);
	fTUNECmd->SetGuidance("Set refenerg[3] parameters");
	fTUNECmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
        fMTUNCmd = new G4UIcmdWithAString("/DRAGON/phys/mtun", this);
	fMTUNCmd->SetGuidance("Set offset[5] parameters (cm)");
        fMTUNCmd->SetParameterName("mtun", false);
	fMTUNCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fMCUTCmd = new G4UIcmdWithAString("/DRAGON/phys/cuts", this);
	fMCUTCmd->SetGuidance("Set secondaries production cut values (cutgam, cutele, cutmuo, cuthad and cutneu) in mm");
        fMCUTCmd->SetParameterName("cuts", false);
	fMCUTCmd->AvailableForStates(G4State_Idle,G4State_PreInit); 	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONPhysicsListMessenger::~DRAGONPhysicsListMessenger()
{
delete fFKINCmd;
delete fMINPUTCmd;
delete fBEAMCmd;
delete fTUNECmd;
delete fMTUNCmd;
delete fMCUTCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fBEAMCmd)
    	{
		fPhys->SetBEAM(fBEAMCmd->GetNewDoubleValue(newValue));
           }
 else if(command ==  fFKINCmd) {
		fPhys->SetFKIN(newValue);
//C     .                                          //From uginit.f
//C     .-->   Define radioactive ion reactions
//C     .

fPhys->ureact();   
G4RunManager::GetRunManager()->PhysicsHasBeenModified();
fPhys->beaminit();
//fPhys->GetDRAGONRunAction()->GetHistoManager()->uhinit();

fPhys->AddIonGasModels();   
}
else if(command ==  fMINPUTCmd)
{
fPhys->SetInputCardname(newValue);
}
	else if(command ==  fTUNECmd) {
		fPhys->SetTUNE(newValue);
	}
else if(command ==  fMTUNCmd) {
		fPhys->SetMTUN(newValue);
	}
else if (command == fMCUTCmd)
        {
         fPhys->SetCUTs(newValue);
         G4RunManager::GetRunManager()->PhysicsHasBeenModified();
         }
else 
	   G4cerr << "Error reading user command: " << newValue << G4endl;
}
	

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
