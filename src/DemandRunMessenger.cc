#include <sstream>
#include <algorithm>
#include <iterator>

#include "DemandRunMessenger.hh"
#include "DemandRunAction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandRunMessenger::DemandRunMessenger(DemandRunAction* run):
	fRun(run)
{
  //Setup a command directory for detector controls with guidance
	fRunDir = new G4UIdirectory("/demand/run/");
	fRunDir->SetGuidance("Run actions");
 
	fGeant3SetupCmd =
		new G4UIcmdWithAString("/demand/run/geant3", this);
	fGeant3SetupCmd->SetGuidance(
		"Set up GEANT3 input");
	fGeant3SetupCmd->SetParameterName("g3fname",false);
	fGeant3SetupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fGeant3SetupCmd->SetToBeBroadcasted(false);
	fGeant3SetupCmd->SetDefaultValue("");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandRunMessenger::~DemandRunMessenger()
{
//  delete fRunDir;
	delete fGeant3SetupCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandRunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fGeant3SetupCmd){
		G4long nevt = fRun->SetupGeant3Input(newValue);
		G4RunManager::GetRunManager()->BeamOn(nevt);
	}
}
