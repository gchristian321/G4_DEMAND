#include <sstream>
#include <algorithm>
#include <iterator>

#include "DemandPrimaryGeneratorMessenger.hh"
#include "DemandPrimaryGeneratorAction.hh"

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
using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandPrimaryGeneratorMessenger::DemandPrimaryGeneratorMessenger(DemandPrimaryGeneratorAction* primary):
	fPrimary(primary)
{
  //Setup a command directory for detector controls with guidance
  fPrimaryDir = new G4UIdirectory("/pterp/generator/");
  fPrimaryDir->SetGuidance("Capture reaction parameters control");
 
  //Various commands for modifying primary generator
	fBeamZ = new G4UIcmdWithAnInteger("/pterp/generator/beam/Z",this);
  fBeamZ->SetGuidance("Beam atomic number (Z).");
  fBeamZ->SetParameterName("Z_beam",false);
  fBeamZ->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamZ->SetToBeBroadcasted(false);
	fBeamZ->SetDefaultValue(0);
	
	fBeamA = new G4UIcmdWithAnInteger("/pterp/generator/beam/A",this);
  fBeamA->SetGuidance("Beam atomic mass (A).");
  fBeamA->SetParameterName("A_beam",false);
  fBeamA->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamA->SetToBeBroadcasted(false);
	fBeamA->SetDefaultValue(1);

	fSourceZ = new G4UIcmdWithAnInteger("/pterp/generator/source/Z",this);
  fSourceZ->SetGuidance("Source atomic number (Z).");
  fSourceZ->SetParameterName("Z_source",false);
  fSourceZ->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSourceZ->SetToBeBroadcasted(false);
	fSourceZ->SetDefaultValue(0);
	
	fSourceA = new G4UIcmdWithAnInteger("/pterp/generator/source/A",this);
  fSourceA->SetGuidance("Source atomic mass (A).");
  fSourceA->SetParameterName("A_source",false);
  fSourceA->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSourceA->SetToBeBroadcasted(false);
	fSourceA->SetDefaultValue(1);

  fBeamEnergy =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/beam/energy",this);
  fBeamEnergy->SetGuidance("Set the beam energy");
  fBeamEnergy->SetParameterName("e_beam",false);
  fBeamEnergy->SetDefaultUnit("MeV");
  fBeamEnergy->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamEnergy->SetToBeBroadcasted(false);
	fBeamEnergy->SetDefaultValue(1.*MeV);

	fReactionCmd =
		new G4UIcmdWithAString("/pterp/generator/reaction", this);
	fReactionCmd->SetGuidance("Set the reaction");
	fReactionCmd->SetParameterName("reaction",false);
	fReactionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fReactionCmd->SetToBeBroadcasted(false);
	fReactionCmd->SetDefaultValue("none");

	fAngularDistributionCmd =
		new G4UIcmdWithAString("/pterp/generator/reaction/angdist", this);
	fAngularDistributionCmd->SetGuidance(
		"Set the reaction angular distribution");
	fAngularDistributionCmd->SetParameterName("ang_dist",false);
	fAngularDistributionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fAngularDistributionCmd->SetToBeBroadcasted(false);
	fAngularDistributionCmd->SetDefaultValue("flat");

  fExRecoil =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/reaction/excitation",this);
  fExRecoil->SetGuidance("Set the recoil excitation energy");
  fExRecoil->SetParameterName("ex_recoil",false);
  fExRecoil->SetDefaultUnit("MeV");
  fExRecoil->AvailableForStates(G4State_PreInit,G4State_Idle);
  fExRecoil->SetToBeBroadcasted(false);
	fExRecoil->SetDefaultValue(0);

  fBeamSigmaX =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/beam/sigma_x",this);
  fBeamSigmaX->SetGuidance("Set the beam x-position spread");
  fBeamSigmaX->SetParameterName("sig_x",false);
  fBeamSigmaX->SetDefaultUnit("mm");
  fBeamSigmaX->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamSigmaX->SetToBeBroadcasted(false);
	fBeamSigmaX->SetDefaultValue(0);

  fBeamSigmaY =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/beam/sigma_y",this);
  fBeamSigmaY->SetGuidance("Set the beam y-position spread");
  fBeamSigmaY->SetParameterName("sig_y",false);
  fBeamSigmaY->SetDefaultUnit("mm");
  fBeamSigmaY->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamSigmaY->SetToBeBroadcasted(false);
	fBeamSigmaY->SetDefaultValue(0);

  fBeamSigmaThetaX =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/beam/sigma_thx",this);
  fBeamSigmaThetaX->SetGuidance("Set the beam x-angle spread");
  fBeamSigmaThetaX->SetParameterName("sig_thx",false);
  fBeamSigmaThetaX->SetDefaultUnit("mrad");
  fBeamSigmaThetaX->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamSigmaThetaX->SetToBeBroadcasted(false);
	fBeamSigmaThetaX->SetDefaultValue(0);
	
  fBeamSigmaThetaY =
    new G4UIcmdWithADoubleAndUnit("/pterp/generator/beam/sigma_thy",this);
  fBeamSigmaThetaY->SetGuidance("Set the beam y-angle spread");
  fBeamSigmaThetaY->SetParameterName("sig_thy",false);
  fBeamSigmaThetaY->SetDefaultUnit("mrad");
  fBeamSigmaThetaY->AvailableForStates(G4State_PreInit,G4State_Idle);
  fBeamSigmaThetaY->SetToBeBroadcasted(false);
	fBeamSigmaThetaY->SetDefaultValue(0);

	fSourceEnergies =
		new G4UIcmdWith3VectorAndUnit("/pterp/generator/source/energies",this);
  fSourceEnergies->SetGuidance("Set the range of source energies");
  fSourceEnergies->SetParameterName("source_elow","source_ehigh","dummy",false);
  fSourceEnergies->SetDefaultUnit("MeV");
  fSourceEnergies->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSourceEnergies->SetToBeBroadcasted(false);
	fSourceEnergies->SetDefaultValue(G4ThreeVector(0,10,0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandPrimaryGeneratorMessenger::~DemandPrimaryGeneratorMessenger()
{
  delete fPrimaryDir;
	delete fBeamZ;
	delete fBeamA;
	delete fSourceZ;
	delete fSourceA;
	delete fBeamEnergy;
	delete fReactionCmd;
	delete fAngularDistributionCmd;
	delete fExRecoil;
	delete fBeamSigmaX;
	delete fBeamSigmaY;
	delete fBeamSigmaThetaX;
	delete fBeamSigmaThetaY;
	delete fSourceEnergies;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(0) {
	}
  else if (command == fBeamZ || command == fSourceZ){
		G4int Z = fBeamZ->GetNewIntValue(newValue);
		fPrimary->SetBeamZ(Z);
  }
  else if (command == fBeamA || command == fSourceA){
		G4int A = fBeamA->GetNewIntValue(newValue);
		fPrimary->SetBeamA(A);
  }
	else if(command == fBeamEnergy){
    G4double ebeam = fBeamEnergy->GetNewDoubleValue(newValue);
		fPrimary->SetBeamEnergy(ebeam);
  }
	else if(command == fReactionCmd){
		fPrimary->SetReaction(newValue);
	}
	else if(command == fAngularDistributionCmd){
		fPrimary->SetReactionAngdist(newValue);
	}
	else if(command == fExRecoil){
    G4double ex = fExRecoil->GetNewDoubleValue(newValue);
		fPrimary->SetExRecoil(ex);
  }
	else if(command == fBeamSigmaX){
    G4double sig = fBeamSigmaX->GetNewDoubleValue(newValue);
		fPrimary->SetBeamSigmaX(sig);
  }
	else if(command == fBeamSigmaY){
    G4double sig = fBeamSigmaY->GetNewDoubleValue(newValue);
		fPrimary->SetBeamSigmaY(sig);
  }
	else if(command == fBeamSigmaThetaX){
    G4double sig = fBeamSigmaThetaX->GetNewDoubleValue(newValue);
		fPrimary->SetBeamSigmaThetaX(sig);
  }
	else if(command == fBeamSigmaThetaY){
    G4double sig = fBeamSigmaThetaY->GetNewDoubleValue(newValue);
		fPrimary->SetBeamSigmaThetaY(sig);
  }
	else if(command == fSourceEnergies){
		auto e = fSourceEnergies->GetNew3VectorValue(newValue);
		fPrimary->SetSourceEnergyLimits(e[0],e[1]);
	}
}
