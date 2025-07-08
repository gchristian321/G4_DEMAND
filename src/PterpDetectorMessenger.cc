#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4Scintillation.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "PterpDetectorMessenger.hh"
#include "PterpDetectorConstruction.hh"
#include "PterpPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpDetectorMessenger::PterpDetectorMessenger(
	PterpDetectorConstruction* detector)
	: fPterpDetector(detector)
{
  //Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/pterp/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");
 
  //Various commands for modifying detector geometry
  fDimensionsCmdX =
    new G4UIcmdWithADoubleAndUnit("/pterp/detector/wx",this);
  fDimensionsCmdX->SetGuidance("Set the x-dimensions of the detector volume.");
  fDimensionsCmdX->SetParameterName("wx",false);
  fDimensionsCmdX->SetDefaultUnit("cm");
  fDimensionsCmdX->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDimensionsCmdX->SetToBeBroadcasted(false);
	fDimensionsCmdX->SetDefaultValue(2.54*cm);

	fDimensionsCmdY =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/wy",this);
  fDimensionsCmdY->SetGuidance("Set the y-dimensions of the detector volume.");
  fDimensionsCmdY->SetParameterName("wy",false);
  fDimensionsCmdY->SetDefaultUnit("cm");
  fDimensionsCmdY->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDimensionsCmdY->SetToBeBroadcasted(false);
	fDimensionsCmdY->SetDefaultValue(2.54*cm);
	
  fDimensionsCmdZ =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/wz",this);
  fDimensionsCmdZ->SetGuidance("Set the z-dimensions of the detector volume.");
  fDimensionsCmdZ->SetParameterName("wz",false);
  fDimensionsCmdZ->SetDefaultUnit("cm");
  fDimensionsCmdZ->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDimensionsCmdZ->SetToBeBroadcasted(false);
	fDimensionsCmdZ->SetDefaultValue(2.54*cm);

  fPosCmdX =
    new G4UIcmdWithADoubleAndUnit("/pterp/detector/xpos",this);
  fPosCmdX->SetGuidance("Set the x-position of the detector volume.");
  fPosCmdX->SetParameterName("px",false);
  fPosCmdX->SetDefaultUnit("cm");
  fPosCmdX->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPosCmdX->SetToBeBroadcasted(false);
	fPosCmdX->SetDefaultValue(0*cm);

	fPosCmdY =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/ypos",this);
  fPosCmdY->SetGuidance("Set the y-position of the detector volume.");
  fPosCmdY->SetParameterName("py",false);
  fPosCmdY->SetDefaultUnit("cm");
  fPosCmdY->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPosCmdY->SetToBeBroadcasted(false);
	fPosCmdY->SetDefaultValue(0*cm);
	
  fPosCmdZ =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/zpos",this);
  fPosCmdZ->SetGuidance("Set the z-position of the detector volume.");
  fPosCmdZ->SetParameterName("pz",false);
  fPosCmdZ->SetDefaultUnit("cm");
  fPosCmdZ->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPosCmdZ->SetToBeBroadcasted(false);
	fPosCmdZ->SetDefaultValue(0*cm);
	
  fDistanceCmd =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/dist",this);
  fDistanceCmd->SetGuidance(
		"Set the detector distance to front face");
  fDistanceCmd->SetParameterName("dist",false);
  fDistanceCmd->SetDefaultUnit("cm");
  fDistanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDistanceCmd->SetToBeBroadcasted(false);
	fDistanceCmd->SetDefaultValue(100*cm);

	
	fThetaXCmd =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/angle",this);
  fThetaXCmd->SetGuidance(
		"Set the angle of the detector relative to beam.");
  fThetaXCmd->SetParameterName("angle",false);
  fThetaXCmd->SetDefaultUnit("deg");
  fThetaXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fThetaXCmd->SetToBeBroadcasted(false);
	fThetaXCmd->SetDefaultValue(0*deg);
	
	fThetaCmd =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/theta",this);
  fThetaCmd->SetGuidance(
		"Set the theta angle of the detector relative to beam.");
  fThetaCmd->SetParameterName("theta",false);
  fThetaCmd->SetDefaultUnit("deg");
  fThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fThetaCmd->SetToBeBroadcasted(false);
	fThetaCmd->SetDefaultValue(0*deg);

	fPhiCmd =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/phi",this);
  fPhiCmd->SetGuidance(
		"Set the phi angle of the detector relative to beam.");
  fPhiCmd->SetParameterName("phi",false);
  fPhiCmd->SetDefaultUnit("deg");
  fPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPhiCmd->SetToBeBroadcasted(false);
	fPhiCmd->SetDefaultValue(0*deg);


	fThresholdCmd =
		new G4UIcmdWithADoubleAndUnit("/pterp/detector/thresh",this);
  fThresholdCmd->SetGuidance(
		"Set the detector threshold");
  fThresholdCmd->SetParameterName("thresh",false);
  fThresholdCmd->SetDefaultUnit("keV");
  fThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fThresholdCmd->SetToBeBroadcasted(false);
	fThresholdCmd->SetDefaultValue(50*keV);

	
  fNxCmd = new G4UIcmdWithAnInteger("/pterp/detector/nx",this);
  fNxCmd->SetGuidance("Set the number of detectors along the x-dimension.");
  fNxCmd->SetParameterName("nx",false);
  fNxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNxCmd->SetToBeBroadcasted(false);
	fNxCmd->SetDefaultValue(3);

  fNyCmd = new G4UIcmdWithAnInteger("/pterp/detector/ny",this);
  fNyCmd->SetGuidance("Set the number of detectors along the y-dimension.");
  fNyCmd->SetParameterName("ny",false);
  fNyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNyCmd->SetToBeBroadcasted(false);
	fNyCmd->SetDefaultValue(1);

  fNzCmd = new G4UIcmdWithAnInteger("/pterp/detector/nz",this);
  fNzCmd->SetGuidance("Set the number of detectors along the z-dimension.");
  fNzCmd->SetParameterName("nz",false);
  fNzCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNzCmd->SetToBeBroadcasted(false);
	fNzCmd->SetDefaultValue(9);

	fReadoutTypeCmd =
		new G4UIcmdWithAString("/pterp/detector/readout", this);
	fReadoutTypeCmd->SetGuidance("Set the readout method (cube or bar)");
	fReadoutTypeCmd->SetParameterName("readout_type",false);
	fReadoutTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fReadoutTypeCmd->SetToBeBroadcasted(false);
	fReadoutTypeCmd->SetDefaultValue("cube");

	fRotateCmd =
		new G4UIcmdWithABool("/pterp/detector/rotate", this);
	fRotateCmd->SetGuidance("Tell whether to rotate detector so that it faces the target (true) or in beam direction (false)");
	fRotateCmd->SetParameterName("do_rotate",false);
	fRotateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fRotateCmd->SetToBeBroadcasted(false);
	fRotateCmd->SetDefaultValue(true);

  fTargetThicknessCmd =
    new G4UIcmdWithADoubleAndUnit(
			"/pterp/detector/target/thickness",this);
  fTargetThicknessCmd->SetGuidance("Set the target thickness");
  fTargetThicknessCmd->SetParameterName("trgt_thick",false);
  fTargetThicknessCmd->SetDefaultUnit("micrometer");
  fTargetThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTargetThicknessCmd->SetToBeBroadcasted(false);
	fTargetThicknessCmd->SetDefaultValue(0);
	
	fTargetMaterialCmd =
		new G4UIcmdWithAString("/pterp/detector/target/material", this);
	fTargetMaterialCmd->SetGuidance("Set the target material");
	fTargetMaterialCmd->SetParameterName("target_material",false);
	fTargetMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	fTargetMaterialCmd->SetToBeBroadcasted(false);
	fTargetMaterialCmd->SetDefaultValue("CD2");

  fAddModuleCmd = new G4UIcmdWithoutParameter("/pterp/detector/module",this);
  fAddModuleCmd->SetGuidance("Add additional detector module.");
//  fAddModuleCmd->SetParameterName("add_module",false);
  fAddModuleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAddModuleCmd->SetToBeBroadcasted(false);

	fCurrentModuleNumber = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpDetectorMessenger::~PterpDetectorMessenger()
{
  delete fDimensionsCmdX;
  delete fDimensionsCmdY;
  delete fDimensionsCmdZ;
  delete fPosCmdX;
  delete fPosCmdY;
  delete fPosCmdZ;
	delete fDistanceCmd;
	delete fThetaXCmd;
	delete fThetaCmd;
	delete fPhiCmd;
	delete fThresholdCmd;
  delete fNxCmd;
  delete fNyCmd;
  delete fNzCmd;
	delete fReadoutTypeCmd;
	delete fRotateCmd;
	delete fTargetThicknessCmd;
	delete fTargetMaterialCmd;
	delete fAddModuleCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PterpDetectorMessenger::SetNewValue(
	G4UIcommand* command, G4String newValue)
{
	auto module = fPterpDetector->
		GetModuleByOriginalSpecification(fCurrentModuleNumber);
	if(!module) {
		throw std::logic_error(
			"PterpDetectorMessenger::SetNewValue -- "
			"module does not exist");
	}
  if(0) {
	}
	else if(command == fDimensionsCmdX){
    module->m_Voxsize[0] =
			fDimensionsCmdX->GetNewDoubleValue(newValue);
  }
	else if(command == fDimensionsCmdY){
    module->m_Voxsize[1] =
			fDimensionsCmdY->GetNewDoubleValue(newValue);
  }
	else if(command == fDimensionsCmdZ){
    module->m_Voxsize[2] =
			fDimensionsCmdZ->GetNewDoubleValue(newValue);
  }
	else if(command == fPosCmdX){
    module->m_Position[0] =
			fPosCmdX->GetNewDoubleValue(newValue);
		module->m_HavePosition[0] = 1;
  }
	else if(command == fPosCmdY){
    module->m_Position[1] =
			fPosCmdY->GetNewDoubleValue(newValue);
		module->m_HavePosition[1] = 1;
  }
	else if(command == fPosCmdZ){
    module->m_Position[2] =
			fPosCmdZ->GetNewDoubleValue(newValue);
		module->m_HavePosition[2] = 1;
  }
  else if (command == fNxCmd){
		module->m_Voxnum[0] =
 			fNxCmd->GetNewIntValue(newValue);
  }
  else if (command == fNyCmd){
		module->m_Voxnum[1] =
			fNyCmd->GetNewIntValue(newValue);
  }
  else if (command == fNzCmd){
		module->m_Voxnum[2] =
			fNzCmd->GetNewIntValue(newValue);
  }
	else if(command == fDistanceCmd){
    module->m_Dist = 
			fDistanceCmd->GetNewDoubleValue(newValue);
	}
	else if(command == fThetaXCmd){
		module->m_Theta =
			fThetaXCmd->GetNewDoubleValue(newValue);
		module->m_Phi = 0;
	}
	else if(command == fThetaCmd){
		module->m_Theta =
			fThetaCmd->GetNewDoubleValue(newValue);
	}
	else if(command == fPhiCmd){
		module->m_Phi =
			fPhiCmd->GetNewDoubleValue(newValue);
	}
	else if(command == fThresholdCmd){
		module->m_Threshold =
			fThresholdCmd->GetNewDoubleValue(newValue);
	}
	else if(command == fRotateCmd){
		module->m_Rotate =
			fRotateCmd->GetNewBoolValue(newValue);
	}
	else if(command == fReadoutTypeCmd){
		if(newValue == "cube") {
			module->m_ReadoutType = PterpDetectorConstruction::kCube;
		}
		else if(newValue == "bar") {
			module->m_ReadoutType = PterpDetectorConstruction::kBar;
		}
		else {
			char buf[4096];
			sprintf(buf,"PterpDetectorMessenger: invalid argument for readout type: \"%s\""
							". Valid entries are \"bar\" or \"cube\".", newValue.c_str());
			throw std::invalid_argument(buf);
		}
	}
	else if(command == fTargetMaterialCmd){
		G4String mat = newValue;
		fPterpDetector->SetTargetMaterial(mat);
	}
	else if(command == fTargetThicknessCmd){
    G4double thick = fTargetThicknessCmd->GetNewDoubleValue(newValue);
		fPterpDetector->SetTargetThickness(thick);
  }
	else if(command == fAddModuleCmd){
		fPterpDetector->AddModule();
		++fCurrentModuleNumber;
	}
}
