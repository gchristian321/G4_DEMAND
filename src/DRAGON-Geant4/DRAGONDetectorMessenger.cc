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
/// \file DRAGON/src/DRAGONDetectorMessenger.cc
/// \brief Implementation of the DRAGON::DRAGONDetectorMessenger class

#include "G4UIdirectory.hh"                 //Geant4
#include "G4UIcmdWithAnInteger.hh"          
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

#include "DRAGONDetectorMessenger.hh"       //local
#include "DRAGONDetectorConstruction.hh"


namespace DRAGON
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONDetectorMessenger::DRAGONDetectorMessenger(DRAGONDetectorConstruction* det)
 :G4UImessenger(),fDetectorConstruction(det)
{
	fDirectory = new G4UIdirectory("/DRAGON/");
	fDirectory->SetGuidance("UI commands specific to DRAGON.");

	fDetDirectory = new G4UIdirectory("/DRAGON/det/");
	fDetDirectory->SetGuidance("Detector construction control");

	fTUBECmd = new G4UIcmdWithAnInteger("/DRAGON/det/tube", this);
	fTUBECmd->SetGuidance("Set tubetype parameter");
	fTUBECmd->SetParameterName("tubetype", false);
	fTUBECmd->SetDefaultValue(0);
	fTUBECmd->SetRange("tubetype>=0");
	fTUBECmd->AvailableForStates(G4State_Idle,G4State_PreInit);

	fTARGCmd = new G4UIcmdWithAnInteger("/DRAGON/det/targ", this);
	fTARGCmd->SetGuidance("Set targtype parameter");
	fTARGCmd->SetParameterName("targtype", false);
	fTARGCmd->SetDefaultValue(0);
	fTARGCmd->SetRange("targtype>=0 && targtype<=1");
	fTARGCmd->AvailableForStates(G4State_Idle,G4State_PreInit);

	fMPMTCmd = new G4UIcmdWithAnInteger("/DRAGON/det/mpmt", this);
	fMPMTCmd->SetGuidance("Set mtype_pmt parameter");
	fMPMTCmd->SetParameterName("mpmt", false);
	fMPMTCmd->SetDefaultValue(1);
	fMPMTCmd->SetRange("mpmt=1 || mpmt=2");
	fMPMTCmd->AvailableForStates(G4State_Idle,G4State_PreInit);

	fPMTRCmd = new G4UIcmdWithAString("/DRAGON/det/pmtr", this);
	fPMTRCmd->SetGuidance("Set pmt_size and pmt_length parameters (cm)");
    fPMTRCmd->SetParameterName("pmtr", true);
	fPMTRCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fWALLCmd = new G4UIcmdWithAString("/DRAGON/det/wall", this);
	fWALLCmd->SetGuidance("Set wall[3] array parameters (cm)");
	fWALLCmd->SetParameterName("wall", false);
	fWALLCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fBGAPCmd = new G4UIcmdWithADouble("/DRAGON/det/bgap", this);
	fBGAPCmd->SetGuidance("Set box_width parameter (cm)");
	fBGAPCmd->SetParameterName("bgap", false);
    fBGAPCmd->SetDefaultValue(5.08);
	fBGAPCmd->AvailableForStates(G4State_Idle,G4State_PreInit); 
	
	fFSIDCmd = new G4UIcmdWithAString("/DRAGON/det/fsid", this);
	fFSIDCmd->SetGuidance("Set s_finger, z_finger, air_gap, d_air[0], d_air[1] and d_mtl parameters (cm)");
    fFSIDCmd->SetParameterName("fsid", false);
	fFSIDCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fNDMATCmd = new G4UIcmdWithAnInteger("/DRAGON/det/dmat", this);
	fNDMATCmd->SetGuidance("Set n_detmate parameter");
	fNDMATCmd->SetParameterName("dmat", true);
        fNDMATCmd->SetDefaultValue(14);
	fNDMATCmd->SetRange("dmat>=0");
	fNDMATCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fHOLECmd = new G4UIcmdWithADouble("/DRAGON/det/hole", this);
	fHOLECmd->SetGuidance("Set aprt parameter (cm)");
	fHOLECmd->SetParameterName("hole", true);
	fHOLECmd->SetDefaultValue(4.496);
    fHOLECmd->SetRange("hole>=0");
	fHOLECmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fBLKACmd = new G4UIcmdWithADouble("/DRAGON/det/blka", this);
	fBLKACmd->SetGuidance("Set bulk_absorption parameter (%)");
	fBLKACmd->SetParameterName("bulk_absorption", true);
	fHOLECmd->SetDefaultValue(100.0);
	fBLKACmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fREFLCmd = new G4UIcmdWithADouble("/DRAGON/det/refl", this);
	fREFLCmd->SetGuidance("Set paint_absorption parameter (%)");
	fREFLCmd->SetParameterName("paint_absorption", false);
	fREFLCmd->SetDefaultValue(0.11);
	fREFLCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fOVRLCmd = new G4UIcmdWithABool("/DRAGON/det/ovrl", this);
	fOVRLCmd->SetGuidance("Set checkOverlaps parameter (true or false)");
	fOVRLCmd->SetParameterName("ovrl", false);
	fOVRLCmd->SetDefaultValue(true);
	fOVRLCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
	
	fMAXSCmd = new G4UIcmdWithAString("/DRAGON/det/maxs", this);
	fMAXSCmd->SetGuidance("Set max_step and len_max parameters (cm)");
    fMAXSCmd->SetParameterName("maxs", false);
	fMAXSCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
		
	fMASKCmd = new G4UIcmdWithAString("/DRAGON/det/mask", this);
	fMASKCmd->SetGuidance("Set MASKING file name");
    fMASKCmd->SetParameterName("mask", true);
    fMASKCmd->SetDefaultValue("rayfile.txt");
	fMASKCmd->AvailableForStates(G4State_Idle,G4State_PreInit);

    fSHLDCmd = new G4UIcmdWithAString("/DRAGON/det/shld", this);
    fSHLDCmd->SetGuidance("Set shield_end[0] and shield_end[1] parameters (cm)");
    fSHLDCmd->SetParameterName("shld", false);
	fSHLDCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONDetectorMessenger::~DRAGONDetectorMessenger()
{
	delete fMAXSCmd;
	delete fDirectory;
	delete fDetDirectory;
	delete fTUBECmd;
	delete fTARGCmd;
	delete fMPMTCmd;
	delete fPMTRCmd;
	delete fWALLCmd;
	delete fBGAPCmd;
	delete fFSIDCmd;
	delete fNDMATCmd;
	delete fHOLECmd;
	delete fBLKACmd;
	delete fREFLCmd;
	delete fOVRLCmd;
	delete fMASKCmd;
    delete fSHLDCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	if(command == fTUBECmd)
	  {fDetectorConstruction->SetTUBE(fTUBECmd->GetNewIntValue(newValue));}
	else if(command == fTARGCmd)
	       {fDetectorConstruction->SetTARG(fTARGCmd->GetNewIntValue(newValue));}
	else if(command == fMPMTCmd)
	       {fDetectorConstruction->SetMPMT(fMPMTCmd->GetNewIntValue(newValue));}
	else if(command == fPMTRCmd)
	       {fDetectorConstruction->SetPMTR(newValue);}
	else if(command == fWALLCmd)
   	       {fDetectorConstruction->SetWALL(newValue);}
	else if(command == fBGAPCmd)
	       {fDetectorConstruction->SetBGAP(fBGAPCmd->GetNewDoubleValue(newValue));}
	else if(command == fFSIDCmd)
	       {fDetectorConstruction->SetFSID(newValue);}
	else if(command == fNDMATCmd)
	       {fDetectorConstruction->SetDMAT(fNDMATCmd->GetNewIntValue(newValue));}
	else if(command == fHOLECmd)
	       {fDetectorConstruction->SetHOLE(fHOLECmd->GetNewDoubleValue(newValue));}
	else if(command == fBLKACmd)
	       {fDetectorConstruction->SetBLKA(fBLKACmd->GetNewDoubleValue(newValue));}
	else if(command == fREFLCmd)
	       {fDetectorConstruction->SetREFL(fREFLCmd->GetNewDoubleValue(newValue));}
	else if(command == fOVRLCmd) 
		   {fDetectorConstruction->SetCheckOverlaps(fOVRLCmd->GetNewBoolValue(newValue));}
	else if(command == fMAXSCmd)
	       {fDetectorConstruction->SetMaxStepsTrackLength(newValue);}
	else if(command == fMASKCmd)
	       {fDetectorConstruction->SetMASK(newValue);}
	else if(command == fSHLDCmd)
           {fDetectorConstruction->SetSHLD(newValue);}
    else
	    {G4cerr << "Error reading user command: " << newValue << G4endl;}	
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
