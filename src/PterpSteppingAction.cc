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
// $Id$
//
/// \file PterpSteppingAction.cc
/// \brief Implementation of the PterpSteppingAction class

#include "PterpSteppingAction.hh"
#include "PterpEventAction.hh"
#include "PterpDetectorConstruction.hh"
#include "PterpAnalysis.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpSteppingAction::PterpSteppingAction()
	: G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpSteppingAction::~PterpSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PterpSteppingAction::UserSteppingAction(const G4Step* step)
{
	if(step->GetPreStepPoint() &&
		 step->GetPreStepPoint()->GetProcessDefinedStep() &&
		 step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessType() == fTransportation)
	{
		// get volume of the current step
		G4LogicalVolume* volume 
			= step->GetPreStepPoint()->GetTouchableHandle()
			->GetVolume()->GetLogicalVolume();

		if(volume->GetName() == "logic_Pterp_LocalBox_" &&
			 step->GetTrack()->GetTrackID() == 1)
		{
			PterpAnalysis::Instance()->AddEventCrossingDetector();
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

