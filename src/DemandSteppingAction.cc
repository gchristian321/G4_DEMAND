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
/// \file DemandSteppingAction.cc
/// \brief Implementation of the DemandSteppingAction class

#include "DemandSteppingAction.hh"
#include "DemandEventAction.hh"
#include "DemandDetectorConstruction.hh"
#include "DemandAnalysis.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandSteppingAction::DemandSteppingAction()
	: G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandSteppingAction::~DemandSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandSteppingAction::UserSteppingAction(const G4Step* step)
{
	// Save location of all energy-depositions in all materials
	// for primary track
	if (step->GetTotalEnergyDeposit() > 0){
    auto* track = step->GetTrack();

		// check primary track
    if (track->GetParentID() == 0){
			const G4ThreeVector& pos = step->GetPostStepPoint()->GetPosition();
			auto touch = step->GetPreStepPoint()->GetTouchableHandle();
			auto pv    = touch->GetVolume();
			const G4String& vname = pv ? pv->GetName() : "NULL";

			auto pre = step->GetPreStepPoint();
			G4double time   = pre->GetGlobalTime();      // time before scatter
			G4double energy = pre->GetKineticEnergy();   // KE before scatter


			DemandAnalysis::Instance()->AddPrimaryScatter(
				pos, time, energy, vname );
		}
	}

	// Check if step intersects detector
	if(step->GetPreStepPoint() &&
		 //step->GetPreStepPoint()->GetProcessDefinedStep() &&
		 //step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessType() == fTransportation
		 step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary
		)
	{
		// get volume of the current step
		G4LogicalVolume* volume 
			= step->GetPreStepPoint()->GetTouchableHandle()
			->GetVolume()->GetLogicalVolume();

		if((volume->GetName() == "logic_Demand_LocalBox_" ||
				volume->GetName() == "DEMAND_scintLV")
			 && step->GetTrack()->GetTrackID() == 1)
		{
			DemandAnalysis::Instance()->AddEventCrossingDetector(
				step->GetPreStepPoint()->GetKineticEnergy()
				);
		}
	}
	// PrintWorldLocationOfVolume(step, "PDA1");
	// PrintWorldLocationOfVolume(step, "PDD1");
	// PrintWorldLocationOfVolume(step, "PDD2");
//	PrintWorldLocationOfVolume(step, "DEMAND_scintPV");
}

void DemandSteppingAction::PrintWorldLocationOfVolume(
	const G4Step* step, const G4String& pvName, bool abort)
{
	auto pv = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	if (!pv) return;

	if (pv->GetName() == pvName) {
		auto touch = step->GetPreStepPoint()->GetTouchableHandle();
		auto T = touch->GetHistory()->GetTopTransform().Inverse();

		G4ThreeVector world = T.NetTranslation();
		
		G4cout << "\n=== GEOMETRY PROBE ===\n"
					 << "PV: " << pv->GetName()
					 << " copy " << touch->GetCopyNumber()
					 << "\nWorld origin [mm]: " << world/CLHEP::mm << G4endl;
		auto p = step->GetPreStepPoint()->GetPosition();
		G4cout << "step pos [mm] = " << p/CLHEP::mm << G4endl;
		touch->GetVolume()
			->GetLogicalVolume()
			->GetSolid()
			->DumpInfo();
		G4cout << "=====================\n";

		if(abort){
			G4RunManager::GetRunManager()->AbortRun(true);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// === GEOMETRY PROBE ===
// PV: PDA1 copy 0
// World origin [mm]: (0,0,90.065)
// step pos [mm] = (17.4752,4.799,85.345)
// -----------------------------------------------------------
//     *** Dump for solid - PDA1 ***
//     ===================================================
//  Solid type: G4Tubs
//  Parameters:
//     inner radius : 0 mm
//     outer radius : 19.05 mm
//     half length Z: 4.72 mm
//     starting phi : 0 degrees
//     delta phi    : 360 degrees
// -----------------------------------------------------------

// === GEOMETRY PROBE ===
// PV: PDD2 copy 0
// World origin [mm]: (0,0,128.325)
// step pos [mm] = (4.19778,-1.39985,123.525)
// -----------------------------------------------------------
//     *** Dump for solid - PDD2 ***
//     ===================================================
//  Solid type: G4Tubs
//  Parameters:
//     inner radius : 0 mm
//     outer radius : 5.2 mm
//     half length Z: 4.8 mm
//     starting phi : 0 degrees
//     delta phi    : 360 degrees
// -----------------------------------------------------------

//// ---> 33.5 mm between (?)
////  --> but frame thickness is 
